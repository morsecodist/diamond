/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <list>
#include <atomic>
#include <thread>
#include <numeric>
#include <limits.h>
#include "../dp.h"
#include "../../util/log_stream.h"
#include "../../data/sequence_set.h"
#include "../score_vector.h"
#include "../score_vector_int16.h"
#include "../score_vector_int8.h"
#include "cell_update.h"
#include "banded_matrix.h"
#include "full_matrix.h"
#include "full_swipe.h"
#include "banded_swipe.h"

using std::list;
using std::atomic;
using std::thread;
using std::array;
using std::atomic_size_t;

template<bool tb, typename RC, typename C, typename IdM>
struct SwipeConfig {
	static constexpr bool traceback = tb;
	using RowCounter = RC;
	using Cell = C;
	using IdMask = IdM;
};

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

static void sort(const vector<DpTarget>::iterator begin, const vector<DpTarget>::iterator end) {
	std::sort(begin, end);
}

static void sort(const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end) {
}

static unsigned bin(int x) {
	return x <= UCHAR_MAX ? 0 : (x <= USHRT_MAX ? 1 : 2);
}

unsigned bin(HspValues v, int query_len, int score, int ungapped_score, size_t dp_size, unsigned score_width) {
#if !defined(__SSE4_1__) && !defined(__SSE2__)
	return 2;
#endif
	unsigned b = 0;
	if (flag_any(v, HspValues::TRANSCRIPT | HspValues::QUERY_COORDS | HspValues::IDENT | HspValues::MISMATCHES | HspValues::LENGTH))
		b = std::max(b, bin(query_len));
	if (flag_any(v, HspValues::TRANSCRIPT) && dp_size > config.max_swipe_dp)
		b = 2;
	b = std::max(b, bin(score));
	if (ungapped_score > config.cutoff_score_8bit)
		b = std::max(b, 1u);
	b = std::max(b, score_width);
#ifdef __SSE4_1__
	return b;
#else
	return std::max(b, 1u);
#endif
}

static const HspValues NO_TRACEBACK = HspValues::COORDS | HspValues::IDENT | HspValues::LENGTH | HspValues::MISMATCHES | HspValues::GAP_OPENINGS;

static bool reversed(const HspValues v) {
	return flag_only(v, NO_TRACEBACK)
		&& flag_any(v, HspValues::QUERY_START | HspValues::TARGET_START | HspValues::MISMATCHES | HspValues::GAP_OPENINGS);
}

template<typename Sv, typename Cbs, typename Cfg>
static list<Hsp> dispatch_swipe(
	const Sequence& query,
	const Frame frame,
	const vector<DpTarget>::const_iterator subject_begin,
	const vector<DpTarget>::const_iterator subject_end,
	Cbs composition_bias,
	vector<DpTarget>& overflow,
	Statistics& stat)
{
	return ::DP::BandedSwipe::DISPATCH_ARCH::swipe<Sv, Cbs, Cfg>(query, frame, subject_begin, subject_end, composition_bias, overflow, stat);
}

template<typename Sv, typename Cbs, typename Cfg>
static list<Hsp> dispatch_swipe(
	const Sequence& query,
	const Frame frame,
	const SequenceSet::ConstIterator subject_begin,
	const SequenceSet::ConstIterator subject_end,
	Cbs composition_bias,
	vector<DpTarget>& overflow,
	Statistics& stat)
{
	return {};
}

template<typename Sv, typename Cbs, typename It, typename Cfg>
static list<Hsp> dispatch_swipe(const Sequence& query,
	const It begin,
	const It end,
	atomic_size_t* const next,
	const Frame frame,
	Cbs composition_bias,
	const Flags flags,
	vector<DpTarget>& overflow,
	Statistics& stat)
{
	constexpr auto CHANNELS = vector<DpTarget>::const_iterator::difference_type(::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS);
	if (flag_any(flags, Flags::FULL_MATRIX))
		return ::DP::Swipe::DISPATCH_ARCH::swipe<Sv, Cbs, It, Cfg>(query, frame, begin, end, next, composition_bias, overflow, stat);
	else {
		list<Hsp> out;
		for (It i = begin; i < end; i += std::min(CHANNELS, end - i))
			out.splice(out.end(), dispatch_swipe<Sv, Cbs, Cfg>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, overflow, stat));
		return out;
	}
}

template<typename Sv, typename It, typename Cfg>
static list<Hsp> dispatch_swipe(
	const Sequence&query,
	const It begin,
	const It end,
	atomic_size_t* const next,
	const Frame frame,	
	const int8_t* composition_bias,
	const Flags flags,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	if (composition_bias == nullptr)
		return dispatch_swipe<Sv, NoCBS, It, Cfg>(query, begin, end, next, frame, NoCBS(), flags, overflow, stat);
	else
		return dispatch_swipe<Sv, const int8_t*, It, Cfg>(query, begin, end, next, frame, composition_bias, flags, overflow, stat);
}

template<typename Sv, typename It>
static list<Hsp> dispatch_swipe(const Sequence& query,
	const It begin,
	const It end,
	atomic_size_t* const next,
	const Frame frame,
	const int8_t *composition_bias,
	const Flags flags,
	const HspValues v,
	vector<DpTarget> &overflow,
	Statistics &stat,
	const int round)
{
	if (v == HspValues::NONE) {
		using Cfg = SwipeConfig<false, DummyRowCounter, Sv, DummyIdMask<Sv>>;
		return dispatch_swipe<Sv, It, Cfg>(query, begin, end, next, frame, composition_bias, flags, overflow, stat);
	}
	else if (!flag_only(v, NO_TRACEBACK)) {
		using Cfg = SwipeConfig<true, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
		return dispatch_swipe<Sv, It, Cfg>(query, begin, end, next, frame, composition_bias, flags, overflow, stat);
	}
	else if (round == 0) {
		if (!flag_any(v, HspValues::IDENT | HspValues::LENGTH)) {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(query, begin, end, next, frame, composition_bias, flags, overflow, stat);
		}
		else {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, ForwardCell<Sv>, VectorIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(query, begin, end, next, frame, composition_bias, flags, overflow, stat);
		}
	}
	else if (round == 1) {
		if (!flag_any(v, HspValues::MISMATCHES | HspValues::GAP_OPENINGS)) {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(query, begin, end, next, frame, composition_bias, flags, overflow, stat);
		}
		else {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, BackwardCell<Sv>, VectorIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(query, begin, end, next, frame, composition_bias, flags, overflow, stat);
		}
	}
	throw std::runtime_error("Unreachable");
}

template<typename _sv, typename It>
static void swipe_worker(const Sequence* query,
	const It begin,
	const It end,
	atomic_size_t* const next,
	const Frame frame,
	const int8_t *composition_bias,
	const Flags flags,
	const HspValues v,
	list<Hsp> *out,
	vector<DpTarget> *overflow,
	Statistics *stat,
	const int round)
{
	const size_t CHANNELS = ::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS;
	Statistics stat2;
	size_t pos;
	vector<DpTarget> of;
	if (flag_any(flags, Flags::FULL_MATRIX))
		*out = dispatch_swipe<_sv, It>(*query, begin, end, next, frame, composition_bias, flags, v, of, stat2, round);
	else
		while (begin + (pos = next->fetch_add(CHANNELS)) < end)
			out->splice(out->end(), dispatch_swipe<_sv, It>(*query, begin + pos, std::min(begin + pos + CHANNELS, end), next, frame, composition_bias, flags, v, of, stat2, round));
		
	*overflow = std::move(of);
	*stat += stat2;
}

template<typename _sv, typename It>
static list<Hsp> swipe_threads(const Sequence& query,
	const It begin,
	const It end,
	const Frame frame,
	const int8_t *composition_bias,
	const Flags flags,
	const HspValues v,
	vector<DpTarget> &overflow,
	Statistics &stat,
	const int round) {
	if (begin == end)
		return {};

	atomic_size_t next(0);
	if (flag_any(flags, Flags::PARALLEL)) {
		task_timer timer("Banded swipe (run)", config.target_parallel_verbosity);
		const size_t n = config.threads_align ? config.threads_align : config.threads_;
		vector<thread> threads;
		vector<list<Hsp>> thread_out(n);
		vector<vector<DpTarget>> thread_overflow(n);
		for (size_t i = 0; i < n; ++i)
			threads.emplace_back(
				swipe_worker<_sv, It>,
				&query,
				begin,
				end,
				&next,
				frame,
				composition_bias,
				flags,
				v,
				&thread_out[i],
				&thread_overflow[i],
				&stat,
				round);
		for (auto &t : threads)
			t.join();
		timer.go("Banded swipe (merge)");
		list<Hsp> out;
		for (list<Hsp> &l : thread_out)
			out.splice(out.end(), l);
		overflow.reserve(std::accumulate(thread_overflow.begin(), thread_overflow.end(), (size_t)0, [](size_t n, const vector<DpTarget> &v) { return n + v.size(); }));
		for (const vector<DpTarget> &v : thread_overflow)
			overflow.insert(overflow.end(), v.begin(), v.end());
		return out;
	}
	else
		return dispatch_swipe<_sv, It>(query, begin, end, &next, frame, composition_bias, flags, v, overflow, stat, round);
}

template<typename It>
static pair<list<Hsp>, vector<DpTarget>> swipe_bin(const unsigned bin, const Sequence &query, const It begin, const It end, Frame frame, const int8_t* composition_bias, const Flags flags, const HspValues v, Statistics &stat, const int round) {
	if (end - begin == 0)
		return { {},{} };
	vector<DpTarget> overflow;
	list<Hsp> out;
	auto time_stat = flag_any(v, HspValues::TRANSCRIPT) ? Statistics::TIME_TRACEBACK_SW : Statistics::TIME_SW;
	if (!flag_any(flags, Flags::FULL_MATRIX))
		sort(begin, end);
	stat.inc(Statistics::value(Statistics::EXT8 + bin), end - begin);
	task_timer timer;
	switch (bin) {
#ifdef __SSE4_1__
	case 0:
		out = swipe_threads<::DISPATCH_ARCH::score_vector<int8_t>, It>(query, begin, end, frame, composition_bias, flags, v, overflow, stat, round);
		break;
#endif
#ifdef __SSE2__
	case 1:
		out = swipe_threads<::DISPATCH_ARCH::score_vector<int16_t>, It>(query, begin, end, frame, composition_bias, flags, v, overflow, stat, round);
		break;
#endif
	case 2:
		out = swipe_threads<int32_t, It>(query, begin, end, frame, composition_bias, flags, v, overflow, stat, round);
		break;
	default:
		throw std::runtime_error("Invalid SWIPE bin.");
	}
	if (!flag_any(flags, Flags::PARALLEL)) stat.inc(time_stat, timer.microseconds());
	return { out, overflow };
}

static list<Hsp> recompute_reversed(const Sequence& query, Frame frame, const Bias_correction* composition_bias, Flags flags, HspValues v, Statistics& stat, list<Hsp>::const_iterator begin, list<Hsp>::const_iterator end) {
	Targets dp_targets;
	vector<DpTarget> overflow;
	SequenceSet reversed_targets;
	const int qlen = (int)query.length();

	for (auto i = begin; i != end; ++i)
		reversed_targets.reserve(i->target_seq.length());
	reversed_targets.finish_reserve();

	size_t j = 0;
	for (auto i = begin; i != end; ++i, ++j) {
		std::reverse_copy(i->target_seq.data(), i->target_seq.end(), reversed_targets.ptr(j));
		const int band = flag_any(flags, Flags::FULL_MATRIX) ? qlen : i->d_end - i->d_begin,
			tlen = (int)i->target_seq.length(),
			b = bin(v, band, i->score, 0, 0, 0);
		const DpTarget::CarryOver carry_over{ i->query_range.end_, i->subject_range.end_, i->identities, i->length };
		dp_targets[b].emplace_back(reversed_targets[j], -i->d_end + qlen - tlen + 1, -i->d_begin + qlen - tlen + 1, i->swipe_target, qlen, nullptr, carry_over);
	}

	list<Hsp> out;
	vector<Letter> reversed = query.reverse();
	Bias_correction rev_cbs = composition_bias ? composition_bias->reverse() : Bias_correction();
	const int8_t* cbs = composition_bias ? rev_cbs.int8.data() : nullptr;
	for (unsigned bin = 0; bin < BINS; ++bin)
		out.splice(out.end(), swipe_bin(bin, Sequence(reversed), dp_targets[bin].begin(), dp_targets[bin].end(), frame, cbs, flags, v, stat, 1).first);
	return out;
}

list<Hsp> swipe(const Sequence &query, const array<vector<DpTarget>, BINS> &targets, const Frame frame, const Bias_correction *composition_bias, const DP::Flags flags, const HspValues v, Statistics &stat)
{
	const auto cbs = composition_bias ? composition_bias->int8.data() : nullptr;
	pair<list<Hsp>, vector<DpTarget>> result;
	list<Hsp> out;
	for (unsigned bin = 0; bin < BINS; ++bin) {
		vector<DpTarget> round_targets;
		round_targets.reserve(targets[bin].size() + result.second.size());
		round_targets.insert(round_targets.end(), targets[bin].begin(), targets[bin].end());
		round_targets.insert(round_targets.end(), result.second.begin(), result.second.end());
		result = swipe_bin(bin, query, round_targets.begin(), round_targets.end(), frame, cbs, flags, v, stat, 0);
		out.splice(out.end(), result.first);
	}
	assert(result.second.empty());
	return reversed(v) ? recompute_reversed(query, frame, composition_bias, flags, v, stat, out.begin(), out.end()) : out;
}

list<Hsp> swipe_set(const Sequence& query, const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end, const Frame frame, const Bias_correction* composition_bias, const DP::Flags flags, const HspValues v, Statistics& stat) {
	const auto cbs = composition_bias ? composition_bias->int8.data() : nullptr;
	const unsigned b = bin(v, 0, 0, 0, 0, 0);
	pair<list<Hsp>, vector<DpTarget>> result = swipe_bin(b, query, begin, end, frame, cbs, flags, v, stat, 0);
	if (reversed(v))
		result.first = recompute_reversed(query, frame, composition_bias, flags, v, stat, result.first.begin(), result.first.end());
	if (b < BINS - 1 && !result.second.empty()) {
		array<vector<DpTarget>, BINS> targets;
		targets[b + 1] = std::move(result.second);
		result.first.splice(result.first.end(), swipe(query, targets, frame, composition_bias, flags, v, stat));
	}
	return result.first;
}

}}}
