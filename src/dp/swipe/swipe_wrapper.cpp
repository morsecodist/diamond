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
	return x < UCHAR_MAX ? 0 : (x < USHRT_MAX ? 1 : 2);
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

template<typename Sv>
static size_t matrix_size(const int query_len, const vector<DpTarget>::const_iterator begin, const vector<DpTarget>::const_iterator end) {
	size_t s = 0;
	for (auto i = begin; i != end; ++i)
		s = std::max(s, i->seq.length());
	return (size_t)query_len * s * ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS / 2;
}

template<typename Sv>
static size_t matrix_size(const int query_len, const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end) {
	return 0;
}

static const HspValues NO_TRACEBACK = HspValues::COORDS | HspValues::IDENT | HspValues::LENGTH | HspValues::MISMATCHES | HspValues::GAP_OPENINGS;

static bool reversed(const HspValues v) {
	return flag_only(v, NO_TRACEBACK)
		&& flag_any(v, HspValues::QUERY_START | HspValues::TARGET_START | HspValues::MISMATCHES | HspValues::GAP_OPENINGS);
}

template<typename Sv, typename Cbs, typename Cfg>
static list<Hsp> dispatch_swipe(const vector<DpTarget>::const_iterator subject_begin, const vector<DpTarget>::const_iterator subject_end, Cbs composition_bias, vector<DpTarget>& overflow, Params& p)
{
	return ::DP::BandedSwipe::DISPATCH_ARCH::swipe<Sv, Cbs, Cfg>(subject_begin, subject_end, composition_bias, overflow, p);
}

template<typename Sv, typename Cbs, typename Cfg>
static list<Hsp> dispatch_swipe(const SequenceSet::ConstIterator subject_begin, const SequenceSet::ConstIterator subject_end, Cbs composition_bias, vector<DpTarget>& overflow, Params& p)
{
	return {};
}

template<typename Sv, typename Cbs, typename It, typename Cfg>
static list<Hsp> dispatch_swipe(const It begin, const It end, atomic_size_t* const next, Cbs composition_bias, vector<DpTarget>& overflow, Params& p)
{
	constexpr auto CHANNELS = vector<DpTarget>::const_iterator::difference_type(::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS);
	if (flag_any(p.flags, Flags::FULL_MATRIX))
		return ::DP::Swipe::DISPATCH_ARCH::swipe<Sv, Cbs, It, Cfg>(begin, end, next, composition_bias, overflow, p);
	else {
		list<Hsp> out;
		for (It i = begin; i < end; i += std::min(CHANNELS, end - i))
			out.splice(out.end(), dispatch_swipe<Sv, Cbs, Cfg>(i, i + std::min(CHANNELS, end - i), composition_bias, overflow, p));
		return out;
	}
}

template<typename Sv, typename It, typename Cfg>
static list<Hsp> dispatch_swipe(const It begin, const It end, atomic_size_t* const next, vector<DpTarget> &overflow, Params& p)
{
	if (p.composition_bias == nullptr)
		return dispatch_swipe<Sv, NoCBS, It, Cfg>(begin, end, next, NoCBS(), overflow, p);
	else
		return dispatch_swipe<Sv, const int8_t*, It, Cfg>(begin, end, next, p.composition_bias, overflow, p);
}

template<typename Sv, typename It>
static list<Hsp> dispatch_swipe(const It begin, const It end, atomic_size_t* const next, vector<DpTarget> &overflow, const int round, Params& p)
{
	if (p.v == HspValues::NONE) {
		using Cfg = SwipeConfig<false, DummyRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
		return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
	}
	const size_t s = round == 0 ? matrix_size<Sv>((int)p.query.length(), begin, end) : SIZE_MAX;
	if (!flag_only(p.v, NO_TRACEBACK) || s <= config.max_traceback_matrix_size) {
		using Cfg = SwipeConfig<true, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
		return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
	}
	else if (round == 0) {
		if (!flag_any(p.v, HspValues::IDENT | HspValues::LENGTH)) {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
		else {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, ForwardCell<Sv>, VectorIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
	}
	else if (round == 1) {
		if (!flag_any(p.v, HspValues::MISMATCHES | HspValues::GAP_OPENINGS)) {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, Sv, DummyIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
		else {
			using Cfg = SwipeConfig<false, VectorRowCounter<Sv>, BackwardCell<Sv>, VectorIdMask<Sv>>;
			return dispatch_swipe<Sv, It, Cfg>(begin, end, next, overflow, p);
		}
	}
	throw std::runtime_error("Unreachable");
}

template<typename _sv, typename It>
static void swipe_worker(const It begin, const It end, atomic_size_t* const next, list<Hsp> *out, vector<DpTarget> *overflow, const int round, Params* p)
{
	const ptrdiff_t CHANNELS = ::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS;
	Statistics stat2;
	size_t pos;
	vector<DpTarget> of;
	Params params{
		p->query,
		p->frame,
		p->query_source_len,
		p->composition_bias,
		p->flags,
		p->v,
		stat2
	};
	if (flag_any(p->flags, Flags::FULL_MATRIX))
		*out = dispatch_swipe<_sv, It>(begin, end, next, of, round, params);
	else
		while (begin + (pos = next->fetch_add(CHANNELS)) < end) {
			const auto start = begin + pos;
			out->splice(out->end(), dispatch_swipe<_sv, It>(start, start + std::min(CHANNELS, end - start), next, of, round, params));
		}
		
	*overflow = std::move(of);
	p->stat += stat2;
}

template<typename _sv, typename It>
static list<Hsp> swipe_threads(const It begin, const It end, vector<DpTarget> &overflow, const int round, Params& p) {
	if (begin == end)
		return {};

	atomic_size_t next(0);
	if (flag_any(p.flags, Flags::PARALLEL)) {
		task_timer timer("Banded swipe (run)", config.target_parallel_verbosity);
		const size_t n = config.threads_align ? config.threads_align : config.threads_;
		vector<thread> threads;
		vector<list<Hsp>> thread_out(n);
		vector<vector<DpTarget>> thread_overflow(n);
		for (size_t i = 0; i < n; ++i)
			threads.emplace_back(swipe_worker<_sv, It>, begin, end, &next, &thread_out[i], &thread_overflow[i], round, &p);
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
		return dispatch_swipe<_sv, It>(begin, end, &next, overflow, round, p);
}

template<typename It>
static pair<list<Hsp>, vector<DpTarget>> swipe_bin(const unsigned bin, const It begin, const It end, const int round, Params& p) {
	if (end - begin == 0)
		return { {},{} };
	vector<DpTarget> overflow;
	list<Hsp> out;
	auto time_stat = flag_any(p.v, HspValues::TRANSCRIPT) ? Statistics::TIME_TRACEBACK_SW : Statistics::TIME_SW;
	if (!flag_any(p.flags, Flags::FULL_MATRIX))
		sort(begin, end);
	p.stat.inc(Statistics::value(Statistics::EXT8 + bin), end - begin);
	task_timer timer;
	switch (bin) {
#ifdef __SSE4_1__
	case 0:
		if (flag_any(p.flags, Flags::SEMI_GLOBAL))
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int8_t, 0>, It>(begin, end, overflow, round, p);
		else
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int8_t, SCHAR_MIN>, It>(begin, end, overflow, round, p);
		break;
#endif
#ifdef __SSE2__
	case 1:
		if (flag_any(p.flags, Flags::SEMI_GLOBAL))
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int16_t, 0>, It>(begin, end, overflow, round, p);
		else
			out = swipe_threads<::DISPATCH_ARCH::ScoreVector<int16_t, SHRT_MIN>, It>(begin, end, overflow, round, p);
		break;
#endif
	case 2:
		out = swipe_threads<int32_t, It>(begin, end, overflow, round, p);
		break;
	default:
		throw std::runtime_error("Invalid SWIPE bin.");
	}
	if (!flag_any(p.flags, Flags::PARALLEL)) p.stat.inc(time_stat, timer.microseconds());
	return { out, overflow };
}

static list<Hsp> recompute_reversed(list<Hsp> &hsps, Params& p) {
	Targets dp_targets;
	vector<DpTarget> overflow;
	SequenceSet reversed_targets;
	const int qlen = (int)p.query.length();

	for (const auto& h : hsps)
		if (!h.backtraced)
			reversed_targets.reserve(h.target_seq.length());
	reversed_targets.finish_reserve();

	size_t j = 0;
	list<Hsp> out;
	for (auto i = hsps.begin(); i != hsps.end(); ) {
		if (i->backtraced) {
			auto k = i++;
			out.splice(out.end(), hsps, k);
			continue;
		}
		std::reverse_copy(i->target_seq.data(), i->target_seq.end(), reversed_targets.ptr(j));
		const int band = flag_any(p.flags, Flags::FULL_MATRIX) ? qlen : i->d_end - i->d_begin,
			tlen = (int)i->target_seq.length(),
			b = bin(p.v, band, i->score, 0, 0, 0);
		const DpTarget::CarryOver carry_over{ i->query_range.end_, i->subject_range.end_, i->identities, i->length };
		dp_targets[b].emplace_back(reversed_targets[j], -i->d_end + qlen - tlen + 1, -i->d_begin + qlen - tlen + 1, i->swipe_target, qlen, i->matrix, carry_over);
		++i;
		++j;
	}

	vector<Letter> reversed = p.query.reverse();
	vector<int8_t> rev_cbs = Bias_correction::reverse(p.composition_bias, p.query.length());
	const int8_t* cbs = p.composition_bias ? rev_cbs.data() : nullptr;
	Params params{
		Sequence(reversed),
		p.frame,
		p.query_source_len,
		cbs,
		p.flags,
		p.v,
		p.stat
	};
	for (unsigned bin = 0; bin < BINS; ++bin)
		out.splice(out.end(), swipe_bin(bin, dp_targets[bin].begin(), dp_targets[bin].end(), 1, params).first);
	return out;
}

list<Hsp> swipe(const array<vector<DpTarget>, BINS> &targets, Params& p)
{
	pair<list<Hsp>, vector<DpTarget>> result;
	list<Hsp> out;
	for (unsigned bin = 0; bin < BINS; ++bin) {
		vector<DpTarget> round_targets;
		round_targets.reserve(targets[bin].size() + result.second.size());
		round_targets.insert(round_targets.end(), targets[bin].begin(), targets[bin].end());
		round_targets.insert(round_targets.end(), result.second.begin(), result.second.end());
		result = swipe_bin(bin, round_targets.begin(), round_targets.end(), 0, p);
		out.splice(out.end(), result.first);
	}
	assert(result.second.empty());
	return reversed(p.v) ? recompute_reversed(out, p) : out;
}

list<Hsp> swipe_set(const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end, Params& p) {
	const unsigned b = bin(p.v, 0, 0, 0, 0, 0);
	pair<list<Hsp>, vector<DpTarget>> result = swipe_bin(b, begin, end, 0, p);
	if (reversed(p.v))
		result.first = recompute_reversed(result.first, p);
	if (b < BINS - 1 && !result.second.empty()) {
		array<vector<DpTarget>, BINS> targets;
		targets[b + 1] = std::move(result.second);
		result.first.splice(result.first.end(), swipe(targets, p));
	}
	return result.first;
}

}}}
