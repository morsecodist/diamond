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
#include "../score_vector_int16.h"
#include "../score_vector_int8.h"
#include "../../util/log_stream.h"
#include "../../util/dynamic_iterator.h"
#include "../../data/sequence_set.h"

using std::list;
using std::atomic;
using std::thread;
using std::array;

namespace DP { namespace Swipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(const Sequence& query, Frame frame, DynamicIterator<DpTarget>& targets, _cbs composition_bias, vector<DpTarget>& overflow, Statistics& stats);

}}}

namespace DP { namespace BandedSwipe { namespace DISPATCH_ARCH {

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(
	const Sequence&query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	_cbs composition_bias,
	vector<DpTarget> &overflow,
	Statistics &stat);

template<typename _sv, typename _traceback>
list<Hsp> swipe_dispatch_cbs(
	const Sequence&query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	const int8_t* composition_bias,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	if (composition_bias == nullptr)
		return swipe<_sv, _traceback>(query, frame, subject_begin, subject_end, NoCBS(), overflow, stat);
	else
		return swipe<_sv, _traceback>(query, frame, subject_begin, subject_end, composition_bias, overflow, stat);
}

template<typename _sv, typename _traceback>
list<Hsp> full_swipe_dispatch_cbs(
	const Sequence&query,
	Frame frame,
	DynamicIterator<DpTarget>& targets,
	const int8_t* composition_bias,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	if (composition_bias == nullptr)
		return DP::Swipe::DISPATCH_ARCH::swipe<_sv, _traceback>(query, frame, targets, NoCBS(), overflow, stat);
	else
		return DP::Swipe::DISPATCH_ARCH::swipe<_sv, _traceback>(query, frame, targets, composition_bias, overflow, stat);
}

template<typename _sv>
list<Hsp> swipe_targets(const Sequence&query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	DynamicIterator<DpTarget>* targets,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	constexpr auto CHANNELS = vector<DpTarget>::const_iterator::difference_type(::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS);
	list<Hsp> out;
	if (flags & DP::FULL_MATRIX) {
		if (flags & TRACEBACK)
			return full_swipe_dispatch_cbs<_sv, VectorTraceback>(query, frame, *targets, composition_bias, overflow, stat);
		else if (flags & WITH_COORDINATES)
			return full_swipe_dispatch_cbs<_sv, ScoreWithCoords>(query, frame, *targets, composition_bias, overflow, stat);
		else
			return full_swipe_dispatch_cbs<_sv, ScoreOnly>(query, frame, *targets, composition_bias, overflow, stat);
	}
	else {
		for (vector<DpTarget>::const_iterator i = begin; i < end; i += std::min(CHANNELS, end - i)) {
			if (flags & TRACEBACK)
				out.splice(out.end(), swipe_dispatch_cbs<_sv, VectorTraceback>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, overflow, stat));
			else
				out.splice(out.end(), swipe_dispatch_cbs<_sv, ScoreOnly>(query, frame, i, i + std::min(CHANNELS, end - i), composition_bias, overflow, stat));
		}
	}
	return out;
}

template<typename _sv>
void swipe_worker(const Sequence*query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	DynamicIterator<DpTarget>* targets,
	atomic<size_t> *next,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	list<Hsp> *out,
	vector<DpTarget> *overflow,
	Statistics *stat)
{
	Statistics stat2;
	size_t pos;
	vector<DpTarget> of;
	if (targets == nullptr) {
		while (begin + (pos = next->fetch_add(::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS)) < end)
			out->splice(out->end(), swipe_targets<_sv>(*query, begin + pos, std::min(begin + pos + ::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS, end), nullptr, frame, composition_bias, flags, of, stat2));
	}
	else
		out->splice(out->end(), swipe_targets<_sv>(*query, begin, end, targets, frame, composition_bias, flags, of, stat2));
	*overflow = std::move(of);
	*stat += stat2;
}

template<typename _sv>
list<Hsp> swipe_threads(const Sequence& query,
	vector<DpTarget>::const_iterator begin,
	vector<DpTarget>::const_iterator end,
	DynamicIterator<DpTarget>* targets,
	Frame frame,
	const int8_t *composition_bias,
	int flags,
	vector<DpTarget> &overflow,
	Statistics &stat) {
	if (end <= begin && (!targets || targets->count == 0))
		return {};

	std::unique_ptr<DynamicIterator<DpTarget>> my_targets;
	if (targets == nullptr && (flags & DP::FULL_MATRIX))
		my_targets.reset(new VectorIterator<DpTarget>(begin, end));

	if (flags & PARALLEL) {
		task_timer timer("Banded swipe (run)", config.target_parallel_verbosity);
		const size_t n = config.threads_align ? config.threads_align : config.threads_;
		vector<thread> threads;
		vector<list<Hsp>> thread_out(n);
		vector<vector<DpTarget>> thread_overflow(n);
		atomic<size_t> next(0);
		for (size_t i = 0; i < n; ++i)
			threads.emplace_back(
				swipe_worker<_sv>,
				&query,
				begin,
				end,
				targets ? targets : my_targets.get(),
				&next,
				frame,
				composition_bias,
				flags,
				&thread_out[i],
				&thread_overflow[i],
				&stat);
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
		return swipe_targets<_sv>(query, begin, end, targets ? targets : my_targets.get(), frame, composition_bias, flags, overflow, stat);
}

list<Hsp> recompute_reversed(const Sequence& query, Frame frame, const Bias_correction* composition_bias, Flags flags, Statistics& stat, list<Hsp>::const_iterator begin, list<Hsp>::const_iterator end) {
	array<vector<DpTarget>, 3> dp_targets;
	vector<DpTarget> overflow;
	SequenceSet reversed_targets;
	const int qlen = (int)query.length();
#ifdef __SSE4_1__
	const int min_b = qlen <= UCHAR_MAX ? 0 : 1;
#elif defined(__SSE2__)
	const int min_b = 1;
#else
	const int min_b = 2;
#endif
	for (auto i = begin; i != end; ++i)
		reversed_targets.reserve(i->target_seq.length());
	reversed_targets.finish_reserve();

	size_t j = 0;
	for (auto i = begin; i != end; ++i, ++j) {
		std::reverse_copy(i->target_seq.data(), i->target_seq.end(), reversed_targets.ptr(j));
		int b = i->score <= UCHAR_MAX ? 0 : (i->score <= USHRT_MAX ? 1 : 2);
		b = std::max(b, min_b);
		dp_targets[b].emplace_back(reversed_targets[j], i->swipe_target, i->query_range.end_, i->subject_range.end_);
	}

	list<Hsp> out;
	vector<Letter> reversed = query.reverse();
	Bias_correction rev_cbs = composition_bias ? composition_bias->reverse() : Bias_correction();
#ifdef __SSE4_1__
	out = swipe_threads<::DISPATCH_ARCH::score_vector<int8_t>>(Sequence(reversed), dp_targets[0].begin(), dp_targets[0].end(), nullptr, frame, composition_bias ? rev_cbs.int8.data() : nullptr, flags, overflow, stat);
#endif
#ifdef __SSE2__
	out.splice(out.end(), swipe_threads<::DISPATCH_ARCH::score_vector<int16_t>>(Sequence(reversed), dp_targets[1].begin(), dp_targets[1].end(), nullptr, frame, composition_bias ? rev_cbs.int8.data() : nullptr, flags, overflow, stat));
#endif
	out.splice(out.end(), swipe_threads<int32_t>(Sequence(reversed), dp_targets[2].begin(), dp_targets[2].end(), nullptr, frame, composition_bias ? rev_cbs.int8.data() : nullptr, flags, overflow, stat));
	return out;
}

template<typename It>
pair<list<Hsp>, vector<DpTarget>> swipe_bin(const unsigned bin, const Sequence &query, const It begin, const It end, Frame frame, const int8_t* composition_bias, Flags flags, Statistics &stat) {
	if (end - begin == 0)
		return;
	vector<DpTarget> overflow;
	list<Hsp> out;
	auto time_stat = flag_any(flags, Flags::TRACEBACK) ? Statistics::TIME_TRACEBACK_SW : Statistics::TIME_SW;
	if (!flag_any(flags, Flags::FULL_MATRIX))
		std::sort(begin, end);
	stat.inc(Statistics::value(Statistics::EXT8 + bin), end - begin);
	task_timer timer;
	switch (bin) {
	case 0:
		out = swipe_threads<::DISPATCH_ARCH::score_vector<int8_t>>(query, begin, end, frame, composition_bias, flags, overflow, stat);
		break;
	case 1:
		out = swipe_threads<::DISPATCH_ARCH::score_vector<int16_t>>(query, begin, end, frame, composition_bias, flags, overflow, stat);
		break;
	case 2:
		out = swipe_threads<int32_t>(query, begin, end, frame, composition_bias, flags, overflow, stat);
		break;
	default:
		throw std::runtime_error("Invalid SWIPE bin.");
	}
	if (!flag_any(flags, Flags::PARALLEL)) stat.inc(time_stat, timer.microseconds());
	return { out, overflow };
}

list<Hsp> swipe(const Sequence &query, const array<vector<DpTarget>, BINS> &targets, const Frame frame, const Bias_correction *composition_bias, const DP::Flags flags, Statistics &stat)
{
	const auto cbs = composition_bias ? composition_bias->int8.data() : nullptr;
	pair<list<Hsp>, vector<DpTarget>> result;
	list<Hsp> out;
	for (unsigned bin = 0; bin < BINS; ++bin) {
		vector<DpTarget> round_targets;
		round_targets.reserve(targets[bin].size() + result.second.size());
		round_targets.insert(round_targets.end(), targets[bin].begin(), targets[bin].end());
		round_targets.insert(round_targets.end(), result.second.begin(), result.second.end());
		result = swipe_bin(bin, query, round_targets.begin(), round_targets.end(), frame, cbs, flags, stat);
		out.splice(out.end(), result.first);
	}
	assert(result.second.empty());	
	return flag_any(flags, Flags::WITH_COORDINATES) ? recompute_reversed(query, frame, composition_bias, flags, stat, out.begin(), out.end()) : out;
}

list<Hsp> swipe(const Sequence& query, const array<vector<DpTarget>, BINS>& targets, const Frame frame, const Bias_correction* composition_bias, const DP::Flags flags, Statistics& stat) {
}
		
}}}
