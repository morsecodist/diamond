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

struct DummyRowCounter {
	DummyRowCounter() {}
	DummyRowCounter(int) {}
	void store(void*) {}
	template<typename _sv>
	FORCE_INLINE void inc(const _sv&, const _sv&) {}
	enum { MAX_LEN = INT_MAX };
};

template<typename _sv>
struct VectorRowCounter {
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score Score;
	VectorRowCounter(int i) :
		i(::DISPATCH_ARCH::ScoreTraits<_sv>::zero_score() + Score(i)),
		i_max()
	{
	}
	FORCE_INLINE void inc(const _sv& best, const _sv& current_cell) {
		i_max = blend(i_max, i, best == current_cell);
		i += _sv(Score(1));
	}
	void store(Score* ptr) {
		store_sv(i_max, ptr);
	}
	static constexpr int MAX_LEN = ::DISPATCH_ARCH::ScoreTraits<_sv>::max_int_score();
	_sv i;
	_sv i_max;
};

template<typename _sv>
MSC_INLINE static inline _sv add_cbs(const _sv& v, void*) {
	return v;
}

template<typename _sv>
MSC_INLINE static inline _sv add_cbs(const _sv& v, const _sv& query_bias) {
	return v + query_bias;
}

template<typename _score>
static inline _score add_cbs_scalar(_score x, int8_t b) {
	return x + _score(b);
}

template<typename _score>
static inline _score add_cbs_scalar(_score x, void* b) {
	return x;
}

template<typename _sv>
MSC_INLINE static inline void make_gap_mask(typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask* trace_mask, const _sv& current_cell, const _sv& vertical_gap, const _sv& horizontal_gap) {
	trace_mask->gap = ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask::make(cmp_mask(current_cell, vertical_gap), cmp_mask(current_cell, horizontal_gap));
}

template<typename _sv>
MSC_INLINE static inline void make_gap_mask(std::nullptr_t, const _sv&, const _sv&, const _sv&) {
}

template<typename _sv>
MSC_INLINE static inline void make_open_mask(typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask* trace_mask, const _sv& open, const _sv& vertical_gap, const _sv& horizontal_gap) {
	trace_mask->open = ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask::make(cmp_mask(vertical_gap, open), cmp_mask(horizontal_gap, open));
}

template<typename _sv>
MSC_INLINE static inline void make_open_mask(std::nullptr_t, const _sv&, const _sv&, const _sv&) {
}

template<typename _sv, typename _cbs, typename _trace_mask, typename RowCounter>
MSC_INLINE static inline _sv swipe_cell_update(const _sv& diagonal_cell,
	const _sv& scores,
	_cbs query_bias,
	const _sv& gap_extension,
	const _sv& gap_open,
	_sv& horizontal_gap,
	_sv& vertical_gap,
	_sv& best,
	_trace_mask trace_mask,
	RowCounter& row_counter)
{
	using std::max;

	_sv current_cell = diagonal_cell + add_cbs(scores, query_bias);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);

	make_gap_mask(trace_mask, current_cell, vertical_gap, horizontal_gap);

	best = max(best, current_cell);

	row_counter.inc(best, current_cell);

	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);

	make_open_mask(trace_mask, open, vertical_gap, horizontal_gap);

	return current_cell;
}