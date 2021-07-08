/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

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

#pragma once

static inline uint8_t cmp_mask(int x, int y) {
	return x == y;
}

static inline int blend(int v, int w, int mask) {
	return mask ? w : v;
}

template<typename Sv>
struct DummyIdMask {
	DummyIdMask(const Letter q, const Sv& t)
	{}
};

template<typename Sv>
struct VectorIdMask {
	VectorIdMask(const Letter q, const Sv& t) :
		mask(blend(Sv(0), Sv(1), Sv(typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score(q)) == t))
	{}
	const Sv mask;
};

template<typename Sv>
struct ForwardCell : public Sv {
	ForwardCell() :
		Sv(),
		ident(),
		len()
	{}
	ForwardCell(const Sv& v) :
		Sv(v),
		ident(),
		len()
	{}
	Sv ident, len;
};

template<>
struct ForwardCell<int32_t> {
	int32_t v, ident, len;
	operator int32_t() const {
		return v;
	}
	ForwardCell(const int32_t v) :
		v(v),
		ident(0),
		len(0)
	{}
	ForwardCell() :
		v(0),
		ident(0),
		len(0)
	{}
	ForwardCell& operator-=(int32_t x) {
		v -= x;
		return *this;
	}
	ForwardCell& operator+=(int32_t x) {
		v += x;
		return *this;
	}
	void max(const ForwardCell& x) {
		v = std::max(v, x.v);
	}
};

template<typename Sv>
struct BackwardCell : public Sv {
	BackwardCell() :
		Sv(),
		mismatch(),
		gapopen()
	{}
	BackwardCell(const Sv& v) :
		Sv(v),
		mismatch(),
		gapopen()
	{}
	Sv mismatch, gapopen;
};

template<>
struct BackwardCell<int32_t> {
	int32_t v, mismatch, gapopen;
	operator int32_t() const {
		return v;
	}
	BackwardCell(const int32_t v) :
		v(v),
		mismatch(0),
		gapopen(0)
	{}
	BackwardCell() :
		v(0),
		mismatch(0),
		gapopen(0)
	{}
	BackwardCell& operator-=(int32_t x) {
		v -= x;
		return *this;
	}
	BackwardCell& operator+=(int32_t x) {
		v += x;
		return *this;
	}
	void max(const BackwardCell& x) {
		v = std::max(v, x.v);
	}
};

template<typename Sv>
FORCE_INLINE static int extract_stats(const Sv&, int) {
	return 0;
}

template<typename Sv>
FORCE_INLINE static int extract_stats(const BackwardCell<Sv>& v, int channel) {
	return ::DISPATCH_ARCH::ScoreTraits<Sv>::int_score(extract_channel(v.gapopen, channel));
}

template<typename Sv>
FORCE_INLINE static void update_stats(const Sv&, const Sv&, const Sv&, const DummyIdMask<Sv>&) {
}

template<typename Sv>
FORCE_INLINE static void update_stats(ForwardCell<Sv>& current_cell, ForwardCell<Sv>& horizontal_gap, ForwardCell<Sv>& vertical_gap, const VectorIdMask<Sv>& id_mask) {
	const Sv zero = Sv(), one = Sv(1), zero_mask = current_cell == zero;
	current_cell.ident += id_mask.mask;
	current_cell.ident = blend(current_cell.ident, zero, zero_mask);
	current_cell.len += one;
	current_cell.len = blend(current_cell.len, zero, zero_mask);
	horizontal_gap.len += one;
	vertical_gap.len += one;
}

template<typename Sv>
FORCE_INLINE static void update_stats(BackwardCell<Sv>& current_cell, BackwardCell<Sv>& horizontal_gap, BackwardCell<Sv>& vertical_gap, const VectorIdMask<Sv>& id_mask) {
	const Sv zero = Sv(), one = Sv(1), zero_mask = current_cell == zero;
	current_cell.mismatch += one - id_mask.mask;
	current_cell.mismatch = blend(current_cell.mismatch, zero, zero_mask);
	current_cell.gapopen = blend(current_cell.gapopen, zero, zero_mask);
}

template<typename Sv>
FORCE_INLINE static void update_open(const Sv&) {
}

template<typename Sv>
FORCE_INLINE static void update_open(BackwardCell<Sv>& v) {
	v.gapopen += Sv(1);
}

template<typename Sv>
FORCE_INLINE static void set_max(ForwardCell<Sv>& v, const ForwardCell<Sv>& x) {
	v.max(x);
	const Sv mask = v == x;
	v.ident = blend(v.ident, x.ident, mask);
	v.len = blend(v.len, x.len, mask);
}

template<typename Sv>
FORCE_INLINE static void set_max(BackwardCell<Sv>& v, const BackwardCell<Sv>& x) {
	v.max(x);
	const Sv mask = v == x;
	v.mismatch = blend(v.mismatch, x.mismatch, mask);
	v.gapopen = blend(v.gapopen, x.gapopen, mask);
}
