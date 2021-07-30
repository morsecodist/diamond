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

#pragma once
#include <list>
#include <vector>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../stats/hauser_correction.h"
#include "../basic/statistics.h"
#include "../basic/config.h"
#include "../data/sequence_set.h"
#include "../stats/cbs.h"
#include "flags.h"

struct DpTarget
{
	struct CarryOver {
		CarryOver() :
			i1(0), j1(0), ident(0), len(0)
		{}
		CarryOver(int i1, int j1, int ident, int len) :
			i1(i1), j1(j1), ident(ident), len(len)
		{}
		int i1, j1, ident, len;
	};
	enum { BLANK = -1 };
	int get_cols(int qlen) const {
		int pos = std::max(d_end - 1, 0) - (d_end - 1);
		const int d0 = d_begin;
		const int j1 = std::min(qlen - 1 - d0, (int)(seq.length() - 1)) + 1;
		return j1 - pos;
	}
	DpTarget():
		d_begin(),
		d_end(),
		target_idx(BLANK),
		cols(),
		matrix(nullptr)
	{}
	DpTarget(const Sequence &seq, int d_begin, int d_end, int target_idx, int qlen, const Stats::TargetMatrix* matrix = nullptr, const CarryOver& carry_over = CarryOver()) :
		seq(seq),
		d_begin(d_begin),
		d_end(d_end),
		target_idx(target_idx),
		cols(get_cols(qlen)),
		carry_over(carry_over),
		matrix(matrix)
	{
	}
	DpTarget(const Sequence& seq, int target_idx, const Stats::TargetMatrix* matrix = nullptr, const CarryOver& carry_over = CarryOver()):
		seq(seq),
		d_begin(),
		d_end(),
		target_idx(target_idx),
		cols(),
		carry_over(carry_over),
		matrix(matrix)
	{}
	DpTarget(const std::pair<const Letter*, size_t> seq) :
		seq(seq.first, seq.second),
		d_begin(),
		d_end(),
		target_idx(BLANK),
		cols(),
		matrix(nullptr)
	{
	}
	int left_i1() const
	{
		return std::max(d_end - 1, 0);
	}
	int band() const {
		return d_end - d_begin;
	}
	bool operator<(const DpTarget &x) const
	{
		const int i = left_i1(), j = x.left_i1(), b1 = band(), b2 = x.band(), bin_b1 = b1 / config.band_bin, bin_b2 = b2 / config.band_bin,
			t1 = cols, t2 = x.cols, bin_t1 = t1 / config.col_bin, bin_t2 = t2 / config.col_bin;
		return bin_b1 < bin_b2 || (bin_b1 == bin_b2 && (bin_t1 < bin_t2 || (bin_t1 == bin_t2 && i < j)));
		//return i < j || (i == j && (target_idx < x.target_idx || (target_idx == x.target_idx && d_begin < x.d_begin)));
	}
	bool blank() const {
		return target_idx == BLANK;
	}
	bool adjusted_matrix() const {
		return matrix != nullptr;
	}
	int matrix_scale() const {
		return adjusted_matrix() ? config.cbs_matrix_scale : 1;
	}
	Sequence seq;
	int d_begin, d_end, target_idx, cols;
	CarryOver carry_over;
	const Stats::TargetMatrix* matrix;
};

struct DpStat
{
	DpStat():
		gross_cells(0),
		net_cells(0)
	{}
	DpStat& operator+=(DpStat &x)
	{
		mtx_.lock();
		gross_cells += x.gross_cells;
		net_cells += x.net_cells;
		mtx_.unlock();
		return *this;
	}
	size_t gross_cells, net_cells;
private:
	std::mutex mtx_;
};

extern DpStat dp_stat;

namespace DP {

enum { BINS = 3};

struct Traceback {};
struct ScoreOnly {};

using Targets = std::array<std::vector<DpTarget>, BINS>;

struct NoCBS {
	constexpr void* operator[](int i) const { return nullptr; }
};
	
namespace Swipe {

//DECL_DISPATCH(std::list<Hsp>, swipe, (const sequence &query, const sequence *subject_begin, const sequence *subject_end, int score_cutoff))

}

namespace BandedSwipe {

DECL_DISPATCH(std::list<Hsp>, swipe, (const Sequence& query, const Targets& targets, const Frame frame, const Bias_correction* composition_bias, const DP::Flags flags, const HspValues v, Statistics& stat))
DECL_DISPATCH(std::list<Hsp>, swipe_set, (const Sequence& query, const SequenceSet::ConstIterator begin, const SequenceSet::ConstIterator end, const Frame frame, const Bias_correction* composition_bias, const DP::Flags flags, const HspValues v, Statistics& stat))
DECL_DISPATCH(unsigned, bin, (HspValues v, int query_len, int score, int ungapped_score, size_t dp_size, unsigned score_width))

}

}

DECL_DISPATCH(std::list<Hsp>, banded_3frame_swipe, (const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, DpStat &stat, bool score_only, bool parallel))
