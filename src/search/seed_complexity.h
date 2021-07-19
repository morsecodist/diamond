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
#include <math.h>
#include <stddef.h>
#include "../basic/value.h"
#include "../basic/shape.h"
#include "../basic/reduction.h"
#include "../basic/config.h"
#include "../data/seed_array.h"
#include "../util/data_structures/double_array.h"
#include "../run/config.h"

struct SeedComplexity
{

	static void init(const Reduction &r)
	{
		double p[20];
		for (size_t i = 0; i < 20; ++i)
			p[i] = 0;
		for (size_t i = 0; i < 20; ++i)
			p[r(i)] += score_matrix.background_freqs()[i];
		for (size_t i = 0; i < 20; ++i)
			prob_[i] = log(p[r(i)]);
		for (size_t i = 20; i < AMINO_ACID_COUNT; ++i)
			prob_[i] = 1000;
	}

	static bool complex(const Letter *seq, Shape shape)
	{
		double p = 0;
		for (unsigned i = 0; i < shape.weight_; ++i)
			p += prob_[(size_t)seq[shape.positions_[i]]];
		return p <= -config.freq_treshold;
	}

private:

	static double prob_[AMINO_ACID_COUNT];

};

namespace Search {

bool seed_is_complex(const Letter* seq, const Shape& shape, const double cut);
void mask_seeds(const Shape& shape, const SeedPartitionRange& range, DoubleArray<SeedArray::Entry::Value>* query_seed_hits, DoubleArray<SeedArray::Entry::Value>* ref_seed_hits, Search::Config& cfg);

}