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

#include <array>
#include <atomic>
#include "seed_complexity.h"
#include "../data/block.h"
#include "../util/algo/join_result.h"

using std::array;
using std::endl;

extern double lnfact[];

bool Search::seed_is_complex(const Letter* seq, const Shape& shape, const double cut)
{
	array<unsigned, TRUE_AA> count;
	count.fill(0);
	for (unsigned i = 0; i < shape.weight_; ++i) {
		const Letter l = letter_mask(seq[shape.positions_[i]]);
		++count[Reduction::reduction(l)];
	}
	double entropy = lnfact[shape.weight_];
	for (unsigned i = 0; i < Reduction::reduction.size(); ++i)
		entropy -= lnfact[count[i]];
	return entropy / shape.weight_ >= cut;
}

bool Search::seed_is_complex_unreduced(Letter* seq, const Shape& shape, const double cut, const bool mask_seeds)
{
	array<unsigned, AMINO_ACID_COUNT> count;
	count.fill(0);
	double f = 0.0;
	for (unsigned i = 0; i < shape.weight_; ++i) {
		const Letter l = letter_mask(seq[shape.positions_[i]]);
		const unsigned r = Reduction::reduction(l);
		if(r == MASK_LETTER) {
			if(mask_seeds) *seq |= SEED_MASK;
			return false;
		}
		++count[(int)l];
		f += Reduction::reduction.freq(r);
	}
	double entropy = lnfact[shape.weight_];
	for (unsigned i = 0; i < AMINO_ACID_COUNT; ++i)
		entropy -= lnfact[count[i]];
	if (entropy / shape.weight_ < cut || f / shape.weight_ > config.max_seed_freq) {
		if (mask_seeds) *seq |= SEED_MASK;
		return false;
	}
	return true;	
}

/*void Search::mask_seeds(const Shape& shape, const SeedPartitionRange& range, DoubleArray<SeedArray::Entry::Value>* query_seed_hits, DoubleArray<SeedArray::Entry::Value>* ref_seed_hits, Search::Config& cfg)
{
	task_timer timer("Masking low complexity seeds");
	SequenceSet& query_seqs = cfg.query->seqs();
	std::atomic_uint32_t seedp(range.begin());
	std::atomic_size_t seed_count(0), masked_seed_count(0), query_count(0), target_count(0);
	const double cut = cfg.seed_complexity_cut;

	auto worker = [&] {
		unsigned p;
		size_t sc(0), msc(0), qc(0), tc(0);
		while ((p = seedp++) < range.end()) {
			for (auto it = JoinIterator<SeedArray::Entry::Value>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it;) {
				++sc;
				const Range<SeedArray::Entry::Value*> query_hits = *it.r;
				const Letter* seq = query_seqs.data(*query_hits.begin());
				if (!seed_is_complex(seq, shape, cut)) {
					++msc;
					qc += it.r->size();
					tc += it.s->size();
					for (SeedArray::Entry::Value* i = query_hits.begin(); i < query_hits.end(); ++i) {
						Letter* p = query_seqs.data(*i);
						*p |= SEED_MASK;
					}
					it.erase();
				}
				else
					++it;
			}
		}
		seed_count += sc;
		masked_seed_count += msc;
		query_count += qc;
		target_count += tc;
	};

	vector<std::thread> threads;
	for (size_t i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (auto& i : threads)
		i.join();
	timer.finish();
	verbose_stream << "Masked seeds: " << masked_seed_count << '/' << seed_count << endl;
	verbose_stream << "Masked positions (query): " << query_count << endl;
	verbose_stream << "Masked positions (target): " << target_count << endl;
}
*/