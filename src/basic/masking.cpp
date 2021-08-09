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


#include <math.h>
#include <algorithm>
#include <atomic>
#include <numeric>
#include "masking.h"
#include "../lib/tantan/LambdaCalculator.hh"
#include "../util/tantan.h"
#include "../lib/blast/blast_filter.h"

using std::unique_ptr;
using std::atomic;
using std::thread;
using std::vector;
using std::pair;

const EMap<MaskingAlgo> EnumTraits<MaskingAlgo>::to_string{ {MaskingAlgo::NONE, "None"}, {MaskingAlgo::SEG, "SEG"}, {MaskingAlgo::TANTAN, "tantan"} };
const SEMap<MaskingAlgo> EnumTraits<MaskingAlgo>::from_string{
	{"0", MaskingAlgo::NONE}
};
const SEMap<MaskingMode> EnumTraits<MaskingMode>::from_string{
	{"0", {MaskingMode::NONE, false}},
	{"none", MaskingMode::NONE},
	{"1", {MaskingMode::TANTAN, false}},
	{"tantan", MaskingMode::TANTAN},
	{"seg", MaskingMode::BLAST_SEG}
};

unique_ptr<Masking> Masking::instance;
const int8_t Masking::bit_mask = (int8_t)128;

void MaskingTable::add(const size_t block_id, const int begin, const int end, Letter* seq) {
	{
		std::lock_guard<std::mutex> lock(mtx_);
		entry_.emplace_back(block_id, begin);
		seqs_.push_back(seq + begin, seq + end);
	}
	std::fill(seq + begin, seq + end, MASK_LETTER);
}

void MaskingTable::remove(SequenceSet& seqs) const {
	for (size_t i = 0; i < entry_.size(); ++i) {
		Letter* ptr = seqs.ptr(entry_[i].block_id) + entry_[i].begin;
		std::copy(seqs_[i].data(), seqs_[i].end(), ptr);
	}
}

void MaskingTable::apply(SequenceSet& seqs) const {
	for (size_t i = 0; i < entry_.size(); ++i) {
		Letter* ptr = seqs.ptr(entry_[i].block_id) + entry_[i].begin;
		std::fill(ptr, ptr + seqs_[i].length(), MASK_LETTER);
	}
}

Masking::Masking(const Score_matrix &score_matrix)
{
	const unsigned n = value_traits.alphabet_size;
	int int_matrix[20][20], *int_matrix_ptr[20];
	std::copy(int_matrix, int_matrix + 20, int_matrix_ptr);
	for (size_t i = 0; i < 20; ++i)
		for (size_t j = 0; j < 20; ++j)
			int_matrix[i][j] = score_matrix(i, j);
	cbrc::LambdaCalculator lc;
	lc.calculate(int_matrix_ptr, 20);
	
	const double lambda = lc.lambda(); // 0.324032
	for (size_t i = 0; i < size; ++i) {
		mask_table_x_[i] = value_traits.mask_char;
		mask_table_bit_[i] = (int8_t)i | bit_mask;
		for (size_t j = 0; j < size; ++j)
			if (i < n && j < n) {
				likelihoodRatioMatrixf_[i][j] = (float)exp(lambda * score_matrix(i, j));
			}
	}
	std::copy(likelihoodRatioMatrixf_, likelihoodRatioMatrixf_ + size, probMatrixPointersf_);

	blast_seg_ = SegParametersNewAa();
}

Masking::~Masking() {
	SegParametersFree(blast_seg_);
}

void Masking::operator()(Letter *seq, size_t len, MaskingAlgo algo) const
{
	if(algo == MaskingAlgo::TANTAN)
		Util::tantan::mask(seq, (int)len, (const float**)probMatrixPointersf_, 0.005f, 0.05f, 1.0f / 0.9f, (float)config.tantan_minMaskProb, mask_table_x_);
	else {
		BlastSeqLoc* seg_locs;
		SeqBufferSeg((uint8_t*)seq, len, 0u, blast_seg_, &seg_locs);
		unsigned nMasked = 0;

		if (seg_locs) {
			BlastSeqLoc* l = seg_locs;
			do {
				for (signed i = l->ssr->left; i <= l->ssr->right; i++) {
					nMasked++;
					seq[i] = value_traits.mask_char;
				}
			} while ((l = l->next) != 0);
			BlastSeqLocFree(seg_locs);
		}
	}
}

void Masking::mask_bit(Letter *seq, size_t len) const
{
	Util::tantan::mask(seq, (int)len, (const float**)probMatrixPointersf_, 0.005f, 0.05f, 1.0f / 0.9f, (float)config.tantan_minMaskProb, mask_table_bit_);
}

void Masking::bit_to_hard_mask(Letter *seq, size_t len, size_t &n) const
{
	for (size_t i = 0; i < len; ++i)
		if (seq[i] & bit_mask) {
			seq[i] = value_traits.mask_char;
			++n;
		}
}

void Masking::remove_bit_mask(Letter *seq, size_t len) const
{
	for (size_t i = 0; i < len; ++i)
		if (seq[i] & bit_mask)
			seq[i] &= ~bit_mask;
}

void mask_worker(atomic<size_t> *next, SequenceSet *seqs, const Masking *masking, bool hard_mask, const MaskingAlgo algo)
{
	size_t i;
	while ((i = (*next)++) < seqs->size()) {
		seqs->convert_to_std_alph(i);
		if (hard_mask)
			masking->operator()(seqs->ptr(i), seqs->length(i), algo);
		else
			masking->mask_bit(seqs->ptr(i), seqs->length(i));
	}
}

size_t mask_seqs(SequenceSet &seqs, const Masking &masking, bool hard_mask, const MaskingAlgo algo)
{
	if (algo == MaskingAlgo::NONE)
		return 0;
	vector<thread> threads;
	atomic<size_t> next(0);
	for (size_t i = 0; i < config.threads_; ++i)
		threads.emplace_back(mask_worker, &next, &seqs, &masking, hard_mask, algo);
	for (auto &t : threads)
		t.join();
	size_t n = 0;
	for (size_t i = 0; i < seqs.size(); ++i)
		n += std::count(seqs[i].data(), seqs[i].end(), value_traits.mask_char);
	seqs.alphabet() = Alphabet::STD;
	return n;
}

void mask_motifs(Letter* seq, const size_t len, const size_t block_id, MaskingTable& table) {
	if (len < MOTIF_LEN)
		return;
	vector<pair<ptrdiff_t, ptrdiff_t>> pos;
	KmerIterator<MOTIF_LEN> it(Sequence(seq, len));
	while (it.good()) {
		if (motif_table.find(*it) != motif_table.end()) {
			const ptrdiff_t p = it - seq;
			if (!pos.empty() && p <= pos.back().second)
				pos.back().second = p + MOTIF_LEN;
			else
				pos.emplace_back(p, p + MOTIF_LEN);
		}
	}
	const ptrdiff_t n = std::accumulate(pos.cbegin(), pos.cend(), (ptrdiff_t)0, [](const ptrdiff_t s, const pair<ptrdiff_t, ptrdiff_t>& r) { return s + r.second - r.first; });
	
	if ((double)n / len >= 0.5)
		return;

	for (auto i = pos.cbegin(); i != pos.cend(); ++i)
		table.add(block_id, (int)i->first, (int)i->second, seq);
}
