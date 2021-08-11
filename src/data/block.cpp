/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include "sequence_set.h"
#include "../basic/translate.h"
#include "../util/seq_file_format.h"
#include "block.h"
#include "../util/sequence/sequence.h"
#include "sequence_file.h"

using std::mutex;
using std::lock_guard;
using std::array;

static bool looks_like_dna(const Sequence& seq) {
	array<size_t, AMINO_ACID_COUNT> count;
	count.fill(0);
	for (size_t i = 0; i < seq.length(); ++i)
		++count[(int)seq[i]];
	return count[(int)value_traits.from_char('A')]
		+ count[(int)value_traits.from_char('C')]
		+ count[(int)value_traits.from_char('G')]
		+ count[(int)value_traits.from_char('T')]
		+ count[(int)value_traits.from_char('N')] == seq.length();
}

Block::Block(Alphabet alphabet):
	seqs_(alphabet),
	source_seqs_(Alphabet::STD),
	unmasked_seqs_(alphabet)
{
}

bool Block::empty() const {
	return seqs_.size() == 0;
}

void Block::convert_to_std_alph(size_t block_id)
{
	throw std::runtime_error("");
}

unsigned Block::source_len(unsigned block_id) const
{
	return align_mode.query_translated ? (unsigned)seqs_.reverse_translated_len(block_id * align_mode.query_contexts) : (unsigned)seqs_.length(block_id);
}

TranslatedSequence Block::translated(size_t block_id) const
{
	if (align_mode.query_translated)
		return seqs_.translated_seq(source_seqs_[block_id], block_id * align_mode.query_contexts);
	else
		return TranslatedSequence(seqs_[block_id]);
}

bool Block::long_offsets() const {
	return seqs_.raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
}

static size_t push_seq(SequenceSet& ss, SequenceSet& source_seqs, const vector<Letter>& seq, unsigned frame_mask, SequenceType seq_type)
{
	if (seq_type == SequenceType::amino_acid) {
		ss.push_back(seq.cbegin(), seq.cend());
		return seq.size();
	}
	else {
		source_seqs.push_back(seq.cbegin(), seq.cend());
		if (seq.size() < 2) {
			for (unsigned j = 0; j < 6; ++j)
				ss.fill(0, value_traits.mask_char);
			return 0;
		}
		vector<Letter> proteins[6];
		size_t n = Translator::translate(seq, proteins);

		unsigned bestFrames(Translator::computeGoodFrames(proteins, config.get_run_len((unsigned)seq.size() / 3)));
		for (unsigned j = 0; j < 6; ++j) {
			if ((bestFrames & (1 << j)) && (frame_mask & (1 << j)))
				ss.push_back(proteins[j].cbegin(), proteins[j].cend());
			else
				ss.fill(proteins[j].size(), value_traits.mask_char);
		}
		return n;
	}
}

Block::Block(std::list<TextInputFile>::iterator file_begin,
	std::list<TextInputFile>::iterator file_end,
	const Sequence_file_format& format,
	size_t max_letters,
	const Value_traits& value_traits,
	bool with_quals,
	bool lazy_masking,
	size_t modulo):
	seqs_(Alphabet::STD),
	source_seqs_(Alphabet::STD),
	unmasked_seqs_(Alphabet::STD)
{
	static constexpr size_t CHECK_FOR_DNA_COUNT = 10;
	size_t letters = 0, n = 0;
	vector<Letter> seq;
	string id;
	vector<char> qual;
	string id2;

	unsigned frame_mask = (1 << 6) - 1;
	if (config.query_strands == "plus")
		frame_mask = (1 << 3) - 1;
	else if (config.query_strands == "minus")
		frame_mask = ((1 << 3) - 1) << 3;

	std::list<TextInputFile>::iterator file_it = file_begin;
	bool read_success = true;

	while ((letters < max_letters || (n % modulo != 0)) && (read_success = format.get_seq(id, seq, *file_it, value_traits, with_quals ? &qual : nullptr))) {
		if (seq.size() > 0) {
			ids_.push_back(id.begin(), id.end());
			letters += push_seq(seqs_, source_seqs_, seq, frame_mask, value_traits.seq_type);
			if (with_quals)
				qual_.push_back(qual.begin(), qual.end());
			++n;
			if (seqs_.size() > (size_t)std::numeric_limits<int>::max())
				throw std::runtime_error("Number of sequences in file exceeds supported maximum.");
			if (n <= CHECK_FOR_DNA_COUNT && value_traits.seq_type == amino_acid && looks_like_dna(Sequence(seq)) && !config.ignore_warnings)
				throw std::runtime_error("The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.");
		}
		++file_it;
		if (file_it == file_end)
			file_it = file_begin;
	}
	ids_.finish_reserve();
	if (with_quals)
		qual_.finish_reserve();
	seqs_.finish_reserve();
	source_seqs_.finish_reserve();
	if (file_it != file_begin || (!read_success && ++file_it != file_end && format.get_seq(id, seq, *file_it, value_traits, nullptr)))
		throw std::runtime_error("Unequal number of sequences in paired read files.");
}

bool Block::fetch_seq_if_unmasked(size_t block_id, std::vector<Letter>& seq) {
	if (masked_[block_id])
		return false;
	{
		lock_guard<mutex> lck(mask_lock_);
		if (masked_[block_id])
			return false;
		seq.clear();
		Sequence s = seqs_[block_id];
		std::copy(s.data(), s.end(), std::back_inserter(seq));
		return true;
	}
}

void Block::write_masked_seq(size_t block_id, const std::vector<Letter>& seq) {
	lock_guard<mutex> lck(mask_lock_);
	if (masked_[block_id])
		return;
	std::copy(seq.begin(), seq.end(), seqs_.ptr(block_id));
	masked_[block_id] = true;
}

uint32_t Block::dict_id(size_t block, size_t block_id, SequenceFile& db) const
{
	string t;
	if (has_ids()) {
		const char* title = ids()[block_id];
		if (config.salltitles)
			t = title;
		else if (config.sallseqid)
			t = Util::Seq::all_seqids(title);
		else
			t = Util::Seq::seqid(title, config.short_seqids);
	}
	const Letter* seq = unmasked_seqs().empty() ? nullptr : unmasked_seqs()[block_id].data();
	return db.dict_id(block, block_id, block_id2oid(block_id), seqs().length(block_id), t.c_str(), seq);
}

void Block::soft_mask(const MaskingAlgo algo) {
	if (soft_masked_)
		return;
	if (soft_masking_table_.blank())
		mask_seqs(seqs_, Masking::get(), true, algo, &soft_masking_table_);
	else
		soft_masking_table_.apply(seqs_);
	soft_masked_ = true;
}

void Block::remove_soft_masking() {
	if (!soft_masked_)
		return;
	soft_masking_table_.remove(seqs_);
	soft_masked_ = false;
}

bool Block::soft_masked() const {
	return soft_masked_;
}
