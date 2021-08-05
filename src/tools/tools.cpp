#include <array>
#include <iostream>
#include <vector>
#include "tsv_record.h"
#include "../basic/config.h"
#include "../basic/value.h"
#include "../util/util.h"
#include "../basic/sequence.h"
#include "../util/seq_file_format.h"
#include "../util/sequence/sequence.h"
#include "../stats/cbs.h"
#include "../util/algo/MurmurHash3.h"
#include "../util/sequence/sequence.h"
#include "../data/sequence_file.h"
#include "../search/search.h"
#define _REENTRANT
#include "../lib/ips4o/ips4o.hpp"
#include "../basic/masking.h"

using std::array;
using std::cout;
using std::endl;
using std::vector;
using std::unique_ptr;
using std::pair;

void filter_blasttab() {
	TextInputFile in("");
	TSVRecord r;
	string query;
	size_t query_hit;
	while (in >> r) {
		if (r.qseqid != query) {
			query = r.qseqid;
			query_hit = 0;
		}
		else
			++query_hit;
		if(query_hit < config.max_alignments && r.evalue <= config.max_evalue)
			cout << r << endl;
	}
}

void split() {
	TextInputFile in(config.single_query_file());
	string id;
	vector<Letter> seq;
	size_t n = 0, f = 0, b = (size_t)(config.chunk_size * 1e9);
	OutputFile *out = new OutputFile(std::to_string(f) + ".faa.gz", Compressor::ZLIB);
	while (FASTA_format().get_seq(id, seq, in, value_traits)) {
		if (n >= b) {
			out->close();
			delete out;
			out = new OutputFile(std::to_string(++f) + ".faa.gz", Compressor::ZLIB);
			n = 0;
		}
		string blast_id = Util::Seq::seqid(id.c_str(), false);
		Util::Seq::format(Sequence(seq), blast_id.c_str(), nullptr, *out, "fasta", amino_acid_traits);
		n += seq.size();
	}
	out->close();
	delete out;
}

void composition() {
	TextInputFile in(config.single_query_file());
	string id;
	vector<Letter> seq;
	while (FASTA_format().get_seq(id, seq, in, value_traits)) {
		auto c = Stats::composition(seq);
		for (double x : c)
			std::cout << x << '\t';
		std::cout << endl;
	}
}

void hash_seqs() {
	TextInputFile f(config.query_file.front());
	FASTA_format fmt;
	string id;
	vector<Letter> seq;
	while (fmt.get_seq(id, seq, f, amino_acid_traits)) {
		array<char, 16> hash;
		hash.fill('\0');
		MurmurHash3_x64_128(seq.data(), seq.size(), hash.data(), hash.data());
		cout << Util::Seq::seqid(id.c_str(), false) << '\t' << hex_print(hash.data(), 16) << endl;
	}
	f.close();
}

static double freq(const string& s, const Reduction& r) {
	double f = 0.0;
	for (char c : s) {
		f += r.freq(r(amino_acid_traits.from_char(c)));
	}
	return f / s.length();
}

void list_seeds() {
	struct Callback {
		bool operator()(uint64_t seed, size_t, unsigned, size_t) {
			seeds.push_back(seed);
			return true;
		};
		void finish() {}
		vector<uint64_t>& seeds;
	};
	unique_ptr<SequenceFile> db(SequenceFile::auto_create(config.database));
	unique_ptr<Block> block(db->load_seqs(SIZE_MAX));
	mask_seqs(block->seqs(), Masking::get(), true, MaskingAlgo::TANTAN);
	vector<uint64_t> seeds;
	seeds.reserve(block->seqs().letters());
	PtrVector<Callback> cb;
	cb.push_back(new Callback{ seeds });
	auto parts = block->seqs().partition(1);
	::shapes = ShapeConfig(config.shape_mask.empty() ? shape_codes.at(Sensitivity::DEFAULT) : config.shape_mask, config.shapes);
	Reduction::reduction = Reduction("A R N D C Q E G H I L K M F P S T W Y V");
	enum_seeds(&block->seqs(), cb, parts, 0, 1, &no_filter, SeedEncoding::SPACED_FACTOR, nullptr, false, false);
	ips4o::parallel::sort(seeds.begin(), seeds.end());

	auto it = merge_keys(seeds.begin(), seeds.end(), [](uint64_t seed) {return seed; });
	vector<pair<uint64_t, uint64_t>> counts;
	while (it.good()) {
		counts.push_back({ it.count(), it.key() });
		++it;
	}
	ips4o::parallel::sort(counts.begin(), counts.end());

	auto end = std::min(counts.rbegin() + config.query_count, counts.rend());
	string s;
	for (auto i = counts.rbegin(); i != end; ++i) {
		s = Reduction::reduction.decode_seed(i->second, shapes[0].weight_);
		cout << i->first << '\t' << s << endl;
	}
}