#include <mutex>
#include "../basic/config.h"
#include "../data/sequence_file.h"
#include "../util/string/tokenizer.h"
#include "../util/string/tsv.h"
#include "../dp/dp.h"
#include "../output/output_format.h"
#include "../output/output.h"
#include "../basic/masking.h"

using std::unique_ptr;
using std::endl;
using std::mutex;
using std::lock_guard;
using std::thread;
using std::list;

static TextBuffer* view_query(const string& query_acc, const string& buf, SequenceFile& query_file, SequenceFile& target_file, Search::Config& cfg, Statistics& stats) {
	const vector<string> target_acc = Util::Tsv::extract_column(buf, 1);
	SequenceSet targets = target_file.seqs_by_accession(target_acc.begin(), target_acc.end());
	vector<Letter> query = query_file.seq_by_accession(query_acc);
	if (cfg.query_masking != MaskingAlgo::NONE)
		Masking::get()(query.data(), query.size(), cfg.query_masking);
	if (cfg.target_masking != MaskingAlgo::NONE)
		for (size_t i = 0; i < targets.size(); ++i)
			Masking::get()(targets.ptr(i), targets.length(i), cfg.target_masking);

	const HspValues v = HspValues::COORDS | HspValues::IDENT | HspValues::LENGTH;
	DP::Targets dp_targets;
	for (size_t i = 0; i < target_acc.size(); ++i)
		dp_targets[DP::BandedSwipe::bin(v, query.size(), 0, 0, 0, 0)].emplace_back(targets[i], i);

	list<Hsp> hsp = DP::BandedSwipe::swipe(Sequence(query), dp_targets, Frame(0), nullptr, DP::Flags::FULL_MATRIX, v, stats);
	hsp.sort();

	Blast_tab_format fmt;
	TranslatedSequence query_seq;
	TextBuffer* out = new TextBuffer();
	for (Hsp& h : hsp) {
		h.query_source_range = h.query_range;
		fmt.print_match(HspContext(h,
			0,
			query_seq,
			query_acc.c_str(),
			0,
			0,
			target_acc[h.swipe_target].c_str(),
			0,
			0,
			Sequence()), cfg, *out);
	}
	return out;
}

void view_tsv() {
	if (config.input_ref_file.size() > 1)
		throw std::runtime_error("Too many arguments for --in.");
	if (config.database.empty())
		throw std::runtime_error("Missing argument: database file (-d)");
	if (config.query_file.empty())
		throw std::runtime_error("Missing argument: query file (-q)");
	if (config.query_file.size() > 1)
		throw std::runtime_error("Too many arguments for query file (--query/-q)");

	task_timer timer("Opening the database file");
	unique_ptr<SequenceFile> db(SequenceFile::auto_create(config.database, SequenceFile::Flags::NO_FASTA));
	score_matrix = Score_matrix("blosum62", -1, -1, 1, 0);
	score_matrix.set_db_letters(config.db_size ? config.db_size : db->letters());
	Masking::instance = unique_ptr<Masking>(new Masking(score_matrix));

	timer.go("Opening the query file");
	unique_ptr<SequenceFile> query_file(SequenceFile::auto_create(config.query_file.front(), SequenceFile::Flags::NO_FASTA));

	if (db->type() != SequenceFile::Type::BLAST || query_file->type() != SequenceFile::Type::BLAST)
		throw std::runtime_error("BLAST database required.");

	timer.go("Opening the input file");
	TextInputFile in(config.input_ref_file.front());

	timer.go("Opening the output file");
	OutputFile output_file(config.output_file);
	OutputSink::instance = unique_ptr<OutputSink>(new OutputSink(0, &output_file));

	timer.go("Computing alignments");
	size_t query_idx = 0;
	mutex mtx;
	Search::Config cfg;

	auto worker = [&] {
		string query, buf;
		size_t q;
		Statistics stats;
		for (;;) {
			{
				lock_guard<mutex> lock(mtx);
				query = Util::Tsv::fetch_block(in, buf);
				q = query_idx++;
			}
			if (query.empty())
				return;
			TextBuffer* out = view_query(query, buf, *query_file, *db, cfg, stats);
			OutputSink::instance->push(q, out);
		}
	};

	vector<thread> threads;
	for (size_t i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (auto& t : threads)
		t.join();

	timer.go("Closing the output file");
	output_file.close();
}