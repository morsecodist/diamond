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

#include "global_ranking.h"
#define _REENTRANT
#include "../../lib/ips4o/ips4o.hpp"
#include "../../search/hit.h"
#include "../../util/data_structures/deque.h"
#include "../../util/util.h"
#include "../../util/algo/algo.h"

using std::endl;
using std::thread;
using SeedHits = Search::Config::RankingBuffer;

namespace Extension { namespace GlobalRanking {

static void update_query(SeedHits::Iterator begin, SeedHits::Iterator end, vector<Hit>& hits, vector<Hit>& merged, std::atomic_size_t& merged_count, Search::Config& cfg) {
	const size_t N = config.global_ranking_targets;
	hits.clear();
	merged.clear();
	const SequenceSet& target_seqs = cfg.target->seqs();
	auto get_target = [&target_seqs](const Search::Hit& hit) { return target_seqs.local_position((uint64_t)hit.subject_).first; };
	auto it = merge_keys(begin, end, get_target);
	while (it.good()) {
		uint16_t score = 0;
		for (SeedHits::Iterator i = it.begin(); i != it.end(); ++i)
			score = std::max(score, i->score_);
		hits.emplace_back((uint32_t)cfg.target->block_id2oid(it.key()), score);
		++it;
	}
	std::sort(hits.begin(), hits.end());
	const size_t q = begin->query_;
	vector<Hit>::iterator table_begin = cfg.ranking_table->begin() + q * N, table_end = table_begin + N;
	while (table_end > table_begin && (table_end - 1)->score == 0) --table_end;
	merged_count += Util::Algo::merge_capped(table_begin, table_end, hits.begin(), hits.end(), N, std::back_inserter(merged));
	std::copy(merged.begin(), merged.end(), table_begin);
}

void update_table(Search::Config& cfg) {
	SeedHits& hits = *cfg.global_ranking_buffer;
	log_stream << "Seed hits = " << hits.size() << endl;
	task_timer timer("Sorting seed hits");
	ips4o::parallel::sort(hits.begin(), hits.end(), Search::Hit::CmpQueryTarget(), config.threads_);
	timer.go("Processing seed hits");
	/*AsyncKeyMerger<SeedHits::Iterator, Search::Hit::SourceQuery> it(hits.begin(), hits.end(), { align_mode.query_contexts });
	std::atomic_size_t merged_count(0);
	auto worker = [&it, &cfg, &merged_count] {
		vector<Hit> hits, merged;
		for (;;) {
			auto r = ++it;
			if (r.first == r.second)
				return;
			update_query(r.first, r.second, hits, merged, merged_count, cfg);
		}
	};*/
	std::atomic_size_t merged_count(0);
	auto worker = [&cfg, &merged_count](SeedHits::Iterator begin, SeedHits::Iterator end) {
		auto it = merge_keys(begin, end, ::Search::Hit::Query());
		while (it.good()) {
			vector<Hit> hits, merged;
			update_query(it.begin(), it.end(), hits, merged, merged_count, cfg);
			++it;
		}
	};
	vector<thread> threads;
	auto p = Util::Algo::partition_table(hits.begin(), hits.end(), config.threads_, ::Search::Hit::Query());
	for (size_t i = 0; i < p.size() - 1; ++i)
		threads.emplace_back(worker, p[i], p[i + 1]);
	for (thread& t : threads)
		t.join();
	timer.finish();
	log_stream << "Merged targets = " << merged_count << endl;
}

}}