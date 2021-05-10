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

using std::endl;
using std::thread;
using SeedHits = Deque<Search::Hit, Async>;

namespace Extension { namespace GlobalRanking {

static void update_query(SeedHits::Iterator begin, SeedHits::Iterator end, vector<Hit>& hits, Search::Config& cfg) {
	hits.clear();
	const SequenceSet& target_seqs = cfg.target->seqs();
	auto get_target = [&target_seqs](const Search::Hit& hit) { return target_seqs.local_position((uint64_t)hit.subject_).first; };
	auto it = merge_keys(begin, end, get_target);
	while (it.good()) {
		uint16_t score = 0;
		for (SeedHits::Iterator i = it.begin(); i != it.end(); ++i)
			score = std::max(score, i->score_);
		hits.emplace_back((uint32_t)it.key(), score);
		++it;
	}
	std::sort(hits.begin(), hits.end());
}

void update_table(Search::Config& cfg) {
	SeedHits& hits = *cfg.global_ranking_buffer;
	log_stream << "Seed hits = " << hits.size() << endl;
	task_timer timer("Sorting seed hits");
	ips4o::parallel::sort(hits.begin(), hits.end(), Search::Hit::CmpQueryTarget(), config.threads_);
	timer.go("Processing seed hits");
	AsyncKeyMerger<SeedHits::Iterator, Search::Hit::SourceQuery> it(hits.begin(), hits.end(), { align_mode.query_contexts });
	auto worker = [&it, &cfg] {
		vector<Hit> hits;
		for (;;) {
			auto r = ++it;
			if (r.first == r.second)
				return;
			update_query(r.first, r.second, hits, cfg);
		}
	};
	vector<thread> threads;
	for (size_t i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (thread& t : threads)
		t.join();
}

}}