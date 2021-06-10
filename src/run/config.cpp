#include "config.h"
#include "../basic/config.h"
#include "../data/block.h"
#include "../data/taxonomy_nodes.h"
#include "../data/sequence_file.h"
#include "../search/hit.h"
#include "../util/async_buffer.h"
#include "../search/hit.h"
#include "../util/data_structures/deque.h"
#include "../align/global_ranking/global_ranking.h"

using std::endl;

namespace Search {

Config::Config() :
	self(config.self),
	track_aligned_queries(false),
	db(nullptr),
	out(nullptr),
	query_file(nullptr),
	taxon_nodes(nullptr),
	taxonomy_scientific_names(nullptr)
{
	if (!config.iterate.empty()) {
		if (config.multiprocessing)
			throw std::runtime_error("Iterated search is not compatible with --multiprocessing.");
		for (const string& s : config.iterate)
			sensitivity.push_back(from_string<Sensitivity>(s));
		message_stream << "Running iterated search mode with sensitivity steps:";
		for (Sensitivity s : sensitivity)
			message_stream << ' ' << to_string(s);
		message_stream << endl;
		track_aligned_queries = true;
	}
	else
		sensitivity.push_back(config.sensitivity);
	if (!config.unaligned.empty() || !config.aligned_file.empty())
		track_aligned_queries = true;
}

Config::~Config() {

}

void Config::free()
{
	delete taxon_nodes;
	delete taxonomy_scientific_names;
	taxon_nodes = nullptr;
	taxonomy_scientific_names = nullptr;
}

}