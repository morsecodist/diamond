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

namespace Search {

Config::Config() :
	self(config.self),
	db(nullptr),
	out(nullptr),
	query_file(nullptr),
	taxon_nodes(nullptr),
	taxonomy_scientific_names(nullptr)
{}

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