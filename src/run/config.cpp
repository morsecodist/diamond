#include "workflow.h"
#include "../basic/config.h"
#include "../data/block.h"
#include "../data/taxonomy_nodes.h"
#include "../data/sequence_file.h"
#include "../search/hit.h"
#include "../util/async_buffer.h"

namespace Search {

Config::Config(bool dealloc) :
	dealloc(dealloc),
	self(config.self),
	db(nullptr),
	consumer(nullptr),
	query_file(nullptr),
	taxon_nodes(nullptr),
	taxonomy_scientific_names(nullptr)
{}


void Config::free()
{
	delete taxon_nodes;
	delete taxonomy_scientific_names;
	taxon_nodes = nullptr;
	taxonomy_scientific_names = nullptr;
	if (dealloc)
		delete db;
}

}