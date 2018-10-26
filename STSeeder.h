#ifndef STSeeder_h
#define STSeeder_h

#include <string>
#include <vector>
#include <sdsl/suffix_trees.hpp>
#include "GraphAlignerWrapper.h"
#include "GfaGraph.h"

class STSeeder
{
public:
	STSeeder(const GfaGraph& graph);
	std::vector<SeedHit> getSeeds(const std::string& sequence, size_t minMatchSize) const;
private:
	void constructTree(const GfaGraph& graph);
	sdsl::cst_sada<> tree;
	std::vector<size_t> nodePositions;
	std::vector<int> nodeIDs;
};

#endif
