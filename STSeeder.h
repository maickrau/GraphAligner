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
	std::vector<SeedHit> getMumSeeds(const std::string& sequence) const;
private:
	void addMumSeeds(std::vector<SeedHit>& result, const std::string& sequence) const;
	void constructTree(const GfaGraph& graph);
	size_t getNodeIndex(size_t indexPos) const;
	sdsl::cst_sada<> tree;
	std::vector<size_t> nodePositions;
	std::vector<int> nodeIDs;
	std::vector<bool> nodeReverse;
};

#endif
