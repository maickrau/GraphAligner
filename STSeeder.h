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
	void addFwMumSeeds(std::vector<SeedHit>& result, const std::string& sequence) const;
	void addBwMumSeeds(std::vector<SeedHit>& result, const std::string& sequence) const;
	void constructTree(const GfaGraph& graph);
	std::vector<SeedHit> removeContainedSeeds(const std::vector<SeedHit>& result) const;
	size_t getNodeIndex(size_t indexPos) const;
	sdsl::cst_sada<> tree;
	std::vector<size_t> nodePositions;
	std::vector<int> nodeIDs;
};

#endif
