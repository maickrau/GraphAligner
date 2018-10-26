#ifndef STSeeder_h
#define STSeeder_h

#include <string>
#include <vector>
#include <sdsl/suffix_trees.hpp>
#include "GraphAlignerWrapper.h"
#include "GfaGraph.h"
#include "vg.pb.h"

class STSeeder
{
public:
	STSeeder(const GfaGraph& graph);
	STSeeder(const vg::Graph& graph);
	std::vector<SeedHit> getMumSeeds(const std::string& sequence) const;
private:
	void addMumSeeds(std::vector<SeedHit>& result, const std::string& sequence) const;
	void constructTree(const GfaGraph& graph);
	void constructTree(const vg::Graph& graph);
	size_t getNodeIndex(size_t indexPos) const;
	sdsl::cst_sct3<> tree;
	std::vector<size_t> nodePositions;
	std::vector<int> nodeIDs;
	std::vector<bool> nodeReverse;
};

#endif
