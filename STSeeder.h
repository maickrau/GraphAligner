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
	STSeeder(const GfaGraph& graph, const std::string& cachePrefix);
	STSeeder(const vg::Graph& graph, const std::string& cachePrefix);
	std::vector<SeedHit> getMumSeeds(const std::string& sequence, size_t maxCount) const;
private:
	void addMumSeeds(std::vector<SeedHit>& result, const std::string& sequence) const;
	void constructTree(const GfaGraph& graph);
	void constructTree(const vg::Graph& graph);
	size_t getNodeIndex(size_t indexPos) const;
	void saveTo(const std::string& cachePrefix) const;
	void loadFrom(const std::string& cachePrefix);
	sdsl::cst_sct3<> tree;
	std::vector<size_t> nodePositions;
	std::vector<int> nodeIDs;
	std::vector<bool> nodeReverse;
};

#endif
