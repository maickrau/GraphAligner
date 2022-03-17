#ifndef MEMSeeder_h
#define MEMSeeder_h

#include <vector>
#include <string>
#include "MEMfinder.h"
#include "FMIndex.h"
#include "GfaGraph.h"
#include "GraphAlignerWrapper.h"
#include "vg.pb.h"

class MEMSeeder
{
public:
	MEMSeeder(const GfaGraph& graph, const std::string& cachePrefix);
	MEMSeeder(const vg::Graph& graph, const std::string& cachePrefix);
	std::vector<SeedHit> getMemSeeds(const std::string& sequence, size_t maxCount, size_t minLen) const;
	std::vector<SeedHit> getMumSeeds(const std::string& sequence, size_t maxCount, size_t minLen) const;
private:
	void saveTo(const std::string& filename) const;
	void loadFrom(const std::string& filename);
	SeedHit matchToSeed(MEMfinder::Match match) const;
	size_t getNodeIndex(size_t indexPos) const;
	size_t nodeLength(size_t indexPos) const;
	void initTree(const GfaGraph& graph);
	void initTree(const vg::Graph& graph);
	FMIndex index;
	std::vector<size_t> nodePositions;
	std::vector<size_t> nodeIDs;
};

#endif
