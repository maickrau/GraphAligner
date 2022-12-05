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
	MEMSeeder(const GfaGraph& graph, const std::string& cachePrefix, double uniqueBonusFactor, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree, const size_t windowSize);
	MEMSeeder(const vg::Graph& graph, const std::string& cachePrefix, double uniqueBonusFactor, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree, const size_t windowSize);
	std::vector<SeedHit> getMemSeeds(const std::string& sequence, size_t maxCount, size_t minLen) const;
	std::vector<SeedHit> getMumSeeds(const std::string& sequence, size_t maxCount, size_t minLen) const;
private:
	void saveTo(const std::string& filename) const;
	void loadFrom(const std::string& filename);
	SeedHit matchToSeed(MEMfinder::Match match) const;
	size_t getNodeIndex(size_t indexPos) const;
	size_t nodeLength(size_t indexPos) const;
	void initTree(const GfaGraph& graph, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree);
	void initTree(const vg::Graph& graph, const bool lowMemoryMEMIndexConstruction, const bool useWaveletTree);
	FMIndex index;
	std::vector<std::pair<size_t, size_t>> prefixIndex;
	std::vector<uint64_t> nodePositions;
	std::vector<uint64_t> nodeIDs;
	const double uniqueBonusFactor;
	size_t windowSize;
};

#endif
