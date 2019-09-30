#ifndef MinimizerSeeder_h
#define MinimizerSeeder_h

#include <random>
#include <vector>
#include <string>
#include <MinimalHashmap.h>
#include "AlignmentGraph.h"
#include "GraphAlignerWrapper.h"

class MinimizerSeeder
{
public:
	MinimizerSeeder(const AlignmentGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads);
	std::vector<SeedHit> getSeeds(const std::string& sequence, size_t maxCount, size_t chunkSize) const;
private:
	void initMinimizers(size_t numThreads);
	SeedHit matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const;
	void initMaxCount();
	const AlignmentGraph& graph;
	MinimalHashmap<size_t, size_t> locator;
	std::vector<size_t> starts;
	std::vector<std::pair<int, size_t>> positions;
	size_t minimizerLength;
	size_t windowSize;
	size_t maxCount;
};

#endif
