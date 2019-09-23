#ifndef MinimizerSeeder_h
#define MinimizerSeeder_h

#include <random>
#include <vector>
#include <string>
#include <phmap.h>
#include "GfaGraph.h"
#include "GraphAlignerWrapper.h"
#include "vg.pb.h"

class MinimizerSeeder
{
public:
	MinimizerSeeder(const GfaGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads);
	MinimizerSeeder(const vg::Graph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads);
	std::vector<SeedHit> getSeeds(const std::string& sequence, size_t maxCount, size_t chunkSize) const;
private:
	void initMinimizers(const GfaGraph& graph, size_t numThreads);
	void initMinimizers(const vg::Graph& graph, size_t numThreads);
	SeedHit matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const;
	void initMaxCount();
	phmap::flat_hash_map<size_t, size_t> locator;
	std::vector<size_t> starts;
	std::vector<std::pair<int, size_t>> positions;
	size_t minimizerLength;
	size_t windowSize;
	size_t maxCount;
	std::random_device rd;
	std::mt19937 gen;
	std::uniform_int_distribution<size_t> dis;
};

#endif
