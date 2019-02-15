#ifndef MinimizerSeeder_h
#define MinimizerSeeder_h

#include <random>
#include <vector>
#include <string>
#include "GfaGraph.h"
#include "GraphAlignerWrapper.h"
#include "vg.pb.h"

class MinimizerSeeder
{
public:
	MinimizerSeeder(const GfaGraph& graph, size_t minimizerLength, size_t windowSize);
	MinimizerSeeder(const vg::Graph& graph, size_t minimizerLength, size_t windowSize);
	std::vector<SeedHit> getSeeds(const std::string& sequence, size_t maxCount) const;
private:
	SeedHit matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const;
	void addMinimizers(const std::string& str, int nodeId);
	size_t getOrder(size_t kmer) const;
	size_t getOrAddOrder(size_t kmer);
	void initMaxCount();
	void finalizeOrder();
	std::unordered_map<size_t, std::vector<std::pair<int, size_t>>> minimizerIndex;
	std::unordered_map<size_t, size_t> kmerOrder;
	size_t minimizerLength;
	size_t windowSize;
	size_t maxCount;
	std::random_device rd;
	std::mt19937 gen;
	std::uniform_int_distribution<size_t> dis;
};

#endif
