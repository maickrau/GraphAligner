#ifndef MinimizerSeeder_h
#define MinimizerSeeder_h

#include <random>
#include <vector>
#include <string>
#include <sdsl/int_vector.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <ParallelBB.h>
#include "AlignmentGraph.h"
#include "GraphAlignerWrapper.h"

class MinimizerSeeder
{
public:
	MinimizerSeeder(const AlignmentGraph& graph, size_t minimizerLength, size_t windowSize, size_t numThreads);
	std::vector<SeedHit> getSeeds(const std::string& sequence, size_t maxCount, size_t chunkSize) const;
private:
	size_t getStart(size_t index) const;
	void initMinimizers(size_t numThreads);
	SeedHit matchToSeedHit(int nodeId, size_t nodeOffset, size_t seqPos, int count) const;
	void initMaxCount();
	const AlignmentGraph& graph;
	ParallelBB<size_t> locator;
	sdsl::int_vector<0> kmerCheck;
	sdsl::bit_vector starts;
	sdsl::bit_vector::select_1_type startSelector;
	sdsl::int_vector<0> positions;
	size_t minimizerLength;
	size_t windowSize;
	size_t maxCount;
};

#endif
