#ifndef MinimizerSeeder_h
#define MinimizerSeeder_h

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
	void initOrder();
	void initEmptyIndex();
	void initMaxCount();
	std::vector<std::vector<std::pair<int, size_t>>> minimizerIndex;
	std::vector<size_t> kmerOrder;
	size_t minimizerLength;
	size_t windowSize;
	size_t maxCount;
};

#endif
