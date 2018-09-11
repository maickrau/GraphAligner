//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include <limits>
#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, bool quietMode, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {initialBandwidth, rampBandwidth, graph, std::numeric_limits<size_t>::max(), quietMode, false, lowMemory};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWay(seq_id, sequence, reusableState);
}

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<std::tuple<int, size_t, size_t, size_t, bool>>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {initialBandwidth, rampBandwidth, graph, maxCellsPerSlice, quietMode, sloppyOptimizations, lowMemory};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	std::vector<GraphAlignerCommon<size_t, int32_t, uint64_t>::SeedHit> seeds;
	seeds.reserve(seedHits.size());
	for (auto seed : seedHits)
	{
		seeds.emplace_back(std::get<0>(seed), std::get<1>(seed), std::get<2>(seed), std::get<3>(seed), std::get<4>(seed));
	}
	return aligner.AlignOneWaySubgraph(seq_id, sequence, seeds, reusableState);
	return aligner.AlignOneWay(seq_id, sequence, seeds, reusableState);
}
