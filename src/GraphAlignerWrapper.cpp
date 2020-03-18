//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include <limits>
#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t initialBandwidth, size_t rampBandwidth, bool quietMode, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory, bool forceGlobal, bool preciseClipping)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {initialBandwidth, rampBandwidth, graph, std::numeric_limits<size_t>::max(), quietMode, false, lowMemory, forceGlobal, preciseClipping, 1};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWay(seq_id, sequence, reusableState);
}

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t initialBandwidth, size_t rampBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<SeedHit>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory, bool forceGlobal, bool preciseClipping, size_t minClusterSize)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {initialBandwidth, rampBandwidth, graph, maxCellsPerSlice, quietMode, sloppyOptimizations, lowMemory, forceGlobal, preciseClipping, minClusterSize};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWay(seq_id, sequence, seedHits, reusableState);
}

void AddAlignment(const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, AlignmentGraph::DummyGraph(), 1, true, true, true, false, false, 1};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddAlignment(seq_id, sequence, alignment);
}

void AddGAFLine(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, graph, 1, true, true, true, false, false, 1};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddGAFLine(seq_id, sequence, alignment);
}

void AddCorrected(AlignmentResult::AlignmentItem& alignment)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, AlignmentGraph::DummyGraph(), 1, true, true, true, false, false, 1};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddCorrected(alignment);
}

void OrderSeeds(const AlignmentGraph& graph, std::vector<SeedHit>& seedHits)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, graph, 1, true, true, true, false, false, 1};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.orderSeedsByChaining(seedHits);
}
