//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include <limits>
#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t initialBandwidth, size_t rampBandwidth, bool quietMode, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory, bool forceGlobal, bool nondeterministicOptimizations, double preciseClippingIdentityCutoff, int Xdropcutoff, size_t DPRestartStride)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {initialBandwidth, rampBandwidth, graph, std::numeric_limits<size_t>::max(), quietMode, false, lowMemory, forceGlobal, 1, 0, nondeterministicOptimizations, preciseClippingIdentityCutoff, Xdropcutoff, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWay(seq_id, sequence, reusableState, DPRestartStride);
}

AlignmentResult AlignOneWayDijkstra(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, bool quietMode, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool forceGlobal)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, graph, std::numeric_limits<size_t>::max(), quietMode, false, true, forceGlobal, 1, 0, false, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWayDijkstra(seq_id, sequence, reusableState);
}

AlignmentResult AlignMultiseed(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t initialBandwidth, size_t rampBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<SeedHit>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory, bool forceGlobal, size_t minClusterSize, double seedExtendDensity, bool nondeterministicOptimizations, double preciseClippingIdentityCutoff, int Xdropcutoff, double multimapScoreFraction)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {initialBandwidth, rampBandwidth, graph, maxCellsPerSlice, quietMode, sloppyOptimizations, lowMemory, forceGlobal, minClusterSize, seedExtendDensity, nondeterministicOptimizations, preciseClippingIdentityCutoff, Xdropcutoff, multimapScoreFraction};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignMultiseed(seq_id, sequence, seedHits, reusableState);
}

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t initialBandwidth, size_t rampBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<SeedHit>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory, bool forceGlobal, size_t minClusterSize, double seedExtendDensity, bool nondeterministicOptimizations, double preciseClippingIdentityCutoff, int Xdropcutoff)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {initialBandwidth, rampBandwidth, graph, maxCellsPerSlice, quietMode, sloppyOptimizations, lowMemory, forceGlobal, minClusterSize, seedExtendDensity, nondeterministicOptimizations, preciseClippingIdentityCutoff, Xdropcutoff, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWay(seq_id, sequence, seedHits, reusableState);
}

void AddAlignment(const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, AlignmentGraph::DummyGraph(), 1, true, true, true, false, 1, 0, false, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddAlignment(seq_id, sequence, alignment);
}

void AddGAFLine(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment, bool cigarMatchMismatchMerge)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, graph, 1, true, true, true, false, 1, 0, false, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddGAFLine(seq_id, sequence, alignment, cigarMatchMismatchMerge);
}

void AddCorrected(AlignmentResult::AlignmentItem& alignment)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, AlignmentGraph::DummyGraph(), 1, true, true, true, false, 1, 0, false, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddCorrected(alignment);
}

void OrderSeeds(const AlignmentGraph& graph, std::vector<SeedHit>& seedHits)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, graph, 1, true, true, true, false, 1, 0, false, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.orderSeedsByChaining(seedHits);
}

void PrepareMultiseeds(const AlignmentGraph& graph, std::vector<SeedHit>& seedHits, const size_t seqLen)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, 1, graph, 1, true, true, true, false, 1, 0, false, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	seedHits = aligner.prepareSeedsForMultiseeding(seedHits, seqLen);
}
