//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include <limits>
#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

size_t SeedCluster::size() const
{
	return hits.size();
}

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t alignmentBandwidth, bool quietMode, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, double preciseClippingIdentityCutoff, int Xdropcutoff, size_t DPRestartStride)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {alignmentBandwidth, graph, std::numeric_limits<size_t>::max(), quietMode, preciseClippingIdentityCutoff, Xdropcutoff, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWay(seq_id, sequence, reusableState, DPRestartStride);
}

AlignmentResult AlignClusters(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t alignmentBandwidth, size_t maxCellsPerSlice, bool quietMode, const std::vector<SeedCluster>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, double preciseClippingIdentityCutoff, int Xdropcutoff, double multimapScoreFraction)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {alignmentBandwidth, graph, maxCellsPerSlice, quietMode, preciseClippingIdentityCutoff, Xdropcutoff, multimapScoreFraction};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignClusters(seq_id, sequence, seedHits, reusableState);
}

void AddAlignment(const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, AlignmentGraph::DummyGraph(), 1, true, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddAlignment(seq_id, sequence, alignment);
}

void AddGAFLine(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment, bool cigarMatchMismatchMerge, bool includeCigar)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, graph, 1, true, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddGAFLine(seq_id, sequence, alignment, cigarMatchMismatchMerge, includeCigar);
}

void AddCorrected(AlignmentResult::AlignmentItem& alignment)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, AlignmentGraph::DummyGraph(), 1, true, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	aligner.AddCorrected(alignment);
}

std::vector<SeedCluster> ClusterSeeds(const AlignmentGraph& graph, const std::vector<SeedHit>& seedHits, const size_t seedClusterMinSize)
{
	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params params {1, graph, 1, true, .5, 0, 0};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.clusterSeeds(seedHits, seedClusterMinSize);
}
