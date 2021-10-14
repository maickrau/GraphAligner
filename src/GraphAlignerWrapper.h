//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#ifndef GraphAlignerWrapper_h
#define GraphAlignerWrapper_h

#include <tuple>
#include "vg.pb.h"
#include "GraphAlignerCommon.h"
#include "AlignmentGraph.h"

using ReusableStateType = GraphAlignerCommon<size_t, int64_t, uint64_t>::AlignerGraphsizedState;

class SeedHit
{
public:
	SeedHit() :
	nodeID(std::numeric_limits<int>::min()),
	nodeOffset(std::numeric_limits<size_t>::max()),
	seqPos(std::numeric_limits<size_t>::max()),
	matchLen(std::numeric_limits<size_t>::max()),
	reverse(false),
	rawSeedGoodness(0)
	{
	}
	SeedHit(int nodeID, size_t nodeOffset, size_t seqPos, size_t matchLen, size_t rawSeedGoodness, bool reverse) :
	nodeID(nodeID),
	nodeOffset(nodeOffset),
	seqPos(seqPos),
	matchLen(matchLen),
	reverse(reverse),
	rawSeedGoodness(rawSeedGoodness)
	{
	}
	int nodeID;
	size_t nodeOffset;
	size_t seqPos;
	size_t matchLen;
	bool reverse;
	size_t rawSeedGoodness;
};

class ProcessedSeedHit
{
public:
	ProcessedSeedHit(size_t seqPos, size_t alignmentGraphNodeId) :
	seqPos(seqPos),
	alignmentGraphNodeId(alignmentGraphNodeId)
	{
	}
	size_t seqPos;
	size_t alignmentGraphNodeId;
};

class SeedCluster
{
public:
	size_t size() const;
	std::vector<ProcessedSeedHit> hits;
	double clusterGoodness;
};

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t alignmentBandwidth, bool quietMode, ReusableStateType& reusableState, double preciseClippingIdentityCutoff, int Xdropcutoff, size_t DPRestartStride, int clipAmbiguousEnds);
AlignmentResult AlignClusters(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t alignmentBandwidth, size_t maxCellsPerSlice, bool quietMode, const std::vector<SeedCluster>& seedHits, ReusableStateType& reusableState, double preciseClippingIdentityCutoff, int Xdropcutoff, double multimapScoreFraction, int clipAmbiguousEnds, size_t maxTraceCount);

void AddAlignment(const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment);
void AddGAFLine(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment, bool cigarMatchMismatchMerge, bool includeCigar);
void AddCorrected(AlignmentResult::AlignmentItem& alignment);
std::vector<SeedCluster> ClusterSeeds(const AlignmentGraph& graph, const std::vector<SeedHit>& seedHits, const size_t seedClusterMinSize);

#endif