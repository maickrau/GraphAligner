//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#ifndef GraphAlignerWrapper_h
#define GraphAlignerWrapper_h

#include <tuple>
#include "vg.pb.h"
#include "GraphAlignerCommon.h"
#include "AlignmentGraph.h"

class SeedHit
{
public:
	SeedHit() :
	nodeID(std::numeric_limits<int>::min()),
	nodeOffset(std::numeric_limits<size_t>::max()),
	seqPos(std::numeric_limits<size_t>::max()),
	matchLen(std::numeric_limits<size_t>::max()),
	reverse(false),
	alignmentGraphNodeId(std::numeric_limits<size_t>::max()),
	alignmentGraphNodeOffset(std::numeric_limits<size_t>::max()),
	rawSeedGoodness(0),
	seedGoodness(0),
	seedClusterSize(0)
	{
	}
	SeedHit(int nodeID, size_t nodeOffset, size_t seqPos, size_t matchLen, size_t rawSeedGoodness, bool reverse) :
	nodeID(nodeID),
	nodeOffset(nodeOffset),
	seqPos(seqPos),
	matchLen(matchLen),
	reverse(reverse),
	alignmentGraphNodeId(std::numeric_limits<size_t>::max()),
	alignmentGraphNodeOffset(std::numeric_limits<size_t>::max()),
	rawSeedGoodness(rawSeedGoodness),
	seedGoodness(0),
	seedClusterSize(0)
	{
	}
	int nodeID;
	size_t nodeOffset;
	size_t seqPos;
	size_t matchLen;
	bool reverse;
	size_t alignmentGraphNodeId;
	size_t alignmentGraphNodeOffset;
	size_t rawSeedGoodness;
	size_t seedGoodness;
	size_t seedClusterSize;
};

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t alignmentBandwidth, bool quietMode, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, double preciseClippingIdentityCutoff, int Xdropcutoff, size_t DPRestartStride);
AlignmentResult AlignMultiseed(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t alignmentBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<SeedHit>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, size_t minClusterSize, double seedExtendDensity, double preciseClippingIdentityCutoff, int Xdropcutoff, double multimapScoreFraction);
AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t alignmentBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<SeedHit>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, size_t minClusterSize, double seedExtendDensity, double preciseClippingIdentityCutoff, int Xdropcutoff);

void AddAlignment(const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment);
void AddGAFLine(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, AlignmentResult::AlignmentItem& alignment, bool cigarMatchMismatchMerge);
void AddCorrected(AlignmentResult::AlignmentItem& alignment);
void OrderSeeds(const AlignmentGraph& graph, std::vector<SeedHit>& seedHits);
void PrepareMultiseeds(const AlignmentGraph& graph, std::vector<SeedHit>& seedHits, const size_t seqLen);

#endif