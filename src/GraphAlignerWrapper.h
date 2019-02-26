//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#ifndef GraphAlignerWrapper_h
#define GraphAlignerWrapper_h

#include <tuple>
#include "GraphAlignerCommon.h"
#include "AlignmentGraph.h"
#include "vg.pb.h"

class AlignmentResult
{
public:
	AlignmentResult() :
		alignments(),
		seedsExtended(0)
	{}
	enum TraceMatchType
	{
		//relative to the graph, aka insertion has no graphchar, but has readchar
		MATCH = 1,
		MISMATCH = 2,
		INSERTION = 3,
		DELETION = 4,
		FORWARDBACKWARDSPLIT = 5
	};
	struct TraceItem
	{
		int nodeID;
		size_t offset;
		bool reverse;
		size_t readpos;
		TraceMatchType type;
		char graphChar;
		char readChar;
	};
	class AlignmentItem
	{
	public:
		AlignmentItem() :
		cellsProcessed(0),
		elapsedMilliseconds(0),
		alignmentStart(0),
		alignmentEnd(0)
		{}
		AlignmentItem(std::shared_ptr<vg::Alignment> alignment, size_t cellsProcessed, size_t ms) :
		alignment(alignment),
		cellsProcessed(cellsProcessed),
		elapsedMilliseconds(ms),
		alignmentStart(0),
		alignmentEnd(0)
		{}
		bool alignmentFailed() const
		{
			return alignmentEnd == alignmentStart;
		}
		std::shared_ptr<vg::Alignment> alignment;
		std::vector<TraceItem> trace;
		size_t cellsProcessed;
		size_t elapsedMilliseconds;
		size_t alignmentStart;
		size_t alignmentEnd;
	};
	std::vector<AlignmentItem> alignments;
	size_t seedsExtended;
};

class SeedHit
{
public:
	SeedHit(int nodeID, size_t nodeOffset, size_t seqPos, size_t matchLen, bool reverse) :
	nodeID(nodeID),
	nodeOffset(nodeOffset),
	seqPos(seqPos),
	matchLen(matchLen),
	reverse(reverse)
	{
	}
	int nodeID;
	size_t nodeOffset;
	size_t seqPos;
	size_t matchLen;
	bool reverse;
};

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t initialBandwidth, size_t rampBandwidth, bool quietMode, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory, bool forceGlobal, bool preciseClipping);
AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, size_t initialBandwidth, size_t rampBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<SeedHit>& seedHits, GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState& reusableState, bool lowMemory, bool forceGlobal, bool preciseClipping);

#endif