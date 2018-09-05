//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#ifndef GraphAlignerWrapper_h
#define GraphAlignerWrapper_h

#include <tuple>
#include "AlignmentGraph.h"

class AlignmentResult
{
public:
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
		AlignmentItem(vg::Alignment alignment, size_t cellsProcessed, size_t ms) :
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
		vg::Alignment alignment;
		std::vector<TraceItem> trace;
		size_t cellsProcessed;
		size_t elapsedMilliseconds;
		size_t alignmentStart;
		size_t alignmentEnd;
	};
	std::vector<AlignmentItem> alignments;
};

// AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth);
AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, const std::vector<std::tuple<int, size_t, bool>>& seedHits);

#endif