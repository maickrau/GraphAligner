//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#ifndef GraphAlignerWrapper_h
#define GraphAlignerWrapper_h

#include <tuple>
#include "AlignmentGraph.h"
#include "vg.pb.h"

class AlignmentResult
{
public:
	AlignmentResult()
	{
	}
	AlignmentResult(vg::Alignment alignment, bool alignmentFailed, size_t cellsProcessed, size_t ms) :
	alignment(alignment),
	alignmentFailed(alignmentFailed),
	cellsProcessed(cellsProcessed),
	elapsedMilliseconds(ms),
	alignmentStart(0),
	alignmentEnd(0)
	{
	}
	vg::Alignment alignment;
	bool alignmentFailed;
	size_t cellsProcessed;
	size_t elapsedMilliseconds;
	size_t alignmentStart;
	size_t alignmentEnd;
};

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, size_t dynamicRowStart, bool sqrtSpace);
AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, size_t dynamicRowStart, bool sqrtSpace, const std::vector<std::tuple<int, size_t, bool>>& seedHits);

#endif