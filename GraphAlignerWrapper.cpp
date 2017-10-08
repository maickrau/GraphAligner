//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int dynamicWidth, size_t dynamicRowStart)
{
	GraphAligner<size_t, int32_t, uint64_t> aligner {graph};
	return aligner.AlignOneWay(seq_id, sequence, dynamicWidth, dynamicRowStart);
}

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int dynamicWidth, size_t dynamicRowStart, const std::vector<std::pair<int, size_t>>& seedHits, int startBandwidth)
{
	GraphAligner<size_t, int32_t, uint64_t> aligner {graph};
	return aligner.AlignOneWay(seq_id, sequence, dynamicWidth, dynamicRowStart, seedHits, startBandwidth);
}
