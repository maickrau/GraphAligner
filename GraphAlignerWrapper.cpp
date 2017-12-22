//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, size_t dynamicRowStart)
{
	GraphAligner<size_t, int32_t, uint64_t> aligner {graph, initialBandwidth, rampBandwidth};
	return aligner.AlignOneWay(seq_id, sequence, dynamicRowStart);
}

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, size_t dynamicRowStart, const std::vector<std::tuple<int, size_t, bool>>& seedHits)
{
	GraphAligner<size_t, int32_t, uint64_t> aligner {graph, initialBandwidth, rampBandwidth};
	return aligner.AlignOneWay(seq_id, sequence, dynamicRowStart, seedHits);
}
