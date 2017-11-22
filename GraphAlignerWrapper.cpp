//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int dynamicWidth, size_t dynamicRowStart, bool sqrtSpace)
{
	GraphAligner<size_t, int32_t, uint64_t> aligner {graph};
	return aligner.AlignOneWay(seq_id, sequence, dynamicWidth, dynamicRowStart, sqrtSpace);
}

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int dynamicWidth, size_t dynamicRowStart, bool sqrtSpace, const std::vector<std::tuple<int, size_t, bool>>& seedHits, int startBandwidth)
{
	GraphAligner<size_t, int32_t, uint64_t> aligner {graph};
	return aligner.AlignOneWay(seq_id, sequence, dynamicWidth, dynamicRowStart, sqrtSpace, seedHits, startBandwidth);
}

AlignmentResult CollectStats(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int dynamicWidth, size_t dynamicRowStart, bool sqrtSpace, const std::vector<std::tuple<int, size_t, bool>>& seedHits, int startBandwidth)
{
	GraphAligner<size_t, int32_t, uint64_t> aligner {graph};
	return aligner.GetAlignmentStats(seq_id, sequence, dynamicWidth, dynamicRowStart, sqrtSpace, seedHits, startBandwidth);
}
