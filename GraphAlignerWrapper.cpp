//split this here so modifying GraphAligner.h doesn't require recompiling every cpp file

#include "GraphAlignerWrapper.h"
#include "GraphAligner.h"
#include "ThreadReadAssertion.h"

// AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth)
// {
// 	GraphAlignerParams<size_t, int32_t, uint64_t> params {initialBandwidth, rampBandwidth, graph};
// 	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
// 	return aligner.AlignOneWay(seq_id, sequence);
// }

AlignmentResult AlignOneWay(const AlignmentGraph& graph, const std::string& seq_id, const std::string& sequence, int initialBandwidth, int rampBandwidth, const std::vector<std::tuple<int, size_t, bool>>& seedHits)
{
	GraphAlignerParams<size_t, int32_t, uint64_t> params {initialBandwidth, rampBandwidth, graph};
	GraphAligner<size_t, int32_t, uint64_t> aligner {params};
	return aligner.AlignOneWay(seq_id, sequence, seedHits);
}
