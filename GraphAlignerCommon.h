#ifndef GraphAlignerCommon_h
#define GraphAlignerCommon_h

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerParams
{
public:
	//cutoff for doing the backtrace in the sqrt-slice pass
	//"bulges" in the band are responsible for almost all of the time spent aligning,
	//and this way they don't need to be recalculated, saving about half of the time.
	//semi-arbitrarily hundred thousand, which is about one second of calculation
	static constexpr size_t BacktraceOverrideCutoff = 100000;
	GraphAlignerParams(LengthType initialBandwidth, LengthType rampBandwidth, const AlignmentGraph& graph) :
	initialBandwidth(initialBandwidth),
	rampBandwidth(rampBandwidth),
	graph(graph)
	{
	}
	const LengthType initialBandwidth;
	const LengthType rampBandwidth;
	const AlignmentGraph& graph;
};

#endif
