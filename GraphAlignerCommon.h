#ifndef GraphAlignerCommon_h
#define GraphAlignerCommon_h

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerParams
{
public:
	//band size in bp when the alternate method is used instead of the bitvector method
	//empirically, two hundred thousand is (close to) the fastest cutoff for aligning ONT's to human DBG
	static constexpr size_t AlternateMethodCutoff = 2000000000;
	//cutoff for doing the backtrace in the sqrt-slice pass
	//"bulges" in the band are responsible for almost all of the time spent aligning,
	//and this way they don't need to be recalculated, saving about half of the time.
	//must be the same as AlternateMethodCutoff because of cell existance etc.
	static constexpr size_t BacktraceOverrideCutoff = AlternateMethodCutoff;
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
