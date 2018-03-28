#ifndef GraphAlignerCommon_h
#define GraphAlignerCommon_h

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerParams
{
public:
	//cutoff for doing the backtrace in the sqrt-slice pass
	//"bulges" in the band are responsible for almost all of the time spent aligning,
	//and this way they don't need to be recalculated, saving about half of the time.
	//semi-arbitrarily fifty thousand, empirically a good enough cutoff
	static constexpr size_t BacktraceOverrideCutoff = 50000;
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

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerCommon
{
public:
	static bool characterMatch(char sequenceCharacter, char graphCharacter)
	{
		assert(graphCharacter == 'A' || graphCharacter == 'T' || graphCharacter == 'C' || graphCharacter == 'G');
		switch(sequenceCharacter)
		{
			case 'A':
			case 'a':
			return graphCharacter == 'A';
			break;
			case 'T':
			case 't':
			return graphCharacter == 'T';
			break;
			case 'C':
			case 'c':
			return graphCharacter == 'C';
			break;
			case 'G':
			case 'g':
			return graphCharacter == 'G';
			break;
			case 'N':
			case 'n':
			return true;
			break;
			case 'R':
			case 'r':
			return graphCharacter == 'A' || graphCharacter == 'G';
			break;
			case 'Y':
			case 'y':
			return graphCharacter == 'C' || graphCharacter == 'T';
			break;
			case 'K':
			case 'k':
			return graphCharacter == 'G' || graphCharacter == 'T';
			break;
			case 'M':
			case 'm':
			return graphCharacter == 'C' || graphCharacter == 'A';
			break;
			case 'S':
			case 's':
			return graphCharacter == 'C' || graphCharacter == 'G';
			break;
			case 'W':
			case 'w':
			return graphCharacter == 'A' || graphCharacter == 'T';
			break;
			case 'B':
			case 'b':
			return graphCharacter == 'C' || graphCharacter == 'G' || graphCharacter == 'T';
			break;
			case 'D':
			case 'd':
			return graphCharacter == 'A' || graphCharacter == 'G' || graphCharacter == 'T';
			break;
			case 'H':
			case 'h':
			return graphCharacter == 'A' || graphCharacter == 'C' || graphCharacter == 'T';
			break;
			case 'V':
			case 'v':
			return graphCharacter == 'A' || graphCharacter == 'C' || graphCharacter == 'G';
			break;
			default:
			assert(false);
			std::abort();
			return false;
			break;
		}
	}
};

#endif
