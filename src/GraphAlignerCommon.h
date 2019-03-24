#ifndef GraphAlignerCommon_h
#define GraphAlignerCommon_h

#include <vector>
#include "AlignmentGraph.h"
#include "ArrayPriorityQueue.h"
#include "ComponentPriorityQueue.h"
#include "NodeSlice.h"
#include "WordSlice.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerCommon
{
public:
	class NodeWithPriority
	{
	public:
		NodeWithPriority(LengthType node, size_t offset, size_t endOffset, int priority) : node(node), offset(offset), endOffset(endOffset), priority(priority) {}
		bool operator>(const NodeWithPriority& other) const
		{
		return priority > other.priority;
		}
		bool operator<(const NodeWithPriority& other) const
		{
		return priority < other.priority;
		}
		LengthType node;
		size_t offset;
		size_t endOffset;
		int priority;
	};
	class EdgeWithPriority
	{
	public:
		EdgeWithPriority(LengthType target, int priority, WordSlice<LengthType, ScoreType, Word> incoming, bool skipFirst) : target(target), priority(priority), incoming(incoming), skipFirst(skipFirst) {}
		bool operator>(const EdgeWithPriority& other) const
		{
			return priority > other.priority;
		}
		bool operator<(const EdgeWithPriority& other) const
		{
			return priority < other.priority;
		}
		LengthType target;
		int priority;
		WordSlice<LengthType, ScoreType, Word> incoming;
		bool skipFirst;
	};
	class AlignerGraphsizedState
	{
	public:
		AlignerGraphsizedState(const AlignmentGraph& graph, size_t maxBandwidth, bool lowMemory) :
		sparseComponentQueue(),
		sparseCalculableQueue(),
		denseComponentQueue(),
		denseCalculableQueue(),
		evenNodesliceMap(),
		oddNodesliceMap(),
		currentBand(),
		previousBand()
		{
			if (!lowMemory)
			{
				evenNodesliceMap.resize(graph.NodeSize(), {});
				oddNodesliceMap.resize(graph.NodeSize(), {});
				denseComponentQueue.initialize(graph.ComponentSize());
				denseCalculableQueue.initialize(WordConfiguration<Word>::WordSize * (WordConfiguration<Word>::WordSize + maxBandwidth + 1) + maxBandwidth + 1, graph.NodeSize());
			}
			else
			{
				sparseComponentQueue.initialize(graph.ComponentSize());
				sparseCalculableQueue.initialize(WordConfiguration<Word>::WordSize * (WordConfiguration<Word>::WordSize + maxBandwidth + 1) + maxBandwidth + 1, graph.NodeSize());
			}
			currentBand.resize(graph.NodeSize(), false);
			previousBand.resize(graph.NodeSize(), false);
		}
		void clear()
		{
			evenNodesliceMap.assign(evenNodesliceMap.size(), {});
			oddNodesliceMap.assign(oddNodesliceMap.size(), {});
			sparseComponentQueue.clear();
			sparseCalculableQueue.clear();
			denseComponentQueue.clear();
			denseCalculableQueue.clear();
			currentBand.assign(currentBand.size(), false);
			previousBand.assign(previousBand.size(), false);
		}
		ComponentPriorityQueue<EdgeWithPriority, true> sparseComponentQueue;
		ArrayPriorityQueue<EdgeWithPriority, true> sparseCalculableQueue;
		ComponentPriorityQueue<EdgeWithPriority, true> denseComponentQueue;
		ArrayPriorityQueue<EdgeWithPriority, true> denseCalculableQueue;
		std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem> evenNodesliceMap;
		std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem> oddNodesliceMap;
		std::vector<bool> currentBand;
		std::vector<bool> previousBand;
	};
	using MatrixPosition = AlignmentGraph::MatrixPosition;
	class Params
	{
	public:
		Params(LengthType initialBandwidth, LengthType rampBandwidth, const AlignmentGraph& graph, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, bool lowMemory, bool forceGlobal, bool preciseClipping) :
		initialBandwidth(initialBandwidth),
		rampBandwidth(rampBandwidth),
		graph(graph),
		maxCellsPerSlice(maxCellsPerSlice),
		quietMode(quietMode),
		sloppyOptimizations(sloppyOptimizations),
		lowMemory(lowMemory),
		forceGlobal(forceGlobal),
		preciseClipping(preciseClipping)
		{
		}
		const LengthType initialBandwidth;
		const LengthType rampBandwidth;
		const AlignmentGraph& graph;
		const size_t maxCellsPerSlice;
		const bool quietMode;
		const bool sloppyOptimizations;
		const bool lowMemory;
		const bool forceGlobal;
		const bool preciseClipping;
	};
	struct TraceItem
	{
		TraceItem() :
		DPposition(),
		nodeSwitch(false),
		sequenceCharacter('-'),
		graphCharacter('-')
		{}
		TraceItem(MatrixPosition DPposition, bool nodeSwitch, char sequenceCharacter, char graphCharacter) :
		DPposition(DPposition),
		nodeSwitch(nodeSwitch),
		sequenceCharacter(sequenceCharacter),
		graphCharacter(graphCharacter)
		{}
		TraceItem(MatrixPosition DPposition, bool nodeSwitch, const std::string& seq, const AlignmentGraph& graph) :
		DPposition(DPposition),
		nodeSwitch(nodeSwitch),
		sequenceCharacter(seq[DPposition.seqPos]),
		graphCharacter(graph.NodeSequences(DPposition.node, DPposition.nodeOffset))
		{}
		MatrixPosition DPposition;
		bool nodeSwitch;
		char sequenceCharacter;
		char graphCharacter;
	};
	class OnewayTrace
	{
	public:
		OnewayTrace() :
		trace(),
		score(0)
		{
		}
		static OnewayTrace TraceFailed()
		{
			OnewayTrace result;
			result.score = std::numeric_limits<ScoreType>::max();
			return result;
		}
		bool failed() const
		{
			return score == std::numeric_limits<ScoreType>::max();
		}
		std::vector<TraceItem> trace;
		ScoreType score;
	};
	class Trace
	{
	public:
		OnewayTrace forward;
		OnewayTrace backward;
	};
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static bool characterMatch(char sequenceCharacter, char graphCharacter)
	{
		return (ambiguousMatch(sequenceCharacter, 'A') && ambiguousMatch(graphCharacter, 'A'))
		|| (ambiguousMatch(sequenceCharacter, 'C') && ambiguousMatch(graphCharacter, 'C'))
		|| (ambiguousMatch(sequenceCharacter, 'G') && ambiguousMatch(graphCharacter, 'G'))
		|| (ambiguousMatch(sequenceCharacter, 'T') && ambiguousMatch(graphCharacter, 'T'));
	}
#ifdef NDEBUG
	__attribute__((always_inline))
#endif
	static bool ambiguousMatch(char ambiguousChar, char exactChar)
	{
		assert(exactChar == 'A' || exactChar == 'T' || exactChar == 'C' || exactChar == 'G');
		switch(ambiguousChar)
		{
			case 'A':
			case 'a':
				return exactChar == 'A';
			break;
			case 'u':
			case 'U':
			case 'T':
			case 't':
				return exactChar == 'T';
			break;
			case 'C':
			case 'c':
				return exactChar == 'C';
			break;
			case 'G':
			case 'g':
				return exactChar == 'G';
			break;
			case 'N':
			case 'n':
				return true;
			break;
			case 'R':
			case 'r':
				return exactChar == 'A' || exactChar == 'G';
			break;
			case 'Y':
			case 'y':
				return exactChar == 'C' || exactChar == 'T';
			break;
			case 'K':
			case 'k':
				return exactChar == 'G' || exactChar == 'T';
			break;
			case 'M':
			case 'm':
				return exactChar == 'C' || exactChar == 'A';
			break;
			case 'S':
			case 's':
				return exactChar == 'C' || exactChar == 'G';
			break;
			case 'W':
			case 'w':
				return exactChar == 'A' || exactChar == 'T';
			break;
			case 'B':
			case 'b':
				return exactChar == 'C' || exactChar == 'G' || exactChar == 'T';
			break;
			case 'D':
			case 'd':
				return exactChar == 'A' || exactChar == 'G' || exactChar == 'T';
			break;
			case 'H':
			case 'h':
				return exactChar == 'A' || exactChar == 'C' || exactChar == 'T';
			break;
			case 'V':
			case 'v':
				return exactChar == 'A' || exactChar == 'C' || exactChar == 'G';
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
