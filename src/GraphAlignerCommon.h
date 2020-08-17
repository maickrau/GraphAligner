#ifndef GraphAlignerCommon_h
#define GraphAlignerCommon_h

#include <vector>
#include "AlignmentGraph.h"
#include "ArrayPriorityQueue.h"
#include "ComponentPriorityQueue.h"
#include "NodeSlice.h"
#include "WordSlice.h"
#include "DijkstraQueue.h"

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
		EdgeWithPriority(LengthType target, int priority, WordSlice<LengthType, ScoreType, Word> incoming, bool skipFirst) : target(target), priority(priority), incoming(incoming), skipFirst(skipFirst), slice(0) {}
		EdgeWithPriority(LengthType target, int priority, WordSlice<LengthType, ScoreType, Word> incoming, bool skipFirst, size_t slice) : target(target), priority(priority), incoming(incoming), skipFirst(skipFirst), slice(slice) {}
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
		size_t slice;
	};
	class AlignerGraphsizedState
	{
	public:
		AlignerGraphsizedState(const AlignmentGraph& graph, size_t maxBandwidth, bool lowMemory) :
		componentQueue(),
		calculableQueue(),
		evenNodesliceMap(),
		oddNodesliceMap(),
		currentBand(),
		previousBand()
		{
			if (!lowMemory)
			{
				evenNodesliceMap.resize(graph.NodeSize(), {});
				oddNodesliceMap.resize(graph.NodeSize(), {});
			}
			componentQueue.initialize(graph.ComponentSize());
			calculableQueue.initialize(WordConfiguration<Word>::WordSize * (WordConfiguration<Word>::WordSize + maxBandwidth + 1) + maxBandwidth + 1, graph.NodeSize());
			currentBand.resize(graph.NodeSize(), false);
			previousBand.resize(graph.NodeSize(), false);
		}
		void clear()
		{
			evenNodesliceMap.assign(evenNodesliceMap.size(), {});
			oddNodesliceMap.assign(oddNodesliceMap.size(), {});
			componentQueue.clear();
			calculableQueue.clear();
			dijkstraQueue.clear();
			currentBand.assign(currentBand.size(), false);
			previousBand.assign(previousBand.size(), false);
		}
		ComponentPriorityQueue<EdgeWithPriority, true> componentQueue;
		ArrayPriorityQueue<EdgeWithPriority, true> calculableQueue;
		DijkstraPriorityQueue<EdgeWithPriority> dijkstraQueue;
		std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem> evenNodesliceMap;
		std::vector<typename NodeSlice<LengthType, ScoreType, Word, true>::MapItem> oddNodesliceMap;
		std::vector<bool> currentBand;
		std::vector<bool> previousBand;
	};
	using MatrixPosition = AlignmentGraph::MatrixPosition;
	class Params
	{
	public:
		Params(LengthType initialBandwidth, LengthType rampBandwidth, const AlignmentGraph& graph, size_t maxCellsPerSlice, bool quietMode, bool sloppyOptimizations, bool lowMemory, bool forceGlobal, bool preciseClipping, size_t minSeedClusterSize, double seedExtendDensity, bool nondeterministicOptimizations, double preciseClippingIdentityCutoff, int Xdropcutoff) :
		initialBandwidth(initialBandwidth),
		rampBandwidth(rampBandwidth),
		graph(graph),
		maxCellsPerSlice(maxCellsPerSlice),
		quietMode(quietMode),
		sloppyOptimizations(sloppyOptimizations),
		lowMemory(lowMemory),
		forceGlobal(forceGlobal),
		preciseClipping(preciseClipping),
		minSeedClusterSize(minSeedClusterSize),
		seedExtendDensity(seedExtendDensity),
		nondeterministicOptimizations(nondeterministicOptimizations),
		XscoreErrorCost(preciseClippingIdentityCutoff / (1.0 - preciseClippingIdentityCutoff) + 1.0),
		Xdropcutoff(Xdropcutoff)
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
		const size_t minSeedClusterSize;
		const double seedExtendDensity;
		const bool nondeterministicOptimizations;
		const double XscoreErrorCost;
		const int Xdropcutoff;
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
		sequenceCharacter(DPposition.seqPos < seq.size() ? seq[DPposition.seqPos] : '-'),
		graphCharacter(graph.NodeSequences(DPposition.node, DPposition.nodeOffset))
		{}
		TraceItem(MatrixPosition DPposition, bool nodeSwitch, const std::string_view& seq, const AlignmentGraph& graph) :
		DPposition(DPposition),
		nodeSwitch(nodeSwitch),
		sequenceCharacter(DPposition.seqPos < seq.size() ? seq[DPposition.seqPos] : '-'),
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
		// force move semantics because copying is very slow and unnecessary
		OnewayTrace(const OnewayTrace& other) = delete;
		OnewayTrace(OnewayTrace&& other) = default;
		OnewayTrace& operator=(const OnewayTrace& other) = delete;
		OnewayTrace& operator=(OnewayTrace&& other) = default;
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
		if (sequenceCharacter == graphCharacter) return true;
		switch(sequenceCharacter)
		{
			case 'a':
			case 'A':
				return ambiguousMatch(graphCharacter, 'A');
			case 'c':
			case 'C':
				return ambiguousMatch(graphCharacter, 'C');
			case 'g':
			case 'G':
				return ambiguousMatch(graphCharacter, 'G');
			case 't':
			case 'T':
				return ambiguousMatch(graphCharacter, 'T');
			case '-':
				return false;
		}
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
			case '-':
				return false;
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

class AlignmentResult
{
public:
	AlignmentResult() :
		alignments(),
		seedsExtended(0)
	{}
	class AlignmentItem
	{
	public:
		AlignmentItem() :
		corrected(),
		alignment(),
		trace(),
		seedGoodness(0),
		cellsProcessed(0),
		elapsedMilliseconds(0),
		alignmentStart(0),
		alignmentEnd(0),
		alignmentScore(std::numeric_limits<size_t>::max())
		{}
		AlignmentItem(GraphAlignerCommon<size_t, int32_t, uint64_t>::OnewayTrace&& trace, size_t cellsProcessed, size_t ms) :
		corrected(),
		alignment(),
		trace(),
		cellsProcessed(cellsProcessed),
		elapsedMilliseconds(ms),
		alignmentStart(0),
		alignmentEnd(0),
		alignmentScore(std::numeric_limits<size_t>::max())
		{
			this->trace = std::make_shared<GraphAlignerCommon<size_t, int32_t, uint64_t>::OnewayTrace>();
			*this->trace = std::move(trace);
		}
		bool alignmentFailed() const
		{
			return alignmentEnd == alignmentStart;
		}
		std::string corrected;
		std::string GAFline;
		std::shared_ptr<vg::Alignment> alignment;
		std::shared_ptr<GraphAlignerCommon<size_t, int32_t, uint64_t>::OnewayTrace> trace;
		size_t seedGoodness;
		size_t cellsProcessed;
		size_t elapsedMilliseconds;
		size_t alignmentStart;
		size_t alignmentEnd;
		size_t alignmentScore;
	};
	std::vector<AlignmentItem> alignments;
	size_t seedsExtended;
	std::string readName;
};

#endif
