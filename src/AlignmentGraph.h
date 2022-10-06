#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <functional>
#include <vector>
#include <tuple>
#include <unordered_set>
#include <phmap.h>
#include "ThreadReadAssertion.h"


class AlignmentGraph
{
public:
	//determines extra band size, shouldn't be too high because of extra slices
	//should be 0 mod (wordsize/2 == 32), otherwise storage has overhead
	//64 is the fastest out of 32, 64, 96
	static constexpr int SPLIT_NODE_SIZE = 64;
	static constexpr size_t BP_IN_CHUNK = sizeof(size_t) * 8 / 2;
	static constexpr size_t CHUNKS_IN_NODE = (SPLIT_NODE_SIZE + BP_IN_CHUNK - 1) / BP_IN_CHUNK;

	struct NodeChunkSequence
	{
		size_t& operator[](size_t pos)
		{
			return s[pos];
		}
		size_t operator[](size_t pos) const
		{
			return s[pos];
		}
		size_t s[CHUNKS_IN_NODE];
	};
	struct AmbiguousChunkSequence
	{
		static_assert(SPLIT_NODE_SIZE == sizeof(size_t)*8);
		//weird interface because it should behave like NodeChunkSequence, which is just a number
		AmbiguousChunkSequence operator[](size_t pos) const
		{
			AmbiguousChunkSequence result = *this;
			result.A >>= pos * BP_IN_CHUNK;
			result.C >>= pos * BP_IN_CHUNK;
			result.G >>= pos * BP_IN_CHUNK;
			result.T >>= pos * BP_IN_CHUNK;
			return result;
		}
		//weird interface because it should behave like NodeChunkSequence, which is just a number
		AmbiguousChunkSequence operator>>=(size_t amount)
		{
			assert(amount % 2 == 0);
			A >>= amount / 2;
			T >>= amount / 2;
			C >>= amount / 2;
			G >>= amount / 2;
			return *this;
		}
		//weird interface because it should behave like NodeChunkSequence, which is just a number
		AmbiguousChunkSequence operator&(size_t val)
		{
			return *this;
		}
		size_t A;
		size_t T;
		size_t C;
		size_t G;
	};

	struct MatrixPosition
	{
		MatrixPosition();
		MatrixPosition(size_t node, size_t nodeOffset, size_t seqPos);
		bool operator==(const MatrixPosition& other) const;
		bool operator!=(const MatrixPosition& other) const;
		size_t node;
		size_t nodeOffset;
		size_t seqPos;
	};

	class SeedHit
	{
	public:
		SeedHit(size_t seqPos, int nodeId, size_t nodePos) : sequencePosition(seqPos), nodeId(nodeId), nodePos(nodePos) {};
		size_t sequencePosition;
		int nodeId;
		size_t nodePos;
	};
	AlignmentGraph();
	void ReserveNodes(size_t numNodes, size_t numSplitNodes);
	void AddNode(int nodeId, const std::string& sequence, const std::string& name, bool reverseNode, const std::vector<size_t>& breakpoints);
	void AddEdgeNodeId(int node_id_from, int node_id_to, size_t startOffset);
	void Finalize(int wordSize);
	std::pair<int, size_t> GetReversePosition(int nodeId, size_t offset) const;
	size_t GetReverseNode(size_t node) const;
	size_t NodeSize() const;
	size_t SizeInBP() const;
	size_t BigraphNodeID(size_t directedNodeId) const;
	size_t NodeLength(size_t nodeIndex) const;
	char NodeSequences(size_t node, size_t offset) const;
	NodeChunkSequence NodeChunks(size_t node) const;
	AmbiguousChunkSequence AmbiguousNodeChunks(size_t node) const;
	size_t GetUnitigNode(int nodeId, size_t offset) const;
	std::string OriginalNodeName(int nodeId) const;
	size_t OriginalNodeSize(int nodeId) const;
	size_t ComponentSize() const;
	static AlignmentGraph DummyGraph();
	std::vector<std::string> originalNodeName;

private:
	void fixChainApproxPos(const size_t start);
	std::pair<bool, size_t> findBubble(const size_t start, const std::vector<bool>& ignorableTip);
	void chainBubble(const size_t start, const std::vector<bool>& ignorableTip, std::vector<size_t>& rank);
	phmap::flat_hash_map<size_t, std::unordered_set<size_t>> chainTips(std::vector<size_t>& rank, std::vector<bool>& ignorableTip);
	void chainCycles(std::vector<size_t>& rank, std::vector<bool>& ignorableTip);
	void findChains();
	void findLinearizable();
	void AddNode(int nodeId, int offset, const std::string& sequence, bool reverseNode);
	void RenumberAmbiguousToEnd();
	void doComponentOrder();
	std::vector<uint8_t> nodeLength;
	std::vector<std::vector<size_t>> nodeLookup;
	std::vector<size_t> originalNodeSize;
	std::vector<size_t> nodeOffset;
	std::vector<int> nodeIDs;
	std::vector<std::vector<size_t>> inNeighbors;
	std::vector<std::vector<size_t>> outNeighbors;
	std::vector<bool> reverse;
	std::vector<bool> linearizable;
	std::vector<NodeChunkSequence> nodeSequences;
	size_t bpSize;
	std::vector<AmbiguousChunkSequence> ambiguousNodeSequences;
	std::vector<bool> ambiguousNodes;
	std::vector<size_t> componentNumber;
	std::vector<size_t> chainNumber;
	std::vector<size_t> chainApproxPos;
	size_t firstAmbiguous;
	bool finalized;

	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAligner;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerVGAlignment;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerGAFAlignment;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerBitvectorBanded;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerBitvectorCommon;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerBitvectorDijkstra;
	friend class DirectedGraph;
	friend class MinimizerSeeder;
};


#endif
