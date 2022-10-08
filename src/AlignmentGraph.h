#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <iterator>
#include <functional>
#include <vector>
#include <tuple>
#include <unordered_set>
#include <set>
#include <phmap.h>
#include "ThreadReadAssertion.h"
#include "DNAString.h"

class AlignmentGraph
{
public:
	class EdgeIterator : public std::iterator<std::input_iterator_tag, size_t>
	{
	public:
		EdgeIterator(size_t implicitEdge, const size_t* vecPointer);
		EdgeIterator& operator=(const EdgeIterator& other) = default;
		bool operator!=(const EdgeIterator& other) const;
		bool operator==(const EdgeIterator& other) const;
		size_t operator*() const;
		EdgeIterator& operator++();
		EdgeIterator operator++(int);
	private:
		size_t implicitEdge;
		const size_t* vecPointer;
	};
	class NodeEdgeIterator
	{
	public:
		using iterator = EdgeIterator;
		NodeEdgeIterator(size_t implicitEdge, const size_t* startPointer, const size_t* endPointer);
		NodeEdgeIterator& operator=(const NodeEdgeIterator& other) = default;
		EdgeIterator begin() const;
		EdgeIterator end() const;
		size_t size() const;
		size_t operator[](size_t index) const;
	private:
		size_t implicitEdge;
		const size_t* startPointer;
		const size_t* endPointer;
	};
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
	void AddNode(int nodeId, const DNAString& sequence, const std::string& name, bool reverseNode, const std::vector<size_t>& breakpoints);
	void AddEdgeNodeId(int node_id_from, int node_id_to, size_t startOffset);
	void Finalize(int wordSize);
	std::pair<int, size_t> GetReversePosition(int nodeId, size_t offset) const;
	size_t GetReverseNode(size_t node) const;
	size_t NodeSize() const;
	size_t SizeInBP() const;
	size_t BigraphNodeID(size_t directedNodeId) const;
	size_t NodeLength(size_t nodeIndex) const;
	size_t NodeOffset(size_t directedNodeId) const;
	char NodeSequences(size_t node, size_t offset) const;
	NodeChunkSequence NodeChunks(size_t node) const;
	AmbiguousChunkSequence AmbiguousNodeChunks(size_t node) const;
	size_t GetUnitigNode(int nodeId, size_t offset) const;
	std::string BigraphNodeName(int nodeId) const;
	size_t BigraphNodeSize(int nodeId) const;
	size_t BigraphNodeCount() const;
	size_t ComponentSize() const;
	size_t ChainApproxPos(size_t bigraphNodeId) const;
	size_t ChainNumber(size_t bigraphNodeId) const;
	size_t ComponentNumber(size_t digraphNodeId) const;
	NodeEdgeIterator OutNeighbors(size_t nodeId) const;
	NodeEdgeIterator InNeighbors(size_t nodeId) const;
	bool Reverse(size_t digraphNodeId) const;
	bool Linearizable(size_t digraphNodeId) const;
	bool Finalized() const;
	size_t FirstAmbiguous() const;
	static AlignmentGraph DummyGraph();
private:
	void buildEdges();
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
	std::vector<std::string> originalNodeName;
	std::vector<std::vector<size_t>> nodeLookup;
	std::vector<size_t> originalNodeSize;
	std::vector<size_t> componentNumber;
	std::vector<size_t> chainNumber;
	std::vector<size_t> chainApproxPos;
	std::vector<uint8_t> nodeLength;
	std::vector<size_t> nodeOffset;
	std::vector<int> nodeIDs;
	std::vector<bool> hasImplicitOutEdge;
	std::vector<size_t> explicitEdges;
	std::vector<size_t> edgeStorage;
	std::vector<bool> reverse;
	std::vector<bool> linearizable;
	std::vector<NodeChunkSequence> nodeSequences;
	std::vector<AmbiguousChunkSequence> ambiguousNodeSequences;
	std::set<std::tuple<size_t, size_t, size_t>> tempConstructionOutEdges;
	std::vector<bool> ambiguousNodes;
	size_t bpSize;
	size_t firstAmbiguous;
	bool finalized;
};


#endif
