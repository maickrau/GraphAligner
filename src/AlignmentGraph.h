#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <iterator>
#include <functional>
#include <vector>
#include <tuple>
#include <unordered_set>
#include <set>
#include <phmap.h>
#include "RankBitvector.h"
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
	void ReserveNodes(size_t numBigraphNodes, size_t numDigraphNodes);
	void AddNode(size_t bigraphNodeId, const DNAString& sequence, const std::string& name, bool reverseNode, const std::vector<size_t>& breakpoints);
	void AddEdgeNodeId(size_t bigraphIdFrom, size_t bigraphIdTo, size_t startOffset);
	void Finalize(int wordSize);
	std::pair<size_t, size_t> GetReversePosition(size_t bigraphNodeId, size_t offset) const;
	std::pair<size_t, size_t> GetReverseDigraphPosition(size_t digraphNodeId, size_t offset) const;
	size_t NodeSize() const;
	size_t SizeInBP() const;
	size_t BigraphNodeID(size_t digraphNodeId) const;
	size_t NodeLength(size_t digraphNodeId) const;
	size_t NodeOffset(size_t digraphNodeId) const;
	char NodeSequences(size_t digraphNodeId, size_t offset) const;
	NodeChunkSequence NodeChunks(size_t digraphNodeId) const;
	AmbiguousChunkSequence AmbiguousNodeChunks(size_t digraphNodeId) const;
	size_t GetDigraphNode(size_t bigraphNodeId, size_t offset) const;
	std::string BigraphNodeName(size_t bigraphNodeId) const;
	size_t BigraphNodeSize(size_t bigraphNodeId) const;
	size_t BigraphNodeCount() const;
	size_t ComponentSize() const;
	size_t ChainApproxPos(size_t bigraphNodeId) const;
	size_t ChainNumber(size_t bigraphNodeId) const;
	size_t ComponentNumber(size_t digraphNodeId) const;
	NodeEdgeIterator OutNeighbors(size_t digraphNodeId) const;
	NodeEdgeIterator InNeighbors(size_t digraphNodeId) const;
	bool Reverse(size_t digraphNodeId) const;
	bool Linearizable(size_t digraphNodeId) const;
	bool Finalized() const;
	size_t FirstAmbiguous() const;
	std::string BigraphNodeSeq(size_t bigraphNodeId) const;
	static AlignmentGraph DummyGraph();
	bool AllNodeNamesAreNumbers() const;
private:
	void makeDinodeIntermediateMapping();
	void sparsenComponentNumbers();
	void replaceIntermediateEdgesWithDinodes();
	size_t addIntermediateNodes(size_t bigraphNodeId, const DNAString& sequence, size_t start, size_t end);
	void fixChainApproxPos(const size_t start);
	std::pair<bool, size_t> findBubble(const size_t start, const std::vector<bool>& ignorableTip);
	void chainBubble(const size_t start, const std::vector<bool>& ignorableTip, std::vector<size_t>& rank);
	phmap::flat_hash_map<size_t, std::unordered_set<size_t>> chainTips(std::vector<size_t>& rank, std::vector<bool>& ignorableTip);
	void chainCycles(std::vector<size_t>& rank, std::vector<bool>& ignorableTip);
	size_t intermediateNodeLength(size_t intermediateId) const;
	void findChains();
	void AddAmbiguousDinode(const std::string& sequence);
	void AddNormalDinode(const std::string& sequence);
	void RenumberAmbiguousToEnd();
	void doComponentOrder();
	size_t intermediateNodeCount() const;
	size_t digraphToIntermediate(size_t digraphNodeId) const;
	size_t intermediateLastDinode(size_t intermediate) const;
	size_t intermediateDinodesCount(size_t intermediate) const;
	size_t bpSize;
	size_t firstAmbiguous;
	bool finalized;
	bool allNodeNamesAreNumbers;
	// bigraph
	std::vector<std::string> originalNodeName;
	std::vector<std::vector<size_t>> bigraphIntermediateList;
	std::vector<size_t> originalNodeSize;
	std::vector<size_t> chainNumber;
	std::vector<size_t> chainApproxPos;
	std::vector<bool> reverse;
	// intermediates
	std::vector<bool> partOfStronglyConnectedComponent;
	std::vector<size_t> componentNumber;
	std::vector<uint8_t> lastDinodeLength;
	std::vector<size_t> firstDinodeOffset;
	std::vector<size_t> intermediateBigraphNodeIDs;
	std::vector<size_t> intermediateDinodesStart;
	std::vector<std::vector<size_t>> intermediateInEdges; // during construction points to intermediates, after points to dinodes
	std::vector<std::vector<size_t>> intermediateOutEdges; // during construction points to intermediates, after points to dinodes
	// digraph
	RankBitvector firstOfIntermediates;
	std::vector<NodeChunkSequence> nodeSequences;
	std::vector<AmbiguousChunkSequence> ambiguousNodeSequences;
};


#endif
