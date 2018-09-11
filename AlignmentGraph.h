#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <functional>
#include <vector>
#include <set>
#include <unordered_map>
#include <tuple>
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

	struct MatrixPosition
	{
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
	void AddNode(int nodeId, const std::string& sequence, bool reverseNode);
	void AddEdgeNodeId(int node_id_from, int node_id_to);
	void Finalize(int wordSize);
	AlignmentGraph GetSubgraph(const std::unordered_map<size_t, size_t>& nodeMapping) const;
	std::pair<int, size_t> GetReversePosition(int nodeId, size_t offset) const;
	size_t GetReverseNode(size_t node) const;
	size_t NodeSize() const;
	size_t NodeLength(size_t nodeIndex) const;
	char NodeSequences(size_t node, size_t offset) const;
	NodeChunkSequence NodeChunks(size_t node) const;
	size_t GetUnitigNode(int nodeId, size_t offset) const;
	// size_t MinDistance(size_t pos, const std::vector<size_t>& targets) const;
	// std::set<size_t> ProjectForward(const std::set<size_t>& startpositions, size_t amount) const;
	int DBGOverlap;

private:
	void AddNode(int nodeId, int offset, const std::string& sequence, bool reverseNode);
	std::vector<size_t> nodeLength;
	std::unordered_map<int, std::vector<size_t>> nodeLookup;
	std::unordered_map<int, size_t> unitigStartNode;
	std::unordered_map<int, size_t> originalNodeSize;
	std::vector<size_t> nodeOffset;
	std::vector<int> nodeIDs;
	std::vector<std::vector<size_t>> inNeighbors;
	std::vector<std::vector<size_t>> outNeighbors;
	std::vector<bool> reverse;
	std::vector<NodeChunkSequence> nodeSequences;
	bool finalized;

	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAligner;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerVGAlignment;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerBitvectorBanded;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class SubgraphExtractor;
};


#endif
