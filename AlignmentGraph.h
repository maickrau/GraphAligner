#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <vector>
#include <map>
#include <tuple>
#include "ThreadReadAssertion.h"

class AlignmentGraph
{
public:
	typedef std::pair<size_t, size_t> MatrixPosition;
	class SeedHit
	{
	public:
		SeedHit(size_t seqPos, int nodeId, size_t nodePos) : sequencePosition(seqPos), nodeId(nodeId), nodePos(nodePos) {};
		size_t sequencePosition;
		int nodeId;
		size_t nodePos;
	};
	AlignmentGraph();
	void AddNode(int nodeId, std::string sequence, bool reverseNode);
	void AddEdgeNodeId(int node_id_from, int node_id_to);
	void Finalize(int wordSize);
	size_t SizeInBp() const;
	std::vector<MatrixPosition> GetSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const;

private:
	void getDAGFromIdenticalSubtrees(const std::vector<size_t>& nodes, const std::vector<size_t>& parents, std::vector<size_t>& resultNodes, std::vector<std::vector<size_t>>& resultPredecessors);
	void getCycleCutterTreeRec(size_t node, size_t parent, int wordSize, int lengthLeft, std::vector<size_t>& nodes, std::vector<size_t>& parents);
	void calculateCycleCutters(size_t cycleStart, int wordSize);
	std::vector<bool> notInOrder;
	std::vector<size_t> nodeStart;
	std::vector<size_t> nodeEnd;
	std::vector<size_t> indexToNode;
	std::map<int, size_t> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::vector<size_t>> inNeighbors;
	std::vector<std::vector<size_t>> outNeighbors;
	std::vector<bool> reverse;
	std::vector<std::vector<size_t>> cycleCuttingNodes;
	std::vector<std::vector<std::vector<size_t>>> cycleCuttingNodePredecessor;
	std::string nodeSequences;
	size_t dummyNodeStart;
	size_t dummyNodeEnd;
	bool finalized;
	size_t firstInOrder;

	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAligner;
};

#endif
