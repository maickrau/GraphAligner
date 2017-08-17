#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <functional>
#include <vector>
#include <unordered_map>
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
	void ReserveNodes(size_t numNodes, size_t totalSequenceLength);
	void AddNode(int nodeId, const std::string& sequence, bool reverseNode);
	void AddEdgeNodeId(int node_id_from, int node_id_to);
	void Finalize(int wordSize, std::string cutFilename);
	size_t SizeInBp() const;
	std::vector<MatrixPosition> GetSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const;

private:
	bool loadCycleCut(std::string filename);
	void saveCycleCut(std::string filename);
	void iterateOverCycleCuttingTreeRec(size_t cycleStart, size_t node, int sizeLeft, std::vector<size_t>& currentStack, std::function<void(const std::vector<size_t>&)> function);
	void iterateOverCycleCuttingTree(size_t cycleStart, int sizeLeft, std::function<void(const std::vector<size_t>&)> function);
	void getCycleCuttersSupersequence(size_t cycleStart, int sizeLeft, std::vector<size_t>& supersequence, std::vector<std::set<size_t>>& supersequencePredecessors, std::vector<bool>& previousCut);
	void calculateCycleCutters(size_t cycleStart, int wordSize);
	std::vector<bool> notInOrder;
	std::vector<size_t> nodeStart;
	std::vector<size_t> nodeEnd;
	std::vector<size_t> indexToNode;
	std::unordered_map<int, size_t> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::set<size_t>> inNeighbors;
	std::vector<std::set<size_t>> outNeighbors;
	std::vector<bool> reverse;
	std::vector<std::vector<size_t>> cycleCuttingNodes;
	std::vector<std::vector<std::set<size_t>>> cycleCuttingNodePredecessor;
	std::vector<std::vector<bool>> cycleCutPreviousCut;
	std::string nodeSequences;
	size_t dummyNodeStart;
	size_t dummyNodeEnd;
	bool finalized;
	size_t firstInOrder;

	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAligner;
};

#endif
