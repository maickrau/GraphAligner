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
	void Finalize();
	size_t SizeInBp() const;
	std::vector<MatrixPosition> GetSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const;

private:
	std::vector<bool> notInOrder;
	std::vector<size_t> nodeStart;
	std::vector<size_t> nodeEnd;
	std::vector<size_t> indexToNode;
	std::map<int, size_t> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::vector<size_t>> inNeighbors;
	std::vector<std::vector<size_t>> outNeighbors;
	std::vector<bool> reverse;
	std::string nodeSequences;
	size_t dummyNodeStart;
	size_t dummyNodeEnd;
	bool finalized;
	size_t firstInOrder;

	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAligner;
};

#endif
