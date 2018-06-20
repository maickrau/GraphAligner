#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <functional>
#include <vector>
#include <set>
#include <unordered_map>
#include <tuple>
#include "ThreadReadAssertion.h"
#include "DNAString.h"

class AlignmentGraph
{
public:
	//arbitrarily 64
	//determines extra band size, shouldn't be too high because of extra slices
	//should be 0 mod wordsize, otherwise storage has overhead
	static constexpr int SPLIT_NODE_SIZE = 64;

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
	void ReserveNodes(size_t numNodes, size_t numSplitNodes, size_t totalSequenceLength);
	void AddNode(int nodeId, const std::string& sequence, bool reverseNode);
	void AddEdgeNodeId(int node_id_from, int node_id_to);
	void Finalize(int wordSize);
	size_t GetReversePosition(size_t position) const;
	size_t SizeInBp() const;
	size_t IndexToNode(size_t index) const;
	size_t NodeSize() const;
	size_t NodeStart(size_t nodeIndex) const;
	size_t NodeEnd(size_t nodeIndex) const;
	size_t NodeLength(size_t nodeIndex) const;
	char NodeSequences(size_t index) const;
	size_t NodeSequencesI(size_t index) const;
	size_t NodeSequencesSize() const;
	size_t MinDistance(size_t pos, const std::vector<size_t>& targets) const;
	std::set<size_t> ProjectForward(const std::set<size_t>& startpositions, size_t amount) const;
	std::vector<MatrixPosition> GetSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const;
	int DBGOverlap;

private:
	void AddNode(int nodeId, int offset, std::string sequence, bool reverseNode);
	std::vector<size_t> nodeStart;
	std::unordered_map<int, std::vector<size_t>> nodeLookup;
	std::vector<size_t> nodeOffset;
	std::vector<int> nodeIDs;
	std::vector<std::vector<size_t>> inNeighbors;
	std::vector<std::vector<size_t>> outNeighbors;
	std::vector<bool> reverse;
	DNAString<size_t> nodeSequences;
	bool finalized;

	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAligner;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerVGAlignment;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAlignerBitvectorBanded;
};


#endif
