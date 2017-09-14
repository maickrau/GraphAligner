#ifndef AlignmentGraph_h
#define AlignmentGraph_h

#include <functional>
#include <vector>
#include <set>
#include <unordered_map>
#include <tuple>
#include "ThreadReadAssertion.h"

class CycleCutCalculation;

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
	class CycleCut
	{
	public:
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & nodes;
			ar & predecessors;
			ar & previousCut;
		}
		std::vector<size_t> nodes;
		std::vector<std::set<size_t>> predecessors;
		std::vector<bool> previousCut;
	};
	AlignmentGraph();
	void ReserveNodes(size_t numNodes, size_t totalSequenceLength);
	void AddNode(int nodeId, const std::string& sequence, bool reverseNode);
	void AddEdgeNodeId(int node_id_from, int node_id_to);
	void Finalize(int wordSize, std::string cutFilename);
	size_t GetReversePosition(size_t position) const;
	size_t SizeInBp() const;
	std::set<size_t> ProjectForward(const std::set<size_t>& startpositions, size_t amount) const;
	std::vector<MatrixPosition> GetSeedHitPositionsInMatrix(const std::string& sequence, const std::vector<SeedHit>& seedHits) const;

private:
	void calculateCycleCuts(int wordSize);
	bool loadCycleCut(std::string filename);
	void saveCycleCut(std::string filename);
	void calculateCycleCutters(const CycleCutCalculation& cutCalculator, size_t cycleStart, int wordSize);
	std::vector<bool> notInOrder;
	std::vector<size_t> nodeStart;
	std::vector<size_t> nodeEnd;
	std::vector<size_t> indexToNode;
	std::unordered_map<int, size_t> nodeLookup;
	std::vector<int> nodeIDs;
	std::vector<std::set<size_t>> inNeighbors;
	std::vector<std::set<size_t>> outNeighbors;
	std::vector<bool> reverse;
	std::vector<CycleCut> cuts;
	std::string nodeSequences;
	size_t dummyNodeStart;
	size_t dummyNodeEnd;
	bool finalized;
	size_t firstInOrder;

	friend class CycleCutCalculation;
	template <typename LengthType, typename ScoreType, typename Word>
	friend class GraphAligner;
};


#endif
