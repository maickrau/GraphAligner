#ifndef BigraphToDigraph_H
#define BigraphToDigraph_H

#include <vector>
#include <string>
#include "vg.pb.h"

class DirectedGraph
{
public:
	struct SeedHit
	{
	public:
		SeedHit(int nodeId, size_t nodePos, size_t seqPos);
		int nodeId;
		size_t nodePos;
		size_t seqPos;
	};
	struct Node
	{
		Node(int nodeId, int originalNodeId, bool rightEnd, std::string sequence);
		int nodeId;
		int originalNodeId;
		bool rightEnd;
		std::string sequence;
	};
	struct Edge
	{
		Edge(size_t from, size_t to);
		size_t fromIndex;
		size_t toIndex;
	};
	DirectedGraph();
	DirectedGraph(const vg::Graph& bigraph);
	std::vector<Node> nodes;
	std::vector<Edge> edges;
	void ReorderByNodeIds(const std::vector<int>& nodeIdOrder);
	void RemoveNodes(const std::set<int>& nodeIndices);
	void AddSubgraph(const DirectedGraph& subgraph);
	void ConnectComponents(const std::vector<int>& previousSinks, const std::vector<int>& nextSources);
	std::vector<SeedHit> GetSeedHits(const std::string& sequence, const std::vector<std::pair<int, size_t>>& hitsOriginalNodeIds) const;
	void PruneByReachability(const std::vector<int>& startNodeIds);
	size_t totalSequenceLength;
private:
	void addReachable(std::vector<bool>& reachable, const std::vector<std::vector<size_t>>& outNeighbors, size_t current);
	std::tuple<size_t, size_t, size_t> getLocalAlignmentSsw(std::string sequence, size_t nodeIndex) const;
	bool edgesPointToValidNodes();
	bool nodeIdsAreValid();
};

#endif
