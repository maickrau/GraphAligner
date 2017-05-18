#ifndef BigraphToDigraph_H
#define BigraphToDigraph_H

#include <vector>
#include <string>
#include "vg.pb.h"

class DirectedGraph
{
public:
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
	DirectedGraph(const vg::Graph& bigraph);
	std::vector<Node> nodes;
	std::vector<Edge> edges;
	void ReorderByNodeIds(const std::vector<int>& nodeIdOrder);
	void RemoveNodes(const std::set<int>& nodeIndices);
private:
	bool edgesPointToValidNodes();
	bool nodeIdsAreValid();
};

#endif
