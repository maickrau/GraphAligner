#ifndef BigraphToDigraph_H
#define BigraphToDigraph_H

#include <tuple>
#include <vector>
#include <string>
#include "AlignmentGraph.h"
#include "vg.pb.h"
#include "GfaGraph.h"

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
		size_t fromId;
		size_t toId;
	};
	static std::pair<Node, Node> ConvertVGNodeToNodes(const vg::Node& node);
	static std::pair<Edge, Edge> ConvertVGEdgeToEdges(const vg::Edge& edge);
	static std::pair<Node, Node> ConvertGFANodeToNodes(int id, const std::string& seq);
	static std::pair<Edge, Edge> ConvertGFAEdgeToEdges(int from, const std::string& fromStart, int to, const std::string& toEnd);
	static AlignmentGraph BuildFromVG(const vg::Graph& graph);
	static AlignmentGraph BuildFromGFA(const GfaGraph& graph);
	static AlignmentGraph StreamVGGraphFromFile(std::string filename);
private:
};

#endif
