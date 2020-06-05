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
		Node(int nodeId, int originalNodeId, bool rightEnd, std::string sequence, std::string name);
		int nodeId;
		int originalNodeId;
		bool rightEnd;
		std::string sequence;
		std::string name;
	};
	struct Edge
	{
		Edge(size_t from, size_t to, size_t overlap);
		size_t fromId;
		size_t toId;
		size_t overlap;
	};
	static std::pair<Node, Node> ConvertVGNodeToNodes(const vg::Node& node);
	static std::pair<Edge, Edge> ConvertVGEdgeToEdges(const vg::Edge& edge);
	static std::pair<Node, Node> ConvertGFANodeToNodes(int id, const std::string& seq, const std::string& name);
	static std::pair<Edge, Edge> ConvertGFAEdgeToEdges(int from, const std::string& fromStart, int to, const std::string& toEnd, size_t overlap);
	static AlignmentGraph BuildFromVG(const vg::Graph& graph);
	static AlignmentGraph BuildFromGFA(const GfaGraph& graph);
	static AlignmentGraph StreamVGGraphFromFile(std::string filename);
private:
};

#endif
