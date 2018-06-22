#ifndef BigraphToDigraph_H
#define BigraphToDigraph_H

#include <tuple>
#include <vector>
#include <string>
#include "AlignmentGraph.h"

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
	static std::pair<Node, Node> ConvertGFANodeToNodes(int id, const std::string& seq, int edgeOverlap);
	static std::pair<Node, Node> ConvertGFANodeToNodes(const std::string& line, int edgeOverlap);
	static std::pair<Edge, Edge> ConvertGFAEdgeToEdges(int from, const std::string& fromStart, int to, const std::string& toEnd);
	static std::pair<Edge, Edge> ConvertGFAEdgeToEdges(const std::string& line);
	static AlignmentGraph StreamVGGraphFromFile(std::string filename);
	static AlignmentGraph StreamGFAGraphFromFile(std::string filename);
private:
};

#endif
