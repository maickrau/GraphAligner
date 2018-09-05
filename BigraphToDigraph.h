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
	static std::pair<Node, Node> ConvertGFANodeToNodes(const std::string& line, int edgeOverlap);
	static std::pair<Edge, Edge> ConvertGFAEdgeToEdges(const std::string& line);
	static AlignmentGraph StreamGFAGraphFromFile(std::string filename);
private:
};

#endif
