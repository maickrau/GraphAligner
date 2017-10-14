#ifndef GfaGraph_h
#define GfaGraph_h

#include <vector>
#include <string>

class GfaGraph
{
public:
	class Edge
	{
	public:
		int from;
		int to;
		bool fromStart;
		bool toEnd;
	};
	class Node
	{
	public:
		int id;
		std::string sequence;
	};
	static GfaGraph LoadFromFile(std::string filename);
	std::vector<Node> nodes;
	std::vector<Edge> edges;
	int edgeOverlap;
private:
	GfaGraph();
};

#endif