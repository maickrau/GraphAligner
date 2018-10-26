#include <fstream>
#include <sstream>
#include <cassert>
#include <unordered_map>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "stream.hpp"

DirectedGraph::Node::Node(int nodeId, int originalNodeId, bool rightEnd, std::string sequence) :
nodeId(nodeId),
originalNodeId(originalNodeId),
rightEnd(rightEnd),
sequence(sequence)
{
}

DirectedGraph::Edge::Edge(size_t from, size_t to) :
fromId(from),
toId(to)
{
}

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertVGNodeToNodes(const vg::Node& node)
{
	return std::make_pair(DirectedGraph::Node { node.id() * 2, node.id(), true, node.sequence() }, DirectedGraph::Node { node.id() * 2 + 1, node.id(), false, CommonUtils::ReverseComplement(node.sequence()) });
}

std::pair<DirectedGraph::Edge, DirectedGraph::Edge> DirectedGraph::ConvertVGEdgeToEdges(const vg::Edge& edge)
{
	size_t fromLeft, fromRight, toLeft, toRight;
	if (edge.from_start())
	{
		fromLeft = edge.from() * 2;
		fromRight = edge.from() * 2 + 1;
	}
	else
	{
		fromLeft = edge.from() * 2 + 1;
		fromRight = edge.from() * 2;
	}
	if (edge.to_end())
	{
		toLeft = edge.to() * 2;
		toRight = edge.to() * 2 + 1;
	}
	else
	{
		toLeft = edge.to() * 2 + 1;
		toRight = edge.to() * 2;
	}
	return std::make_pair(DirectedGraph::Edge { fromRight, toRight }, DirectedGraph::Edge { toLeft, fromLeft });
}

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertGFANodeToNodes(const std::string& node)
{
	std::stringstream str { node };
	std::string dummy;
	std::string sequence;
	int id;
	str >> dummy >> id >> sequence;
	assert(dummy == "S");
	return ConvertGFANodeToNodes(id, sequence);
}

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertGFANodeToNodes(int id, const std::string& sequence)
{
	return std::make_pair(DirectedGraph::Node { id * 2, id, true, sequence }, DirectedGraph::Node { id * 2 + 1, id, false, CommonUtils::ReverseComplement(sequence) });
}

std::pair<DirectedGraph::Edge, DirectedGraph::Edge> DirectedGraph::ConvertGFAEdgeToEdges(const std::string& edge)
{
	std::stringstream str { edge };
	std::string dummy;
	int from;
	int to;
	std::string fromstart;
	std::string toend;
	str >> dummy >> from >> fromstart >> to >> toend;
	assert(dummy == "L");
	return ConvertGFAEdgeToEdges(from, fromstart, to, toend);
}

std::pair<DirectedGraph::Edge, DirectedGraph::Edge> DirectedGraph::ConvertGFAEdgeToEdges(int from, const std::string& fromstart, int to, const std::string& toend)
{
	assert(fromstart == "+" || fromstart == "-");
	assert(toend == "+" || toend == "-");
	size_t fromLeft, fromRight, toLeft, toRight;
	if (fromstart == "-")
	{
		fromLeft = from * 2;
		fromRight = from * 2 + 1;
	}
	else
	{
		fromLeft = from * 2 + 1;
		fromRight = from * 2;
	}
	if (toend == "-")
	{
		toLeft = to * 2;
		toRight = to * 2 + 1;
	}
	else
	{
		toLeft = to * 2 + 1;
		toRight = to * 2;
	}
	return std::make_pair(DirectedGraph::Edge { fromRight, toRight }, DirectedGraph::Edge { toLeft, fromLeft });
}

AlignmentGraph DirectedGraph::StreamVGGraphFromFile(std::string filename)
{
	AlignmentGraph result;
	{
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Graph&)> lambda = [&result](vg::Graph& g) {
			for (int i = 0; i < g.node_size(); i++)
			{
				auto nodes = ConvertVGNodeToNodes(g.node(i));
				result.AddNode(nodes.first.nodeId, nodes.first.sequence, !nodes.first.rightEnd);
				result.AddNode(nodes.second.nodeId, nodes.second.sequence, !nodes.second.rightEnd);
			}
		};
		stream::for_each(graphfile, lambda);
	}
	{
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Graph&)> lambda = [&result](vg::Graph& g) {
			for (int i = 0; i < g.edge_size(); i++)
			{
				auto edges = ConvertVGEdgeToEdges(g.edge(i));
				result.AddEdgeNodeId(edges.first.fromId, edges.first.toId);
				result.AddEdgeNodeId(edges.second.fromId, edges.second.toId);
			}
		};
		stream::for_each(graphfile, lambda);
	}
	result.Finalize(64);
	return result;
}

AlignmentGraph DirectedGraph::StreamGFAGraphFromFile(std::string filename)
{
	AlignmentGraph result;
	result.DBGOverlap = 0;
	std::vector<Edge> edges;
	std::vector<std::pair<int, std::string>> nodes;
	{
		std::ifstream graphfile { filename, std::ios::in };
		while (graphfile.good())
		{
			std::string line;
			std::getline(graphfile, line);
			if (line[0] == 'L')
			{
				int from, to;
				std::string dummy, fromstart, toend;
				std::string overlapstr;
				std::stringstream str { line };
				str >> dummy >> from >> fromstart >> to >> toend >> overlapstr;
				assert(dummy == "L");
				auto pair = ConvertGFAEdgeToEdges(from, fromstart, to, toend);
				edges.push_back(pair.first);
				edges.push_back(pair.second);
				int overlap = std::stoi(overlapstr.substr(0, overlapstr.size()-1));
				assert(result.DBGOverlap == 0 || result.DBGOverlap == overlap);
				result.DBGOverlap = overlap;
			}
			else if (line[0] == 'S')
			{
				int id;
				std::string dummy;
				std::string sequence;
				std::stringstream str { line };
				str >> dummy >> id >> sequence;
				assert(dummy == "S");
				nodes.emplace_back(id, sequence);
			}
		}
	}
	for (auto node : nodes)
	{
		auto nodes = ConvertGFANodeToNodes(node.first, node.second);
		result.AddNode(nodes.first.nodeId, nodes.first.sequence, !nodes.first.rightEnd);
		result.AddNode(nodes.second.nodeId, nodes.second.sequence, !nodes.second.rightEnd);
	}
	for (auto edge : edges)
	{
		result.AddEdgeNodeId(edge.fromId, edge.toId);
	}
	result.Finalize(64);
	return result;
}

AlignmentGraph DirectedGraph::BuildFromVG(const vg::Graph& graph)
{
	AlignmentGraph result;
	for (int i = 0; i < graph.node_size(); i++)
	{
		auto nodes = ConvertVGNodeToNodes(graph.node(i));
		result.AddNode(nodes.first.nodeId, nodes.first.sequence, !nodes.first.rightEnd);
		result.AddNode(nodes.second.nodeId, nodes.second.sequence, !nodes.second.rightEnd);
	}
	for (int i = 0; i < graph.edge_size(); i++)
	{
		auto edges = ConvertVGEdgeToEdges(graph.edge(i));
		result.AddEdgeNodeId(edges.first.fromId, edges.first.toId);
		result.AddEdgeNodeId(edges.second.fromId, edges.second.toId);
	}
	result.Finalize(64);
	return result;
}

AlignmentGraph DirectedGraph::BuildFromGFA(const GfaGraph& graph)
{
	AlignmentGraph result;
	result.DBGOverlap = graph.edgeOverlap;
	for (auto node : graph.nodes)
	{
		auto nodes = ConvertGFANodeToNodes(node.first, node.second);
		result.AddNode(nodes.first.nodeId, nodes.first.sequence, !nodes.first.rightEnd);
		result.AddNode(nodes.second.nodeId, nodes.second.sequence, !nodes.second.rightEnd);
	}
	for (auto edge : graph.edges)
	{
		for (auto target : edge.second)
		{
			auto pair = ConvertGFAEdgeToEdges(edge.first.id, edge.first.end ? "+" : "-", target.id, target.end ? "+" : "-");
			result.AddEdgeNodeId(pair.first.fromId, pair.first.toId);
			result.AddEdgeNodeId(pair.second.fromId, pair.second.toId);
		}
	}
	result.Finalize(64);
	return result;
}
