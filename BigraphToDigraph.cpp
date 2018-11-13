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
	assert(node.id() < std::numeric_limits<int>::max() / 2);
	assert(node.id()+1 < std::numeric_limits<int>::max() / 2);
	return std::make_pair(DirectedGraph::Node { (int)node.id() * 2, (int)node.id(), true, node.sequence() }, DirectedGraph::Node { (int)node.id() * 2 + 1, (int)node.id(), false, CommonUtils::ReverseComplement(node.sequence()) });
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

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertGFANodeToNodes(int id, const std::string& sequence)
{
	return std::make_pair(DirectedGraph::Node { id * 2, id, true, sequence }, DirectedGraph::Node { id * 2 + 1, id, false, CommonUtils::ReverseComplement(sequence) });
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
