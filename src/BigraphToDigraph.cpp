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

static std::vector<bool> getAllowedNucleotides()
{
	std::vector<bool> result;
	result.resize(256, false);
	result['a'] = true;
	result['A'] = true;
	result['c'] = true;
	result['C'] = true;
	result['g'] = true;
	result['G'] = true;
	result['t'] = true;
	result['T'] = true;
	result['y'] = true;
	result['Y'] = true;
	result['r'] = true;
	result['R'] = true;
	result['w'] = true;
	result['W'] = true;
	result['s'] = true;
	result['S'] = true;
	result['k'] = true;
	result['K'] = true;
	result['m'] = true;
	result['M'] = true;
	result['d'] = true;
	result['D'] = true;
	result['v'] = true;
	result['V'] = true;
	result['h'] = true;
	result['H'] = true;
	result['b'] = true;
	result['B'] = true;
	result['n'] = true;
	result['N'] = true;
	return result;
}

auto allowed = getAllowedNucleotides();

DirectedGraph::Node::Node(int nodeId, int originalNodeId, bool rightEnd, std::string sequence, std::string name) :
nodeId(nodeId),
originalNodeId(originalNodeId),
rightEnd(rightEnd),
sequence(sequence),
name(name)
{
}

DirectedGraph::Edge::Edge(size_t from, size_t to, size_t overlap) :
fromId(from),
toId(to),
overlap(overlap)
{
}

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertVGNodeToNodes(const vg::Node& node)
{
	assert(node.id() < std::numeric_limits<int>::max() / 2);
	assert(node.id()+1 < std::numeric_limits<int>::max() / 2);
	return std::make_pair(DirectedGraph::Node { (int)node.id() * 2, (int)node.id(), true, node.sequence(), node.name() }, DirectedGraph::Node { (int)node.id() * 2 + 1, (int)node.id(), false, CommonUtils::ReverseComplement(node.sequence()), node.name() });
}

std::pair<DirectedGraph::Edge, DirectedGraph::Edge> DirectedGraph::ConvertVGEdgeToEdges(const vg::Edge& edge)
{
	assert(edge.overlap() == 0);
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
	return std::make_pair(DirectedGraph::Edge { fromRight, toRight, 0 }, DirectedGraph::Edge { toLeft, fromLeft, 0 });
}

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertGFANodeToNodes(int id, const std::string& sequence, const std::string& name)
{
	return std::make_pair(DirectedGraph::Node { id * 2, id, true, sequence, name }, DirectedGraph::Node { id * 2 + 1, id, false, CommonUtils::ReverseComplement(sequence), name });
}

std::pair<DirectedGraph::Edge, DirectedGraph::Edge> DirectedGraph::ConvertGFAEdgeToEdges(int from, const std::string& fromstart, int to, const std::string& toend, size_t overlap)
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
	return std::make_pair(DirectedGraph::Edge { fromRight, toRight, overlap }, DirectedGraph::Edge { toLeft, fromLeft, overlap });
}

AlignmentGraph DirectedGraph::StreamVGGraphFromFile(std::string filename)
{
	AlignmentGraph result;
	{
		std::vector<size_t> breakpointsFw;
		std::vector<size_t> breakpointsBw;
		breakpointsFw.push_back(0);
		breakpointsBw.push_back(0);
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Graph&)> lambda = [&result, &breakpointsFw, &breakpointsBw](vg::Graph& g) {
			for (int i = 0; i < g.node_size(); i++)
			{
				for (size_t j = 0; j < g.node(i).sequence().size(); j++)
				{
					if (!allowed[g.node(i).sequence()[j]])
					{
						throw CommonUtils::InvalidGraphException("Invalid sequence character: " + g.node(i).sequence()[j]);
					}
				}
				auto nodes = ConvertVGNodeToNodes(g.node(i));
				assert(nodes.first.sequence.size() == nodes.second.sequence.size());
				breakpointsFw.push_back(g.node(i).sequence().size());
				breakpointsBw.push_back(g.node(i).sequence().size());
				result.AddNode(nodes.first.nodeId, nodes.first.sequence, nodes.first.name, !nodes.first.rightEnd, breakpointsFw);
				result.AddNode(nodes.second.nodeId, nodes.second.sequence, nodes.second.name, !nodes.second.rightEnd, breakpointsBw);
				breakpointsFw.erase(breakpointsFw.begin()+1, breakpointsFw.end());
				breakpointsBw.erase(breakpointsBw.begin()+1, breakpointsBw.end());
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
				result.AddEdgeNodeId(edges.first.fromId, edges.first.toId, edges.first.overlap);
				result.AddEdgeNodeId(edges.second.fromId, edges.second.toId, edges.second.overlap);
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
	std::vector<size_t> breakpointsFw;
	std::vector<size_t> breakpointsBw;
	breakpointsFw.push_back(0);
	breakpointsBw.push_back(0);
	for (int i = 0; i < graph.node_size(); i++)
	{
		for (size_t j = 0; j < graph.node(i).sequence().size(); j++)
		{
			if (!allowed[graph.node(i).sequence()[j]])
			{
				throw CommonUtils::InvalidGraphException("Invalid sequence character: " + graph.node(i).sequence()[j]);
			}
		}
		auto nodes = ConvertVGNodeToNodes(graph.node(i));
		breakpointsFw.push_back(graph.node(i).sequence().size());
		breakpointsBw.push_back(graph.node(i).sequence().size());
		result.AddNode(nodes.first.nodeId, nodes.first.sequence, nodes.first.name, !nodes.first.rightEnd, breakpointsFw);
		result.AddNode(nodes.second.nodeId, nodes.second.sequence, nodes.second.name, !nodes.second.rightEnd, breakpointsBw);
		breakpointsFw.erase(breakpointsFw.begin()+1, breakpointsFw.end());
		breakpointsBw.erase(breakpointsBw.begin()+1, breakpointsBw.end());
	}
	for (int i = 0; i < graph.edge_size(); i++)
	{
		auto edges = ConvertVGEdgeToEdges(graph.edge(i));
		result.AddEdgeNodeId(edges.first.fromId, edges.first.toId, edges.first.overlap);
		result.AddEdgeNodeId(edges.second.fromId, edges.second.toId, edges.second.overlap);
	}
	result.Finalize(64);
	return result;
}

AlignmentGraph DirectedGraph::BuildFromGFA(const GfaGraph& graph)
{
	AlignmentGraph result;
	result.DBGoverlap = graph.edgeOverlap;
	std::unordered_map<int, std::vector<size_t>> breakpoints;
	for (auto pair : graph.varyingOverlaps)
	{
		int to = pair.first.second.id * 2;
		if (!pair.first.second.end) to += 1;
		int from = pair.first.first.Reverse().id * 2;
		if (!pair.first.first.Reverse().end) from += 1;
		breakpoints[from].push_back(pair.second);
		breakpoints[to].push_back(pair.second);
	}
	for (auto node : graph.nodes)
	{
		for (size_t j = 0; j < node.second.size(); j++)
		{
			if (!allowed[node.second[j]])
			{
				throw CommonUtils::InvalidGraphException("Invalid sequence character: " + node.second[j]);
			}
		}
		std::string name = graph.OriginalNodeName(node.first);
		auto nodes = ConvertGFANodeToNodes(node.first, node.second, name);
		std::vector<size_t> breakpointsFw = breakpoints[node.first * 2];
		std::vector<size_t> breakpointsBw = breakpoints[node.first * 2 + 1];
		breakpointsFw.push_back(0);
		breakpointsFw.push_back(node.second.size());
		breakpointsBw.push_back(0);
		breakpointsBw.push_back(node.second.size());
		std::sort(breakpointsFw.begin(), breakpointsFw.end());
		std::sort(breakpointsBw.begin(), breakpointsBw.end());
		result.AddNode(nodes.first.nodeId, nodes.first.sequence, nodes.first.name, !nodes.first.rightEnd, breakpointsFw);
		result.AddNode(nodes.second.nodeId, nodes.second.sequence, nodes.second.name, !nodes.second.rightEnd, breakpointsBw);
	}
	for (auto edge : graph.edges)
	{
		for (auto target : edge.second)
		{
			auto overlap = graph.edgeOverlap;
			if (graph.varyingOverlaps.count(std::make_pair(edge.first, target)) == 1)
			{
				overlap = graph.varyingOverlaps.at(std::make_pair(edge.first, target));
			}
			auto pair = ConvertGFAEdgeToEdges(edge.first.id, edge.first.end ? "+" : "-", target.id, target.end ? "+" : "-", overlap);
			result.AddEdgeNodeId(pair.first.fromId, pair.first.toId, pair.first.overlap);
			result.AddEdgeNodeId(pair.second.fromId, pair.second.toId, pair.second.overlap);
		}
	}
	result.Finalize(64);
	return result;
}
