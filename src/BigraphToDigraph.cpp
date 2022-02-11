#include <fstream>
#include <sstream>
#include <cassert>
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

DirectedGraph::Node::Node(size_t nodeId, bool rightEnd, std::string sequence, std::string name) :
nodeId(nodeId),
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

size_t getName(std::unordered_map<int, size_t>& nameMapping, int name)
{
	size_t result = nameMapping.size();
	if (nameMapping.count(name) == 1)
	{
		return nameMapping.at(name);
	}
	nameMapping[name] = result;
	return result;
}

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertVGNodeToNodes(const vg::Node& node, std::unordered_map<int, size_t>& nameMapping)
{
	size_t id = getName(nameMapping, (int)node.id());
	std::string name = std::to_string(node.id());
	return std::make_pair(DirectedGraph::Node { id * 2, true, node.sequence(), name }, DirectedGraph::Node { id * 2 + 1, false, CommonUtils::ReverseComplement(node.sequence()), name });
}

std::pair<DirectedGraph::Edge, DirectedGraph::Edge> DirectedGraph::ConvertVGEdgeToEdges(const vg::Edge& edge, std::unordered_map<int, size_t>& nameMapping)
{
	size_t from = getName(nameMapping, (int)edge.from());
	size_t to = getName(nameMapping, (int)edge.to());
	assert(edge.overlap() == 0);
	size_t fromLeft, fromRight, toLeft, toRight;
	if (edge.from_start())
	{
		fromLeft = from * 2;
		fromRight = from * 2 + 1;
	}
	else
	{
		fromLeft = from * 2 + 1;
		fromRight = from * 2;
	}
	if (edge.to_end())
	{
		toLeft = to * 2;
		toRight = to * 2 + 1;
	}
	else
	{
		toLeft = to * 2 + 1;
		toRight = to * 2;
	}
	return std::make_pair(DirectedGraph::Edge { fromRight, toRight, 0 }, DirectedGraph::Edge { toLeft, fromLeft, 0 });
}

std::pair<DirectedGraph::Node, DirectedGraph::Node> DirectedGraph::ConvertGFANodeToNodes(int id, const std::string& sequence, const std::string& name)
{
	return std::make_pair(DirectedGraph::Node { id * 2, true, sequence, name }, DirectedGraph::Node { id * 2 + 1, false, CommonUtils::ReverseComplement(sequence), name });
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
	std::unordered_map<int, size_t> nameMapping;
	AlignmentGraph result;
	size_t nodeCount = 0;
	{
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Graph&)> countlambda = [&nodeCount](vg::Graph& g) {
			nodeCount += g.node_size();
		};
		stream::for_each(graphfile, countlambda);
	}
	result.ReserveNodes(nodeCount*2, nodeCount*2);
	{
		std::vector<size_t> breakpointsFw;
		std::vector<size_t> breakpointsBw;
		breakpointsFw.push_back(0);
		breakpointsBw.push_back(0);
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Graph&)> lambda = [&result, &breakpointsFw, &breakpointsBw, &nameMapping](vg::Graph& g) {
			for (int i = 0; i < g.node_size(); i++)
			{
				for (size_t j = 0; j < g.node(i).sequence().size(); j++)
				{
					if (!allowed[g.node(i).sequence()[j]])
					{
						throw CommonUtils::InvalidGraphException("Invalid sequence character: " + g.node(i).sequence()[j]);
					}
				}
				auto nodes = ConvertVGNodeToNodes(g.node(i), nameMapping);
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
		std::function<void(vg::Graph&)> lambda = [&result, &nameMapping](vg::Graph& g) {
			for (int i = 0; i < g.edge_size(); i++)
			{
				auto edges = ConvertVGEdgeToEdges(g.edge(i), nameMapping);
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
	std::unordered_map<int, size_t> nameMapping;
	AlignmentGraph result;
	result.ReserveNodes(graph.node_size()*2, graph.node_size()*2);
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
		auto nodes = ConvertVGNodeToNodes(graph.node(i), nameMapping);
		breakpointsFw.push_back(graph.node(i).sequence().size());
		breakpointsBw.push_back(graph.node(i).sequence().size());
		result.AddNode(nodes.first.nodeId, nodes.first.sequence, nodes.first.name, !nodes.first.rightEnd, breakpointsFw);
		result.AddNode(nodes.second.nodeId, nodes.second.sequence, nodes.second.name, !nodes.second.rightEnd, breakpointsBw);
		breakpointsFw.erase(breakpointsFw.begin()+1, breakpointsFw.end());
		breakpointsBw.erase(breakpointsBw.begin()+1, breakpointsBw.end());
	}
	for (int i = 0; i < graph.edge_size(); i++)
	{
		auto edges = ConvertVGEdgeToEdges(graph.edge(i), nameMapping);
		result.AddEdgeNodeId(edges.first.fromId, edges.first.toId, edges.first.overlap);
		result.AddEdgeNodeId(edges.second.fromId, edges.second.toId, edges.second.overlap);
	}
	result.Finalize(64);
	return result;
}

AlignmentGraph DirectedGraph::BuildFromGFA(const GfaGraph& graph)
{
	AlignmentGraph result;
	result.ReserveNodes(graph.nodes.size()*2, graph.nodes.size()*2);
	std::vector<std::vector<size_t>> breakpoints;
	breakpoints.resize(graph.nodes.size() * 2);
	for (const auto& t : graph.edges)
	{
		int to = std::get<1>(t).id * 2;
		if (!std::get<1>(t).end) to += 1;
		int from = std::get<0>(t).Reverse().id * 2;
		if (!std::get<0>(t).Reverse().end) from += 1;
		breakpoints[from].push_back(std::get<2>(t));
		breakpoints[to].push_back(std::get<2>(t));
	}
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		for (size_t j = 0; j < graph.nodes[i].size(); j++)
		{
			if (!allowed[graph.nodes[i][j]])
			{
				throw CommonUtils::InvalidGraphException("Invalid sequence character: " + graph.nodes[i][j]);
			}
		}
		std::string name = graph.OriginalNodeName(i);
		auto nodes = ConvertGFANodeToNodes(i, graph.nodes[i], name);
		std::vector<size_t> breakpointsFw = breakpoints[i * 2];
		std::vector<size_t> breakpointsBw = breakpoints[i * 2 + 1];
		breakpointsFw.push_back(0);
		breakpointsFw.push_back(graph.nodes[i].size());
		breakpointsBw.push_back(0);
		breakpointsBw.push_back(graph.nodes[i].size());
		std::sort(breakpointsFw.begin(), breakpointsFw.end());
		std::sort(breakpointsBw.begin(), breakpointsBw.end());
		result.AddNode(nodes.first.nodeId, nodes.first.sequence, nodes.first.name, !nodes.first.rightEnd, breakpointsFw);
		result.AddNode(nodes.second.nodeId, nodes.second.sequence, nodes.second.name, !nodes.second.rightEnd, breakpointsBw);
	}
	for (const auto& t : graph.edges)
	{
		auto pair = ConvertGFAEdgeToEdges(std::get<0>(t).id, std::get<0>(t).end ? "+" : "-", std::get<1>(t).id, std::get<1>(t).end ? "+" : "-", std::get<2>(t));
		result.AddEdgeNodeId(pair.first.fromId, pair.first.toId, pair.first.overlap);
		result.AddEdgeNodeId(pair.second.fromId, pair.second.toId, pair.second.overlap);
	}
	result.Finalize(64);
	return result;
}
