#include <algorithm>
#include <vector>
#include "TopologicalSort.h"
#include "vg.pb.h"
#include "stream.hpp"

std::vector<int> getNodeIdsInTopologicalOrder(const vg::Graph& graph)
{
	auto indexes = topologicalSort(graph);
	std::vector<int> result;
	result.reserve(indexes.size());
	for (size_t i = 0; i < indexes.size(); i++)
	{
		result.push_back(graph.node(indexes[i]).id());
	}
	return result;
}

template <typename Orderable, typename Function>
void order(std::vector<Orderable>& orderables, Function nodeIdGetter, const vg::Graph& graph)
{
	std::map<int, size_t> nodeIdPositions;
	auto nodeIdsInOrder = getNodeIdsInTopologicalOrder(graph);
	for (size_t i = 0; i < nodeIdsInOrder.size(); i++)
	{
		nodeIdPositions[nodeIdsInOrder[i]] = i;
	}
	std::sort(orderables.begin(), orderables.end(), [&nodeIdPositions, nodeIdGetter](const Orderable& left, const Orderable& right) { return nodeIdPositions[nodeIdGetter(left)] < nodeIdPositions[nodeIdGetter(right)];});
}

vg::Graph mergeGraphs(const std::vector<vg::Graph>& parts)
{
	vg::Graph newGraph;
	std::vector<const vg::Node*> allNodes;
	std::vector<const vg::Edge*> allEdges;
	for (size_t i = 0; i < parts.size(); i++)
	{
		for (int j = 0; j < parts[i].node_size(); j++)
		{
			allNodes.push_back(&parts[i].node(j));
		}
		for (int j = 0; j < parts[i].edge_size(); j++)
		{
			allEdges.push_back(&parts[i].edge(j));
		}
	}
	for (size_t i = 0; i < allNodes.size(); i++)
	{
		auto node = newGraph.add_node();
		node->set_id(allNodes[i]->id());
		node->set_sequence(allNodes[i]->sequence());
		node->set_name(allNodes[i]->name());
	}
	for (size_t i = 0; i < allEdges.size(); i++)
	{
		auto edge = newGraph.add_edge();
		edge->set_from(allEdges[i]->from());
		edge->set_to(allEdges[i]->to());
		edge->set_from_start(allEdges[i]->from_start());
		edge->set_to_end(allEdges[i]->to_end());
		edge->set_overlap(allEdges[i]->overlap());
	}
	return newGraph;
}

bool endsWith(const std::string& str, std::string ending)
{
	if (str.size() < ending.size()) return false;
	return strncmp(str.data() + str.size() - ending.size(), ending.data(), ending.size()) == 0;
}

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		std::cout << "needs 3 parameters";
		std::exit(0);
	}
	vg::Graph graph;
	{
		if (is_file_exist(argv[1])){
			std::cout << "load graph from " << argv[1] << std::endl;
			std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
			std::vector<vg::Graph> parts;
			std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
				parts.push_back(g);
			};
			stream::for_each(graphfile, lambda);
			graph = mergeGraphs(parts);
		}
		else{
			std::cout << "No graph file exists" << std::endl;
			std::exit(0);
		}
	}

	std::string filename {argv[2]};

	if (!is_file_exist(argv[2])){
		std::cout << "No input file exists" << std::endl;
		std::exit(0);
	}

	if (endsWith(filename, ".gam"))
	{
		std::vector<vg::Alignment> alignments;
		std::cout << "load alignments from " << argv[2] << std::endl;
		std::ifstream seedfile { argv[2], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> alignmentLambda = [&alignments](vg::Alignment& a) {
			alignments.push_back(a);
		};
		stream::for_each(seedfile, alignmentLambda);
		order(alignments, [](const vg::Alignment& alignment) {return alignment.path().mapping(0).position().node_id();}, graph);
		std::ofstream alignmentOut { argv[3], std::ios::out | std::ios::binary };
		std::cout << "write alignments to " << argv[3] << std::endl;
		stream::write_buffered(alignmentOut, alignments, 0);
	}
	else if (endsWith(filename, ".snarls"))
	{
		std::vector<vg::Snarl> alignments;
		std::cout << "load snarls from " << argv[2] << std::endl;
		std::ifstream seedfile { argv[2], std::ios::in | std::ios::binary };
		std::function<void(vg::Snarl&)> alignmentLambda = [&alignments](vg::Snarl& a) {
			alignments.push_back(a);
		};
		stream::for_each(seedfile, alignmentLambda);
		order(alignments, [](const vg::Snarl& snarl) {return snarl.start().node_id();}, graph);
		std::ofstream alignmentOut { argv[3], std::ios::out | std::ios::binary };
		std::cout << "write snarls to " << argv[3] << std::endl;
		stream::write_buffered(alignmentOut, alignments, 0);
	}
	else if (endsWith(filename, ".trans"))
	{
		std::vector<vg::SnarlTraversal> alignments;
		std::cout << "load trans from " << argv[2] << std::endl;
		std::ifstream seedfile { argv[2], std::ios::in | std::ios::binary };
		std::function<void(vg::SnarlTraversal&)> alignmentLambda = [&alignments](vg::SnarlTraversal& a) {
			alignments.push_back(a);
		};
		stream::for_each(seedfile, alignmentLambda);
		order(alignments, [](const vg::SnarlTraversal& traversal) {return traversal.snarl().start().node_id();}, graph);
		std::ofstream alignmentOut { argv[3], std::ios::out | std::ios::binary };
		std::cout << "write trans to " << argv[3] << std::endl;
		stream::write_buffered(alignmentOut, alignments, 0);
	}
}