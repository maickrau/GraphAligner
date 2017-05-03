#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include "vg.pb.h"
#include "gssw.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "TopologicalSort.h"

class GraphMappingContainer
{
public:
	GraphMappingContainer(gssw_graph_mapping* ptr) : ptr(ptr) {};
	GraphMappingContainer(const GraphMappingContainer& second) = delete;
	GraphMappingContainer& operator=(const GraphMappingContainer& second) = delete;
	GraphMappingContainer(GraphMappingContainer&& second)
	{
		ptr = second.ptr;
		second.ptr = nullptr;
	};
	GraphMappingContainer& operator=(GraphMappingContainer&& second)
	{
		if (&second == this) return *this;
		unload();
		ptr = second.ptr;
		second.ptr = nullptr;
	};
	~GraphMappingContainer()
	{
		unload();
	};
	operator gssw_graph_mapping*()
	{
		return ptr;
	};
	std::string seq_id;
	gssw_graph_mapping* ptr;
private:
	void unload()
	{
		if (ptr != nullptr) gssw_graph_mapping_destroy(ptr);
	};
};

vg::Alignment gsswToVgMapping(GraphMappingContainer& mapping)
{
	vg::Alignment result;
	result.set_name(mapping.seq_id);
	auto path = new vg::Path;
	result.set_allocated_path(path);
	for (int node = 0; node < mapping.ptr->cigar.length; node++)
	{
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(node+1);
		position->set_node_id(mapping.ptr->cigar.elements[node].node->id);
	}
	return result;
}

vg::Graph mergeGraphs(std::vector<vg::Graph> parts)
{
	std::vector<vg::Node> allNodes;
	std::vector<vg::Edge> allEdges;
	for (size_t i = 0; i < parts.size(); i++)
	{
		for (int j = 0; j < parts[i].node_size(); j++)
		{
			allNodes.push_back(parts[i].node(j));
		}
		for (int j = 0; j < parts[i].edge_size(); j++)
		{
			allEdges.push_back(parts[i].edge(j));
		}
	}
	vg::Graph newGraph;
	for (size_t i = 0; i < allNodes.size(); i++)
	{
		auto node = newGraph.add_node();
		node->set_id(allNodes[i].id());
		node->set_sequence(allNodes[i].sequence());
		node->set_name(allNodes[i].name());
	}
	for (size_t i = 0; i < allEdges.size(); i++)
	{
		auto edge = newGraph.add_edge();
		edge->set_from(allEdges[i].from());
		edge->set_to(allEdges[i].to());
		edge->set_from_start(allEdges[i].from_start());
		edge->set_to_end(allEdges[i].to_end());
		edge->set_overlap(allEdges[i].overlap());
	}
	return newGraph;
}


std::vector<GraphMappingContainer> getOptimalPinnedMappings(const vg::Graph& vggraph, const std::vector<FastQ>& reads)
{
	//code mostly from gssw's example.c
	int8_t match = 1, mismatch = 4;
	uint8_t gap_open = 1, gap_extension = 1;
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(match, mismatch);
	std::vector<gssw_node*> gsswnodes;
	std::vector<size_t> order = topologicalSort(vggraph);
	std::vector<vg::Node> nodesToEnter;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		std::cerr << "before sorting: node index " << i << " id " << vggraph.node(i).id() << std::endl;
	}
	for (int i = 0; i< vggraph.node_size(); i++)
	{
		nodesToEnter.push_back(vggraph.node(order[i]));
	}
	gssw_graph* graph = gssw_graph_create(vggraph.node_size());
	//todo check: do these need to be in this order? create nodes, create edges and only then insert to graph
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		gsswnodes.push_back((gssw_node*)gssw_node_create((void*)"", nodesToEnter[i].id(), nodesToEnter[i].sequence().c_str(), nt_table, mat));
	}
	std::map<size_t, int> ids;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		ids[nodesToEnter[i].id()] = i;
		std::cerr << "node index " << i << " id " << nodesToEnter[i].id() << std::endl;
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		gssw_nodes_add_edge(gsswnodes[ids[vggraph.edge(i).from()]], gsswnodes[ids[vggraph.edge(i).to()]]);
	}
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		gssw_graph_add_node(graph, gsswnodes[i]);
	}

	std::vector<GraphMappingContainer> result;
	for (size_t i = 0; i < reads.size(); i++)
	{
		gssw_graph_fill(graph, reads[i].sequence.c_str(), nt_table, mat, gap_open, gap_extension, 0, 0, 15, 2);
		gssw_graph_mapping* gmp = gssw_graph_trace_back (graph,
			reads[i].sequence.c_str(),
			reads[i].sequence.size(),
			nt_table,
			mat,
			gap_open,
			gap_extension,
			0, 0);
		// gssw_graph_mapping* gmp = gssw_graph_trace_back_pinned_qual_adj (graph,
		// 	gsswnodes.back(),
		// 	reads[i].sequence.c_str(),
		// 	reads[i].quality.c_str(),
		// 	reads[i].sequence.size(),
		// 	nt_table,
		// 	mat,
		// 	gap_open,
		// 	gap_extension,
		// 	0, 0);

		result.emplace_back(gmp);
		result.back().seq_id = reads[i].seq_id;
	}

    //todo: does the graph need to exist to use the graph mapping?
	// gssw_graph_destroy(graph);

	return result;
}

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	vg::Graph graph;
	std::cerr << "load graph from " << argv[1] << std::endl;
	std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
	std::vector<vg::Graph> parts;
	std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
		parts.push_back(g);
	};
	stream::for_each(graphfile, lambda);

	graph = mergeGraphs(parts);

	auto fastqs = loadFastqFromFile(argv[2]);
	std::cout << fastqs.size() << " reads" << std::endl;
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		std::cout << fastqs[i].sequence << std::endl;
	}
	auto result = getOptimalPinnedMappings(graph, fastqs);
	std::vector<vg::Alignment> alignments;
	for (size_t i = 0; i < result.size(); i++)
	{
		alignments.push_back(gsswToVgMapping(result[i]));
		gssw_print_graph_mapping(result[i], stdout);
	}

	std::ofstream alignmentOut { argv[3], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, alignments, 0);

	return 0;
}
