#include <iostream>
#include <fstream>
#include "vg.pb.h"
#include "gssw.h"

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
private:
	gssw_graph_mapping* ptr;
	void unload()
	{
		if (ptr != nullptr) gssw_graph_mapping_destroy(ptr);
	};
};

std::vector<GraphMappingContainer> getOptimalPinnedMappings(const vg::Graph& vggraph, const std::vector<std::string>& reads)
{
	//code mostly from gssw's example.c
	int8_t match = 1, mismatch = 4;
	uint8_t gap_open = 6, gap_extension = 1;
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(match, mismatch);
	std::vector<gssw_node*> gsswnodes;
	gssw_graph* graph = gssw_graph_create(vggraph.node_size());
	//todo check: do these need to be in this order? create nodes, create edges and only then insert to graph
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		gsswnodes.push_back((gssw_node*)gssw_node_create((void*)"", i, vggraph.node(i).sequence().c_str(), nt_table, mat));
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		gssw_nodes_add_edge(gsswnodes[vggraph.edge(i).from()], gsswnodes[vggraph.edge(i).to()]);
	}
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		gssw_graph_add_node(graph, gsswnodes[i]);
	}

	std::vector<GraphMappingContainer> result;
	for (int i = 0; i < reads.size(); i++)
	{

		gssw_graph_fill(graph, reads[i].c_str(), nt_table, mat, gap_open, gap_extension, 0, 0, 15, 2);
		gssw_graph_mapping* gmp = gssw_graph_trace_back_pinned (graph,
			gsswnodes.back(),
			reads[i].c_str(),
			reads[i].size(),
			nt_table,
			mat,
			gap_open,
			gap_extension,
			0, 0);

		result.emplace_back(gmp);
	}

    //todo: does the graph need to exist to use the graph mapping?
	gssw_graph_destroy(graph);

	return result;
}

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;
	vg::Graph graph;
	{
		std::cerr << "load graph from " << argv[1] << std::endl;
		std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
		if (!graph.ParseFromIstream(&graphfile))
		{
			std::cerr << "Graph load failed." << std::endl;
			return 1;
		}
	}
	std::vector<std::string> reads;
	for (int i = 2; i < argc; i++)
	{
		reads.emplace_back(argv[i]);
	}
	auto result = getOptimalPinnedMappings(graph, reads);
	return 0;
}
