#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include "vg.pb.h"
#include "gssw.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "TopologicalSort.h"
#include "SubgraphFromSeed.h"

size_t GraphSizeInBp(const vg::Graph& graph)
{
	size_t result = 0;
	for (int i = 0; i < graph.node_size(); i++)
	{
		result += graph.node(i).sequence().size();
	}
	return result;
}

vg::Alignment gsswToVgMapping(gssw_graph_mapping* mapping, std::string seq_id, bool reverse)
{
	vg::Alignment result;
	result.set_name(seq_id);
	auto path = new vg::Path;
	result.set_allocated_path(path);
	for (size_t node = 0; node < mapping->cigar.length; node++)
	{
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(node+1);
		position->set_node_id(mapping->cigar.elements[node].node->id);
		if (reverse) position->set_is_reverse(true);
	}
	result.set_score(mapping->score);
	return result;
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

vg::Alignment getOptimalPinnedMapping(const vg::Graph& vggraph, const FastQ& read)
{
	//code mostly from gssw's example.c
	int8_t match = 1, mismatch = 4;
	uint8_t gap_open = 1, gap_extension = 1;
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(match, mismatch);
	std::vector<gssw_node*> gsswnodes;
	std::vector<size_t> order = topologicalSort(vggraph);
	std::vector<const vg::Node*> nodesToEnter;
	for (int i = 0; i< vggraph.node_size(); i++)
	{
		nodesToEnter.push_back(&vggraph.node(order[i]));
	}
	gssw_graph* graph = gssw_graph_create(vggraph.node_size());
	//todo check: do these need to be in this order? create nodes, create edges and only then insert to graph
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		gsswnodes.push_back((gssw_node*)gssw_node_create((void*)"", nodesToEnter[i]->id(), nodesToEnter[i]->sequence().c_str(), nt_table, mat));
	}
	std::map<size_t, int> ids;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		ids[nodesToEnter[i]->id()] = i;
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		gssw_nodes_add_edge(gsswnodes[ids[vggraph.edge(i).from()]], gsswnodes[ids[vggraph.edge(i).to()]]);
	}
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		gssw_graph_add_node(graph, gsswnodes[i]);
	}

	std::cout << "read " << read.seq_id << std::endl;
	std::cout << "read size " << read.sequence.size() << " bp, subgraph size " << GraphSizeInBp(vggraph) << " bp" << std::endl;
	std::cout << "align forward" << std::endl;
	gssw_graph_fill(graph, read.sequence.c_str(), nt_table, mat, gap_open, gap_extension, 5, 5, 15, 2);
	gssw_graph_mapping* gmpForward = gssw_graph_trace_back (graph,
		read.sequence.c_str(),
		read.sequence.size(),
		nt_table,
		mat,
		gap_open,
		gap_extension,
		0, 0);

	std::cout << "align backwards" << std::endl;
	auto reverseComplement = read.reverseComplement();
	gssw_graph_fill(graph, reverseComplement.sequence.c_str(), nt_table, mat, gap_open, gap_extension, 5, 5, 15, 2);
	gssw_graph_mapping* gmpBackwards = gssw_graph_trace_back (graph,
		reverseComplement.sequence.c_str(),
		reverseComplement.sequence.size(),
		nt_table,
		mat,
		gap_open,
		gap_extension,
		0, 0);

	gssw_graph_mapping* bestMapping;
	bool reverse = false;

	if (gmpForward->score > gmpBackwards->score)
	{
		bestMapping = gmpForward;
	}
	else
	{
		bestMapping = gmpBackwards;
		reverse = true;
	}

	vg::Alignment vgResult = gsswToVgMapping(bestMapping, read.seq_id, reverse);

	gssw_graph_mapping_destroy(gmpForward);
	gssw_graph_mapping_destroy(gmpBackwards);

	gssw_graph_destroy(graph);

	free(nt_table);
	free(mat);

	return vgResult;
}

void runMappings(const vg::Graph& graph, const std::vector<std::pair<FastQ*, vg::Alignment>>& fastqAlignmentPairs, std::map<FastQ*, vg::Alignment>& alignments, int threadnum)
{
	std::string msg;
	for (size_t i = 0; i < fastqAlignmentPairs.size(); i++)
	{
		{
			std::ostringstream writer;
			writer << "thread " << threadnum << " " << i << "/" << fastqAlignmentPairs.size() << std::endl;
			msg = writer.str();
		}
		std::cerr << msg;
		auto fastq = fastqAlignmentPairs[i].first;
		auto alignment = fastqAlignmentPairs[i].second;
		auto seedGraph = ExtractSubgraph(graph, alignment, fastq->sequence.size());
		auto bestMapping = getOptimalPinnedMapping(seedGraph, *fastq);
		if (bestMapping.score() > -1)
		{
			if (alignments.count(fastq) == 0)
			{
				alignments[fastq] = bestMapping;
				{
					std::ostringstream writer;
					writer << "thread " << threadnum << " successfully aligned read " << fastq->seq_id << std::endl;
					msg = writer.str();
				}
				std::cerr << msg;
				continue;
			}
			auto oldScore = alignments[fastq];
			if (bestMapping.score() > oldScore.score())
			{
				alignments[fastq] = bestMapping;
			}
		}
	}
	{
		std::ostringstream writer;
		writer << "thread " << threadnum << " finished with " << alignments.size() << " alignments" << std::endl;
		msg = writer.str();
	}
	std::cerr << msg;
}

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	vg::Graph graph;
	{
		std::cout << "load graph from " << argv[1] << std::endl;
		std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
		std::vector<vg::Graph> parts;
		std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
			parts.push_back(g);
		};
		stream::for_each(graphfile, lambda);
		graph = mergeGraphs(parts);
		std::cout << "graph is " << GraphSizeInBp(graph) << " bp large" << std::endl;
	}

	auto fastqs = loadFastqFromFile(argv[2]);
	std::cout << fastqs.size() << " reads" << std::endl;

	std::map<std::string, std::vector<vg::Alignment>> seeds;
	{
		std::ifstream seedfile { argv[3], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> alignmentLambda = [&seeds](vg::Alignment& a) {
			seeds[a.name()].push_back(a);
		};
		stream::for_each(seedfile, alignmentLambda);
	}

	int numThreads = std::stoi(argv[5]);
	std::vector<std::vector<std::pair<FastQ*, vg::Alignment>>> readAlignmentPairsPerThread;
	std::vector<std::map<FastQ*, vg::Alignment>> alignmentsPerThread;
	readAlignmentPairsPerThread.resize(numThreads);
	alignmentsPerThread.resize(numThreads);
	int currentThread = 0;
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		if (seeds.count(fastqs[i].seq_id) == 0) continue;
		for (size_t j = 0; j < seeds[fastqs[i].seq_id].size(); j++)
		{
			readAlignmentPairsPerThread[currentThread].emplace_back(&fastqs[i], seeds[fastqs[i].seq_id][j]);
			currentThread++;
			currentThread %= numThreads;
		}
	}

	std::vector<std::thread> threads;

	for (int i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&graph, &readAlignmentPairsPerThread, &alignmentsPerThread, i]() { runMappings(graph, readAlignmentPairsPerThread[i], alignmentsPerThread[i], i); });
	}

	for (int i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}

	std::vector<vg::Alignment> alignments;

	for (size_t i = 0; i < fastqs.size(); i++)
	{
		vg::Alignment* bestAlignment = nullptr;
		for (int j = 0; j < numThreads; j++)
		{
			if (alignmentsPerThread[j].count(&fastqs[i]) == 0) continue;
			if (bestAlignment == nullptr) 
			{
				bestAlignment = &alignmentsPerThread[j][&fastqs[i]];
				continue;
			}
			if (alignmentsPerThread[j][&fastqs[i]].score() > bestAlignment->score())
			{
				bestAlignment = &alignmentsPerThread[j][&fastqs[i]];
				continue;
			}
		}
		if (bestAlignment != nullptr) alignments.push_back(*bestAlignment);
	}

	std::cerr << "final result has " << alignments.size() << " alignments" << std::endl;

	std::ofstream alignmentOut { argv[4], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, alignments, 0);

	return 0;
}
