#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "TopologicalSort.h"
#include "SubgraphFromSeed.h"
#include "GraphAligner.h"
#include "mfvs_graph.h"

size_t GraphSizeInBp(const vg::Graph& graph)
{
	size_t result = 0;
	for (int i = 0; i < graph.node_size(); i++)
	{
		result += graph.node(i).sequence().size();
	}
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

vg::Alignment getOptimalMapping(const vg::Graph& graph, const FastQ& read)
{
	vg::Alignment betterAlignment;
	{
		std::cerr << "forward" << std::endl;
		GraphAligner<uint32_t, int32_t> forwardAlignment;
		for (int i = 0; i < graph.node_size(); i++)
		{
			forwardAlignment.AddNode(graph.node(i).id(), graph.node(i).sequence());
		}
		for (int i = 0; i < graph.edge_size(); i++)
		{
			forwardAlignment.AddEdge(graph.edge(i).from(), graph.edge(i).to());
		}
		forwardAlignment.Finalize();
		betterAlignment = forwardAlignment.AlignOneWay(read.seq_id, read.sequence, false);
	}
	{
		std::cerr << "backward" << std::endl;
		GraphAligner<uint32_t, int32_t> backwardAlignment;
		for (int i = graph.node_size()-1; i >= 0; i--)
		{
			backwardAlignment.AddNode(graph.node(i).id(), FastQ::reverseComplement(graph.node(i).sequence()));
		}
		for (int i = 0; i < graph.edge_size(); i++)
		{
			backwardAlignment.AddEdge(graph.edge(i).to(), graph.edge(i).from());
		}
		backwardAlignment.Finalize();
		auto backwards = backwardAlignment.AlignOneWay(read.seq_id, read.sequence, true);
		if (backwards.score() > betterAlignment.score()) betterAlignment = backwards;
	}
	return betterAlignment;
}

int numberOfVerticesOutOfOrder(const vg::Graph& vggraph)
{
	std::map<int, int> ids;
	int result = 0;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		ids[vggraph.node(i).id()] = i;
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		if (ids[vggraph.edge(i).to()] < ids[vggraph.edge(i).from()]) result++;
	}
	return result;
}

vg::Graph OrderByFeedbackVertexset(const vg::Graph& vggraph)
{
	mfvs::Graph mfvsgraph { vggraph.node_size() };
	std::map<int, int> ids;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		mfvsgraph.addVertex(i);
		ids[vggraph.node(i).id()] = i;
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		mfvsgraph.addEdge(vggraph.edge(i).from(), vggraph.edge(i).to());
	}
	auto vertexSetvector = mfvsgraph.minimumFeedbackVertexSet();
	std::set<int> vertexset { vertexSetvector.begin(), vertexSetvector.end() };
	std::cout << "feedback vertex set size: " << vertexSetvector.size() << std::endl;
	vg::Graph graphWithoutVertexset;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		if (vertexset.count(i) > 0) continue;
		auto node = graphWithoutVertexset.add_node();
		node->set_id(vggraph.node(i).id());
		node->set_sequence(vggraph.node(i).sequence());
		node->set_name(vggraph.node(i).name());
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		if (vertexset.count(ids[vggraph.edge(i).from()]) > 0 || vertexset.count(ids[vggraph.edge(i).to()]) > 0) continue;
		auto edge = graphWithoutVertexset.add_edge();
		edge->set_from(vggraph.edge(i).from());
		edge->set_to(vggraph.edge(i).to());
		edge->set_from_start(vggraph.edge(i).from_start());
		edge->set_to_end(vggraph.edge(i).to_end());
		edge->set_overlap(vggraph.edge(i).overlap());
	}
	auto order = topologicalSort(graphWithoutVertexset);
	vg::Graph resultGraph;
	for (size_t i = 0; i < order.size(); i++)
	{
		auto node = resultGraph.add_node();
		node->set_id(vggraph.node(order[i]).id());
		node->set_sequence(vggraph.node(order[i]).sequence());
		node->set_name(vggraph.node(order[i]).name());
	}
	for (size_t i = 0; i < vertexSetvector.size(); i++)
	{
		auto node = resultGraph.add_node();
		node->set_id(vggraph.node(vertexSetvector[i]).id());
		node->set_sequence(vggraph.node(vertexSetvector[i]).sequence());
		node->set_name(vggraph.node(vertexSetvector[i]).name());
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		auto edge = resultGraph.add_edge();
		edge->set_from(vggraph.edge(i).from());
		edge->set_to(vggraph.edge(i).to());
		edge->set_from_start(vggraph.edge(i).from_start());
		edge->set_to_end(vggraph.edge(i).to_end());
		edge->set_overlap(vggraph.edge(i).overlap());
	}
	return resultGraph;
}

void outputGraph(std::string filename, const vg::Graph& graph)
{
	std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {graph};
	stream::write_buffered(alignmentOut, writeVector, 0);
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
		auto seedGraphUnordered = ExtractSubgraph(graph, alignment, fastq->sequence.size());
		std::cout << "graph size " << seedGraphUnordered.node_size() << " nodes ";
		std::cout << "out of order before sorting: " << numberOfVerticesOutOfOrder(seedGraphUnordered) << std::endl;
		auto seedGraph = OrderByFeedbackVertexset(seedGraphUnordered);
		std::cout << "graph size " << seedGraph.node_size() << " nodes ";
		std::cout << "out of order after sorting: " << numberOfVerticesOutOfOrder(seedGraph) << std::endl;
		std::cout << "align " << fastq->sequence.size() << " bp read to " << GraphSizeInBp(seedGraph) << " bp graph" << std::endl;
		outputGraph("outgraph.gam", seedGraph);
		auto bestMapping = getOptimalMapping(seedGraph, *fastq);
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
