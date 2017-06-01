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
#include "BigraphToDigraph.h"


class BufferedWriter : std::ostream
{
public:
	class FlushClass {};
	BufferedWriter(std::ostream& stream) : stream(stream) {};
	template <typename T>
	BufferedWriter& operator<<(T obj)
	{
		stringstream << obj;
		return *this;
	}
	BufferedWriter& operator<<(FlushClass f)
	{
		flush();
		return *this;
	}
	void flush()
	{
		stringstream << std::endl;
		stream << stringstream.str();
		stringstream.str("");
	}
	static FlushClass Flush;
private:
	std::ostream& stream;
	std::stringstream stringstream;
};

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

size_t GraphSizeInBp(const DirectedGraph& graph)
{
	size_t result = 0;
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		result += graph.nodes[i].sequence.size();
	}
	return result;
}

size_t GraphSizeInBp(const vg::Graph& graph)
{
	size_t result = 0;
	for (int i = 0; i < graph.node_size(); i++)
	{
		result += graph.node(i).sequence().size();
	}
	return result;
}

// augment base VG graph with alignments by embedding alignment paths
vg::Graph augmentGraphwithAlignment(const vg::Graph& graph, const std::vector<vg::Alignment> alignments)
{
	vg::Graph augmentedGraph;
	std::vector<const vg::Node*> allNodes;
	std::vector<const vg::Edge*> allEdges;
	
	for (int j = 0; j < graph.node_size(); j++)
	{
		allNodes.push_back(&graph.node(j));
	}
	for (int j = 0; j < graph.edge_size(); j++)
	{
		allEdges.push_back(&graph.edge(j));
	}
	for (size_t i = 0; i < allNodes.size(); i++)
	{
		auto node = augmentedGraph.add_node();
		node->set_id(allNodes[i]->id());
		node->set_sequence(allNodes[i]->sequence());
		node->set_name(allNodes[i]->name());
	}
	for (size_t i = 0; i < allEdges.size(); i++)
	{
		auto edge = augmentedGraph.add_edge();
		edge->set_from(allEdges[i]->from());
		edge->set_to(allEdges[i]->to());
		edge->set_from_start(allEdges[i]->from_start());
		edge->set_to_end(allEdges[i]->to_end());
		edge->set_overlap(allEdges[i]->overlap());
	}
	
	for(int k=0; k < alignments.size()-1; k++)
	{
		for (int i = 0; i < alignments[k].path().mapping_size()-1; i++)
		{
			if (alignments[k].path().mapping(i).position().is_reverse()){
				auto edge = augmentedGraph.add_edge();
				edge->set_from(alignments[k].path().mapping(i+1).position().node_id());
				edge->set_to(alignments[k].path().mapping(i).position().node_id());
				edge->set_from_start(1);
				edge->set_to_end(1);
				edge->set_overlap(0);
			}else{
				auto edge = augmentedGraph.add_edge();
				edge->set_from(alignments[k].path().mapping(i).position().node_id());
				edge->set_to(alignments[k].path().mapping(i+1).position().node_id());
				edge->set_from_start(0);
				edge->set_to_end(0);
				edge->set_overlap(0);  
			}
			i++;
		}
	}
	return augmentedGraph;
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

int numberOfVerticesOutOfOrder(const DirectedGraph& digraph)
{
	std::map<int, int> ids;
	std::set<int> outOfOrderSet;
	for (size_t i = 0; i < digraph.edges.size(); i++)
	{
		if (digraph.edges[i].toIndex <= digraph.edges[i].fromIndex) outOfOrderSet.insert(digraph.edges[i].toIndex);
	}
	return outOfOrderSet.size();
}

void OrderByFeedbackVertexset(DirectedGraph& graph)
{
	BufferedWriter output { std::cout };
	mfvs::Graph mfvsgraph { graph.nodes.size() };
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		mfvsgraph.addVertex(i);
	}
	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		mfvsgraph.addEdge(graph.edges[i].fromIndex, graph.edges[i].toIndex);
	}
	output << "Before the removal of MVFS, is the graph acyclic:" << mfvsgraph.isAcyclic() << "\n";
	auto vertexSetvector = mfvsgraph.minimumFeedbackVertexSet();
	std::set<int> vertexset { vertexSetvector.begin(), vertexSetvector.end() };
	output << "feedback vertex set size: " << vertexSetvector.size() << "\n";
	output << "After the removal of MVFS, is the graph acyclic:" << mfvsgraph.isAcyclic() << BufferedWriter::Flush;
	DirectedGraph graphWithoutVFS { graph };
	graphWithoutVFS.RemoveNodes(vertexset);
	std::vector<size_t> indexOrder = topologicalSort(graphWithoutVFS);
	std::vector<int> nodeIdOrder;
	for (size_t i = 0; i < vertexSetvector.size(); i++)
	{
		nodeIdOrder.push_back(graph.nodes[vertexSetvector[i]].nodeId);
	}
	for (size_t i = 0; i < indexOrder.size(); i++)
	{
		nodeIdOrder.push_back(graphWithoutVFS.nodes[indexOrder[i]].nodeId);
	}
	graph.ReorderByNodeIds(nodeIdOrder);
}

void outputGraph(std::string filename, const vg::Graph& graph)
{
	std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {graph};
	stream::write_buffered(alignmentOut, writeVector, 0);
}

std::vector<int> getSinkNodes(const DirectedGraph& graph)
{
	std::set<int> notSinkNodes;
	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		notSinkNodes.insert(graph.nodes[graph.edges[i].fromIndex].nodeId);
	}
	std::vector<int> result;
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		if (notSinkNodes.count(graph.nodes[i].nodeId) == 0) result.push_back(graph.nodes[i].nodeId);
	}
	return result;
}

std::vector<int> getSourceNodes(const DirectedGraph& graph)
{
	std::set<int> notSourceNodes;
	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		notSourceNodes.insert(graph.nodes[graph.edges[i].toIndex].nodeId);
	}
	std::vector<int> result;
	for (size_t i = 0; i < graph.nodes.size(); i++)
	{
		if (notSourceNodes.count(graph.nodes[i].nodeId) == 0) result.push_back(graph.nodes[i].nodeId);
	}
	return result;
}

bool GraphEqual(const DirectedGraph& first, const DirectedGraph& second)
{
	if (first.nodes.size() != second.nodes.size()) return false;
	if (first.edges.size() != second.edges.size()) return false;
	for (size_t i = 0; i < first.nodes.size(); i++)
	{
		if (first.nodes[i].nodeId != second.nodes[i].nodeId) return false;
		if (first.nodes[i].originalNodeId != second.nodes[i].originalNodeId) return false;
		if (first.nodes[i].sequence != second.nodes[i].sequence) return false;
		if (first.nodes[i].rightEnd != second.nodes[i].rightEnd) return false;
	}
	for (size_t i = 0; i < first.edges.size(); i++)
	{
		if (first.edges[i].fromIndex != second.edges[i].fromIndex) return false;
		if (first.edges[i].toIndex != second.edges[i].toIndex) return false;
	}
	return true;
}

void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const DirectedGraph& graph)
{
	std::map<int, int> idMapper;
	for (size_t j = 0; j < graph.nodes.size(); j++)
	{
		if (idMapper.count(graph.nodes[j].nodeId) > 0 && idMapper[graph.nodes[j].nodeId] != graph.nodes[j].originalNodeId)
		{
			assert(idMapper.count(graph.nodes[j].nodeId) == 0 || idMapper[graph.nodes[j].nodeId] == graph.nodes[j].originalNodeId);
			idMapper[graph.nodes[j].nodeId] = graph.nodes[j].originalNodeId;
		}
		assert(idMapper.count(graph.nodes[j].nodeId) == 0 || idMapper[graph.nodes[j].nodeId] == graph.nodes[j].originalNodeId);
		idMapper[graph.nodes[j].nodeId] = graph.nodes[j].originalNodeId;
	}
	for (int i = 0; i < alignment.path().mapping_size(); i++)
	{
		int digraphNodeId = alignment.path().mapping(i).position().node_id();
		assert(idMapper.count(digraphNodeId) > 0);
		int originalNodeId = idMapper[digraphNodeId];
		alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_node_id(originalNodeId);
	}
}

void runComponentMappings(const vg::Graph& graph, const std::vector<const FastQ*>& fastQs, const std::map<const FastQ*, std::vector<vg::Alignment>>& seedhits, std::vector<vg::Alignment>& alignments, int threadnum, int bandwidth)
{
	BufferedWriter cerroutput {std::cerr};
	BufferedWriter coutoutput {std::cout};
	for (size_t i = 0; i < fastQs.size(); i++)
	{
		const FastQ* fastq = fastQs[i];
		coutoutput << "thread " << threadnum << " " << i << "/" << fastQs.size() << "\n";
		coutoutput << "read " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" << "\n";
		coutoutput << "components: " << seedhits.at(fastq).size() << BufferedWriter::Flush;
		if (seedhits.at(fastq).size() == 0)
		{
			cerroutput << "read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
			continue;
		}
		std::vector<std::tuple<int, DirectedGraph>> components;
		std::vector<std::pair<int, size_t>> graphAlignerSeedHits;
		for (size_t j = 0; j < seedhits.at(fastq).size(); j++)
		{
			auto& seedhit = seedhits.at(fastq)[j];
			//add seed hit to both strands
			if (seedhit.path().mapping(0).position().node_id() > 0)
			{
				graphAlignerSeedHits.emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.query_position());
			}
			coutoutput << "thread " << threadnum << " read " << i << " component " << j << "/" << seedhits.at(fastq).size() << "\n" << BufferedWriter::Flush;
			auto seedGraphUnordered = ExtractSubgraph(graph, seedhit, fastq->sequence.size()*2);
			DirectedGraph seedGraph {seedGraphUnordered};
			coutoutput << "component size " << GraphSizeInBp(seedGraph) << "bp" << BufferedWriter::Flush;
			coutoutput << "component out of order before sorting: " << numberOfVerticesOutOfOrder(seedGraph) << BufferedWriter::Flush;
			int startpos = 0;
			OrderByFeedbackVertexset(seedGraph);
			bool alreadyIn = false;
			startpos = seedhit.query_position();
			for (size_t k = 0; k < components.size(); k++)
			{
				if (startpos == std::get<0>(components[k]) && GraphEqual(seedGraph,std::get<1>(components[k])))
				{
					coutoutput << "already exists" << BufferedWriter::Flush;
					alreadyIn = true;
					break;
				}
			}
			if (alreadyIn) continue;
			coutoutput << "component out of order after sorting: " << numberOfVerticesOutOfOrder(seedGraph) << BufferedWriter::Flush;
			components.emplace_back(startpos, seedGraph);
			coutoutput << "component position: " << startpos << " " << BufferedWriter::Flush;
		}
		std::sort(components.begin(), components.end(), [](auto& left, auto& right) { return std::get<0>(left) < std::get<0>(right); });

		DirectedGraph augmentedGraph;
		std::vector<std::vector<int>> sources;
		std::vector<std::vector<int>> sinks;
		for (size_t i = 0; i < components.size(); i++)
		{
			sources.emplace_back(getSourceNodes(std::get<1>(components[i])));
			sinks.emplace_back(getSinkNodes(std::get<1>(components[i])));
			augmentedGraph.AddSubgraph(std::get<1>(components[i]));
		}
		for (size_t i = 0; i < components.size(); i++)
		{
			for (size_t j = i+1; j < components.size(); j++)
			{
				augmentedGraph.ConnectComponents(sinks[i], sources[j]);
			}
		}
		coutoutput << "thread " << threadnum << " augmented graph is " << GraphSizeInBp(augmentedGraph) << "bp" << BufferedWriter::Flush;
		auto augmentedGraphSeedHits = augmentedGraph.GetSeedHits(fastq->sequence, graphAlignerSeedHits);
		std::vector<int> augmentedGraphReachable;
		for (size_t i = 0; i < augmentedGraphSeedHits.size(); i++)
		{
			augmentedGraphReachable.push_back(augmentedGraphSeedHits[i].nodeId);
		}
		augmentedGraph.PruneByReachability(augmentedGraphReachable);
		coutoutput << "thread " << threadnum << " augmented graph after pruning is " << GraphSizeInBp(augmentedGraph) << "bp" << BufferedWriter::Flush;
		coutoutput << "augmented graph out of order before sorting: " << numberOfVerticesOutOfOrder(augmentedGraph) << "\n";
		if (augmentedGraph.nodes.size() == 0) continue;
		OrderByFeedbackVertexset(augmentedGraph);
		coutoutput << "augmented graph out of order after sorting: " << numberOfVerticesOutOfOrder(augmentedGraph) << BufferedWriter::Flush;

		GraphAligner<uint32_t, int32_t> augmentedGraphAlignment;
		for (size_t j = 0; j < augmentedGraph.nodes.size(); j++)
		{
			augmentedGraphAlignment.AddNode(augmentedGraph.nodes[j].nodeId, augmentedGraph.nodes[j].sequence);
			// std::cerr << "node: " << augmentedGraph.nodes[j].nodeId << std::endl;
		}
		for (size_t j = 0; j < augmentedGraph.edges.size(); j++)
		{
			augmentedGraphAlignment.AddEdgeNodeId(augmentedGraph.nodes[augmentedGraph.edges[j].fromIndex].nodeId, augmentedGraph.nodes[augmentedGraph.edges[j].toIndex].nodeId);
			// std::cerr << "edge: " << augmentedGraph.nodes[augmentedGraph.edges[j].fromIndex].nodeId << " -> " << augmentedGraph.nodes[augmentedGraph.edges[j].toIndex].nodeId << std::endl;
		}
		augmentedGraphAlignment.Finalize();

		std::vector<GraphAligner<uint32_t, int32_t>::SeedHit> alignerSeedHits;
		for (size_t i = 0; i < augmentedGraphSeedHits.size(); i++)
		{
			alignerSeedHits.emplace_back(augmentedGraphSeedHits[i].seqPos, augmentedGraphSeedHits[i].nodeId, augmentedGraphSeedHits[i].nodePos);
		}

		auto alignment = augmentedGraphAlignment.AlignOneWay(fastQs[i]->seq_id, fastQs[i]->sequence, false, bandwidth, alignerSeedHits);

		//failed alignment, don't output
		if (alignment.alignmentFailed)
		{
			cerroutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			continue;
		}
		if (alignment.alignment.score() == std::numeric_limits<decltype(alignment.alignment.score())>::min())
		{
			cerroutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			continue;
		}

		coutoutput << "read " << fastq->seq_id << " max distance from band " << alignment.maxDistanceFromBand << BufferedWriter::Flush;
		coutoutput << "read " << fastq->seq_id << " score " << alignment.alignment.score() << BufferedWriter::Flush;
		if (alignment.maxDistanceFromBand > bandwidth * 0.66)
		{
			cerroutput << "read " << fastq->seq_id << " max distance from band is high: " << alignment.maxDistanceFromBand << BufferedWriter::Flush;
		}
		if (alignment.alignment.score() < fastq->sequence.size() * 0.7)
		{
			cerroutput << "read " << fastq->seq_id << " score is low: " << alignment.alignment.score() << BufferedWriter::Flush;
		}

		replaceDigraphNodeIdsWithOriginalNodeIds(alignment.alignment, augmentedGraph);

		alignments.push_back(alignment.alignment);
		coutoutput << "thread " << threadnum << " successfully aligned read " << fastq->seq_id << BufferedWriter::Flush;
		std::vector<vg::Alignment> alignmentvec;
		alignmentvec.emplace_back(alignments.back());
		std::string filename;
		filename = "alignment_";
		filename += std::to_string(threadnum);
		filename += "_";
		filename += fastq->seq_id;
		filename += ".gam";
		std::replace(filename.begin(), filename.end(), '/', '_');
		coutoutput << "write to " << filename << BufferedWriter::Flush;
		std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
		stream::write_buffered(alignmentOut, alignmentvec, 0);
		coutoutput << "write finished" << BufferedWriter::Flush;
	}
	coutoutput << "thread " << threadnum << " finished with " << alignments.size() << " alignments" << BufferedWriter::Flush;
}

void alignReads(std::string graphFile, std::string fastqFile, std::string seedFile, int numThreads, int bandwidth, std::string alignmentFile, std::string auggraphFile)
{

	vg::Graph graph;
	{
		if (is_file_exist(graphFile)){
			std::cout << "load graph from " << graphFile << std::endl;
			std::ifstream graphfile { graphFile, std::ios::in | std::ios::binary };
			std::vector<vg::Graph> parts;
			std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
				parts.push_back(g);
			};
			stream::for_each(graphfile, lambda);
			graph = mergeGraphs(parts);
			std::cout << "graph is " << GraphSizeInBp(graph) << " bp large" << std::endl;
		}
		else{
			std::cerr << "No graph file exists" << std::endl;
			std::exit(0);
		}
	}
	std::vector<FastQ> fastqs;
	if (is_file_exist(fastqFile)){
		fastqs = loadFastqFromFile(fastqFile);
		std::cout << fastqs.size() << " reads" << std::endl;
	}
	else{
		std::cerr << "No fastq file exists" << std::endl;
		std::exit(0);
	}



	std::map<std::string, std::vector<vg::Alignment>> seeds;
	{
		if (is_file_exist(seedFile)){
			std::ifstream seedfile { seedFile, std::ios::in | std::ios::binary };
			std::function<void(vg::Alignment&)> alignmentLambda = [&seeds](vg::Alignment& a) {
				seeds[a.name()].push_back(a);
			};
			stream::for_each(seedfile, alignmentLambda);
		}
		else{
			std::cerr << "No seeds file exists" << std::endl;
			std::exit(0);
		}
	}
	

	std::vector<std::vector<const FastQ*>> readsPerThread;
	std::map<const FastQ*, std::vector<vg::Alignment>> seedHits;
	std::vector<std::vector<vg::Alignment>> resultsPerThread;
	readsPerThread.resize(numThreads);
	resultsPerThread.resize(numThreads);
	int currentThread = 0;
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		if (seeds.count(fastqs[i].seq_id) == 0) continue;
		seedHits[&(fastqs[i])].insert(seedHits[&(fastqs[i])].end(), seeds[fastqs[i].seq_id].begin(), seeds[fastqs[i].seq_id].end());
		readsPerThread[currentThread].push_back(&(fastqs[i]));
		currentThread++;
		currentThread %= numThreads;
	}

	std::vector<std::thread> threads;

	for (int i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&graph, &readsPerThread, &seedHits, &resultsPerThread, i, bandwidth]() { runComponentMappings(graph, readsPerThread[i], seedHits, resultsPerThread[i], i, bandwidth); });
	}

	for (int i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}

	std::vector<vg::Alignment> alignments;

	for (int i = 0; i < numThreads; i++)
	{
		alignments.insert(alignments.end(), resultsPerThread[i].begin(), resultsPerThread[i].end());
	}

	std::cerr << "final result has " << alignments.size() << " alignments" << std::endl;

	std::ofstream alignmentOut { alignmentFile, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, alignments, 0);
	vg::Graph augmentedGraphAllReads;

	augmentedGraphAllReads = augmentGraphwithAlignment(graph, alignments);
	outputGraph(auggraphFile, augmentedGraphAllReads);
}
