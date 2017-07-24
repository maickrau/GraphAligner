#include <mutex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "TopologicalSort.h"
#include "GraphAligner.h"
#include "mfvs_graph.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"


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
vg::Graph augmentGraphwithAlignment(const vg::Graph& graph, const std::vector<vg::Alignment>& alignments)
{
	vg::Graph augmentedGraph;
	std::vector<const vg::Node*> allNodes;
	std::vector<const vg::Edge*> allEdges;
	
	for (int j = 0; j < graph.node_size(); j++)
	{
		allNodes.push_back(&graph.node(j));
	}
// 	for (int j = 0; j < graph.edge_size(); j++)
// 	{
// 		allEdges.push_back(&graph.edge(j));
// 	}
	for (size_t i = 0; i < allNodes.size(); i++)
	{
		auto node = augmentedGraph.add_node();
		node->set_id(allNodes[i]->id());
		node->set_sequence(allNodes[i]->sequence());
		node->set_name(allNodes[i]->name());
	}
// 	for (size_t i = 0; i < allEdges.size(); i++)
// 	{
// 		auto edge = augmentedGraph.add_edge();
// 		edge->set_from(allEdges[i]->from());
// 		edge->set_to(allEdges[i]->to());
// 		edge->set_from_start(allEdges[i]->from_start());
// 		edge->set_to_end(allEdges[i]->to_end());
// 		edge->set_overlap(allEdges[i]->overlap());
// 	}
	
	for(int k=0; k < alignments.size(); k++)
	{
		for (int i = 0; i < alignments[k].path().mapping_size()-1; i++)
		{
			auto edge = augmentedGraph.add_edge();
			edge->set_from(alignments[k].path().mapping(i).position().node_id());
			edge->set_to(alignments[k].path().mapping(i+1).position().node_id());
			edge->set_overlap(0);
			edge->set_from_start(0);
			edge->set_to_end(0);
			if (alignments[k].path().mapping(i).position().is_reverse()) {  
				edge->set_from_start(1);
			}
			if (alignments[k].path().mapping(i+1).position().is_reverse()){
				edge->set_to_end(1);
			}
		}
	}
	return augmentedGraph;
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

void runComponentMappings(const DirectedGraph& augmentedGraph, const GraphAligner<size_t, int32_t, uint64_t>& graphAligner, const AlignmentGraph& alignmentGraph, std::vector<const FastQ*>& fastQs, std::mutex& fastqMutex, std::vector<vg::Alignment>& alignments, int threadnum, int dynamicWidth, int dynamicRowStart, const std::map<const FastQ*, std::vector<std::pair<int, size_t>>>* graphAlignerSeedHits, int startBandwidth)
{
	assertSetRead("Before any read");
	BufferedWriter cerroutput {std::cerr};
	BufferedWriter coutoutput {std::cout};
	while (true)
	{
		const FastQ* fastq;
		size_t fastqSize;
		{
			std::lock_guard<std::mutex> lock {fastqMutex};
			if (fastQs.size() == 0) break;
			fastq = fastQs.back();
			fastQs.pop_back();
			fastqSize = fastQs.size();
		}
		assertSetRead(fastq->seq_id);
		coutoutput << "thread " << threadnum << " " << fastqSize << " left\n";
		coutoutput << "read " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" << "\n";

		decltype(graphAligner.AlignOneWay("", "", 0, 0)) alignment;

		if (graphAlignerSeedHits == nullptr)
		{
			alignment = graphAligner.AlignOneWay(fastq->seq_id, fastq->sequence, dynamicWidth, dynamicRowStart);
		}
		else
		{
			auto augmentedGraphSeedHits = augmentedGraph.GetSeedHits(fastq->sequence, graphAlignerSeedHits->at(fastq));
			std::vector<std::remove_reference<decltype(alignmentGraph)>::type::SeedHit> alignerSeedHits;
			for (size_t i = 0; i < augmentedGraphSeedHits.size(); i++)
			{
				alignerSeedHits.emplace_back(augmentedGraphSeedHits[i].seqPos, augmentedGraphSeedHits[i].nodeId, augmentedGraphSeedHits[i].nodePos);
			}
			alignment = graphAligner.AlignOneWay(fastq->seq_id, fastq->sequence, dynamicWidth, dynamicRowStart, alignerSeedHits, startBandwidth);
		}

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
		if (alignment.maxDistanceFromBand > dynamicWidth * 0.66)
		{
			cerroutput << "read " << fastq->seq_id << " max distance from band is high: " << alignment.maxDistanceFromBand << BufferedWriter::Flush;
		}
		if (alignment.alignment.score() < fastq->sequence.size() * 0.7)
		{
			cerroutput << "read " << fastq->seq_id << " score is low: " << alignment.alignment.score() << BufferedWriter::Flush;
		}

		replaceDigraphNodeIdsWithOriginalNodeIds(alignment.alignment, augmentedGraph);

		alignments.push_back(alignment.alignment);
		coutoutput << "thread " << threadnum << " successfully aligned read " << fastq->seq_id << " with " << alignment.cellsProcessed << " cells" << BufferedWriter::Flush;
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
	assertSetRead("After all reads");
	coutoutput << "thread " << threadnum << " finished with " << alignments.size() << " alignments" << BufferedWriter::Flush;
}

void alignReads(std::string graphFile, std::string fastqFile, int numThreads, int dynamicWidth, std::string alignmentFile, std::string auggraphFile, int dynamicRowStart, std::string seedFile, int startBandwidth)
{
	assertSetRead("Preprocessing");

	if (is_file_exist(graphFile)){
		std::cout << "load graph from " << graphFile << std::endl;
	}
	else{
		std::cerr << "No graph file exists" << std::endl;
		std::exit(0);
	}
	vg::Graph graph = CommonUtils::LoadVGGraph(graphFile);

	std::vector<FastQ> fastqs;
	if (is_file_exist(fastqFile)){
		fastqs = loadFastqFromFile(fastqFile);
		std::cout << fastqs.size() << " reads" << std::endl;
	}
	else{
		std::cerr << "No fastq file exists" << std::endl;
		std::exit(0);
	}

	const std::map<const FastQ*, std::vector<std::pair<int, size_t>>>* seedHitsToThreads = nullptr;
	std::map<const FastQ*, std::vector<std::pair<int, size_t>>> seedHits;

	if (seedFile != "")
	{
		std::map<std::string, std::vector<vg::Alignment>> seeds;
		{
			if (is_file_exist(seedFile)){
				std::ifstream seedfile { seedFile, std::ios::in | std::ios::binary };
				std::function<void(vg::Alignment&)> alignmentLambda = [&seeds](vg::Alignment& a) {
					seeds[a.name()].push_back(a);
				};
				stream::for_each(seedfile, alignmentLambda);
			}
			else {
				std::cerr << "No seeds file exists" << std::endl;
				std::exit(0);
			}
		}
		for (size_t i = 0; i < fastqs.size(); i++)
		{
			for (size_t j = 0; j < seeds[fastqs[i].seq_id].size(); j++)
			{
				auto& seedhit = seeds[fastqs[i].seq_id][j];
				seedHits[&(fastqs[i])].emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.query_position());
			}
		}
		seedHitsToThreads = &seedHits;
	}

	std::vector<const FastQ*> readPointers;
	std::vector<std::vector<vg::Alignment>> resultsPerThread;
	resultsPerThread.resize(numThreads);
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		readPointers.push_back(&(fastqs[i]));
	}

	std::vector<std::thread> threads;

	assertSetRead("Running alignments");
	std::mutex readMutex;

	DirectedGraph augmentedGraph {graph};
	std::cout << "augmented graph out of order before sorting: " << numberOfVerticesOutOfOrder(augmentedGraph) << std::endl;
	if (augmentedGraph.nodes.size() == 0)
	{
		std::cerr << "No nodes in the graph" << std::endl;
		std::abort();
	}
	OrderByFeedbackVertexset(augmentedGraph);
	std::cout << "augmented graph out of order after sorting: " << numberOfVerticesOutOfOrder(augmentedGraph) << std::endl;

	AlignmentGraph alignmentGraph;
	for (size_t j = 0; j < augmentedGraph.nodes.size(); j++)
	{
		alignmentGraph.AddNode(augmentedGraph.nodes[j].nodeId, augmentedGraph.nodes[j].sequence, !augmentedGraph.nodes[j].rightEnd);
	}
	for (size_t j = 0; j < augmentedGraph.edges.size(); j++)
	{
		alignmentGraph.AddEdgeNodeId(augmentedGraph.nodes[augmentedGraph.edges[j].fromIndex].nodeId, augmentedGraph.nodes[augmentedGraph.edges[j].toIndex].nodeId);
	}
	alignmentGraph.Finalize();
	GraphAligner<size_t, int32_t, uint64_t> aligner { alignmentGraph };

	for (int i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&augmentedGraph, &aligner, &alignmentGraph, &readPointers, &readMutex, &resultsPerThread, i, dynamicWidth, dynamicRowStart, seedHitsToThreads, startBandwidth]() { runComponentMappings(augmentedGraph, aligner, alignmentGraph, readPointers, readMutex, resultsPerThread[i], i, dynamicWidth, dynamicRowStart, seedHitsToThreads, startBandwidth); });
	}

	for (int i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	assertSetRead("Postprocessing");

	std::vector<vg::Alignment> alignments;

	for (int i = 0; i < numThreads; i++)
	{
		alignments.insert(alignments.end(), resultsPerThread[i].begin(), resultsPerThread[i].end());
	}

	std::cerr << "final result has " << alignments.size() << " alignments" << std::endl;

	std::ofstream alignmentOut { alignmentFile, std::ios::out | std::ios::binary };
	if (auggraphFile != "")
	{
		vg::Graph augmentedGraphAllReads;
		augmentedGraphAllReads = augmentGraphwithAlignment(graph, alignments);
		stream::write_buffered(alignmentOut, alignments, 0);
		outputGraph(auggraphFile, augmentedGraphAllReads);
	}
}
