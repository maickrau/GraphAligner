#include <mutex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerWrapper.h"


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

template <typename T>
T loadFromFile(std::string filename)
{
	T result;
	std::ifstream file {filename};
	if (!file.good()) return result;
	boost::archive::text_iarchive ia(file);
	ia >> result;
	return result;
}

template <typename T>
void saveToFile(const T& vec, std::string filename)
{
	std::ofstream file {filename};
	boost::archive::text_oarchive oa(file);
	oa << vec;
}

void outputGraph(std::string filename, const vg::Graph& graph)
{
	std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {graph};
	stream::write_buffered(alignmentOut, writeVector, 0);
}

void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const std::unordered_map<int, int>& idMapper)
{
	for (int i = 0; i < alignment.path().mapping_size(); i++)
	{
		int digraphNodeId = alignment.path().mapping(i).position().node_id();
		assert(idMapper.count(digraphNodeId) > 0);
		int originalNodeId = idMapper.at(digraphNodeId);
		alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_node_id(originalNodeId);
	}
}

void runComponentMappings(const std::unordered_map<int, int>& newIdToOriginalIdMapper, const AlignmentGraph& alignmentGraph, std::vector<const FastQ*>& fastQs, std::mutex& fastqMutex, std::vector<vg::Alignment>& alignments, int threadnum, int dynamicWidth, int dynamicRowStart, const std::map<const FastQ*, std::vector<std::tuple<int, size_t, bool>>>* graphAlignerSeedHits, int startBandwidth)
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
		coutoutput << "read " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" << BufferedWriter::Flush;

		AlignmentResult alignment;

		if (graphAlignerSeedHits == nullptr)
		{
			alignment = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, dynamicWidth, dynamicRowStart);
		}
		else
		{
			if (graphAlignerSeedHits->find(fastq) == graphAlignerSeedHits->end())
			{
				coutoutput << "read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
				cerroutput << "read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
				coutoutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
				cerroutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
				continue;
			}
			alignment = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, dynamicWidth, dynamicRowStart, graphAlignerSeedHits->at(fastq), startBandwidth);
		}

		coutoutput << "read " << fastq->seq_id << " took " << alignment.elapsedMilliseconds << "ms" << BufferedWriter::Flush;

		//failed alignment, don't output
		if (alignment.alignmentFailed)
		{
			coutoutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			continue;
		}
		if (alignment.alignment.score() == std::numeric_limits<decltype(alignment.alignment.score())>::max())
		{
			coutoutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			continue;
		}

		coutoutput << "read " << fastq->seq_id << " score " << alignment.alignment.score() << BufferedWriter::Flush;
		if (alignment.alignment.score() > fastq->sequence.size() * 0.25)
		{
			cerroutput << "read " << fastq->seq_id << " score is poor: " << alignment.alignment.score() << BufferedWriter::Flush;
		}
		coutoutput << "read " << fastq->seq_id << " alignment positions: " << alignment.alignmentStart << "-" << alignment.alignmentEnd << " (read " << fastq->sequence.size() << "bp)" << BufferedWriter::Flush;

		replaceDigraphNodeIdsWithOriginalNodeIds(alignment.alignment, newIdToOriginalIdMapper);

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

std::pair<AlignmentGraph, std::unordered_map<int, int>> getGraphAndIdMapper(std::string graphFile)
{
	if (is_file_exist(graphFile)){
		std::cout << "load graph from " << graphFile << std::endl;
	}
	else{
		std::cerr << "No graph file exists" << std::endl;
		std::exit(0);
	}
	DirectedGraph augmentedGraph;
	if (graphFile.substr(graphFile.size()-3) == ".vg")
	{
		vg::Graph graph = CommonUtils::LoadVGGraph(graphFile);
		augmentedGraph = { graph };
	}
	else if (graphFile.substr(graphFile.size() - 4) == ".gfa")
	{
		GfaGraph graph = GfaGraph::LoadFromFile(graphFile);
		augmentedGraph = { graph };
	}
	else
	{
		std::cerr << "Unknown graph type (" << graphFile << ")" << std::endl;
		std::exit(0);
	}

	AlignmentGraph alignmentGraph;
	alignmentGraph.ReserveNodes(augmentedGraph.nodes.size(), augmentedGraph.totalSequenceLength);
	for (size_t j = 0; j < augmentedGraph.nodes.size(); j++)
	{
		alignmentGraph.AddNode(augmentedGraph.nodes[j].nodeId, augmentedGraph.nodes[j].sequence, !augmentedGraph.nodes[j].rightEnd);
	}
	for (size_t j = 0; j < augmentedGraph.edges.size(); j++)
	{
		alignmentGraph.AddEdgeNodeId(augmentedGraph.nodes[augmentedGraph.edges[j].fromIndex].nodeId, augmentedGraph.nodes[augmentedGraph.edges[j].toIndex].nodeId);
	}
	alignmentGraph.Finalize(64);

	std::unordered_map<int, int> idMapper;
	idMapper.reserve(augmentedGraph.nodes.size());
	for (size_t j = 0; j < augmentedGraph.nodes.size(); j++)
	{
		if (idMapper.count(augmentedGraph.nodes[j].nodeId) > 0 && idMapper[augmentedGraph.nodes[j].nodeId] != augmentedGraph.nodes[j].originalNodeId)
		{
			assert(idMapper.count(augmentedGraph.nodes[j].nodeId) == 0 || idMapper[augmentedGraph.nodes[j].nodeId] == augmentedGraph.nodes[j].originalNodeId);
			idMapper[augmentedGraph.nodes[j].nodeId] = augmentedGraph.nodes[j].originalNodeId;
		}
		assert(idMapper.count(augmentedGraph.nodes[j].nodeId) == 0 || idMapper[augmentedGraph.nodes[j].nodeId] == augmentedGraph.nodes[j].originalNodeId);
		idMapper[augmentedGraph.nodes[j].nodeId] = augmentedGraph.nodes[j].originalNodeId;
	}

	return std::make_pair(alignmentGraph, idMapper);
}

void alignReads(std::string graphFile, std::string fastqFile, int numThreads, int dynamicWidth, std::string alignmentFile, std::string auggraphFile, int dynamicRowStart, std::string seedFile, int startBandwidth)
{
	assertSetRead("Preprocessing");

	std::vector<FastQ> fastqs;
	if (is_file_exist(fastqFile)){
		fastqs = loadFastqFromFile(fastqFile);
		std::cout << fastqs.size() << " reads" << std::endl;
	}
	else{
		std::cerr << "No fastq file exists" << std::endl;
		std::exit(0);
	}

	const std::map<const FastQ*, std::vector<std::tuple<int, size_t, bool>>>* seedHitsToThreads = nullptr;
	std::map<const FastQ*, std::vector<std::tuple<int, size_t, bool>>> seedHits;

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
				seedHits[&(fastqs[i])].emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.query_position(), seedhit.path().mapping(0).position().is_reverse());
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

	auto graphAndMapper = getGraphAndIdMapper(graphFile);

	std::vector<std::thread> threads;

	assertSetRead("Running alignments");
	std::mutex readMutex;

	for (int i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&graphAndMapper, &readPointers, &readMutex, &resultsPerThread, i, dynamicWidth, dynamicRowStart, seedHitsToThreads, startBandwidth]() { runComponentMappings(graphAndMapper.second, graphAndMapper.first, readPointers, readMutex, resultsPerThread[i], i, dynamicWidth, dynamicRowStart, seedHitsToThreads, startBandwidth); });
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
	stream::write_buffered(alignmentOut, alignments, 0);
	if (auggraphFile != "")
	{
		vg::Graph augmentedGraphAllReads;
		vg::Graph graph = CommonUtils::LoadVGGraph(graphFile);
		augmentedGraphAllReads = augmentGraphwithAlignment(graph, alignments);
		outputGraph(auggraphFile, augmentedGraphAllReads);
	}
}
