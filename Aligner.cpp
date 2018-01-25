#include <mutex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include "Aligner.h"
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerWrapper.h"

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
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

void outputGraph(std::string filename, const vg::Graph& graph)
{
	std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {graph};
	stream::write_buffered(alignmentOut, writeVector, 0);
}

void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment)
{
	for (int i = 0; i < alignment.path().mapping_size(); i++)
	{
		int digraphNodeId = alignment.path().mapping(i).position().node_id();
		int originalNodeId = digraphNodeId / 2;
		alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_node_id(originalNodeId);
	}
}

void writeTrace(const std::vector<AlignmentResult::TraceItem>& trace, const std::string& filename)
{
	std::ofstream file { filename };
	for (size_t i = 0; i < trace.size(); i++)
	{
		file << trace[i].nodeID << " " << trace[i].offset << " " << (trace[i].reverse ? 1 : 0) << " " << trace[i].readpos << " " << (int)trace[i].type << " " << trace[i].graphChar << " " << trace[i].readChar << std::endl;
	}
}

void runComponentMappings(const AlignmentGraph& alignmentGraph, std::vector<const FastQ*>& fastQs, std::mutex& fastqMutex, std::vector<vg::Alignment>& results, int threadnum, const std::map<const FastQ*, std::vector<std::tuple<int, size_t, bool>>>* graphAlignerSeedHits, AlignerParams params)
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

		AlignmentResult alignments;

		try
		{
			if (graphAlignerSeedHits == nullptr)
			{
				assert(false);
				// alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth);
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
				alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, graphAlignerSeedHits->at(fastq));
			}
		}
		catch (const ThreadReadAssertion::AssertionFailure& a)
		{
			coutoutput << "read " << fastq->seq_id << "alignment failed (assertion!)" << BufferedWriter::Flush;
			cerroutput << "read " << fastq->seq_id << "alignment failed (assertion!)" << BufferedWriter::Flush;
			continue;
		}

		//failed alignment, don't output
		if (alignments.alignments.size() == 0)
		{
			coutoutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			continue;
		}

		std::string alignmentpositions;
		std::vector<vg::Alignment> alignmentvec;
		size_t timems = 0;
		size_t totalcells = 0;
		for (auto& alignment : alignments.alignments)
		{
			try
			{
				assert(!alignment.alignmentFailed());
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				continue;
			}
			replaceDigraphNodeIdsWithOriginalNodeIds(alignment.alignment);
			results.push_back(alignment.alignment);
			alignmentvec.emplace_back(alignment.alignment);
			alignmentpositions += std::to_string(alignment.alignmentStart) + "-" + std::to_string(alignment.alignmentEnd) + ", ";
			timems += alignment.elapsedMilliseconds;
			totalcells += alignment.cellsProcessed;
		}
		alignmentpositions.pop_back();
		alignmentpositions.pop_back();

		coutoutput << "read " << fastq->seq_id << " took " << timems << "ms" << BufferedWriter::Flush;
		coutoutput << "read " << fastq->seq_id << " alignment positions: " << alignmentpositions << " (read " << fastq->sequence.size() << "bp)" << BufferedWriter::Flush;

		coutoutput << "thread " << threadnum << " aligned read " << fastq->seq_id << " with " << totalcells << " cells" << BufferedWriter::Flush;
		std::string filename;
		filename = "alignment_";
		filename += std::to_string(threadnum);
		filename += "_";
		filename += fastq->seq_id;
		filename += ".gam";
		std::replace(filename.begin(), filename.end(), '/', '_');
		std::replace(filename.begin(), filename.end(), ':', '_');
		coutoutput << "write alignment to " << filename << BufferedWriter::Flush;
		std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
		stream::write_buffered(alignmentOut, alignmentvec, 0);
		coutoutput << "alignment write finished" << BufferedWriter::Flush;
		// std::string tracefilename;
		// tracefilename = "trace_";
		// tracefilename += std::to_string(threadnum);
		// tracefilename += "_";
		// tracefilename += fastq->seq_id;
		// tracefilename += ".trace";
		// std::replace(tracefilename.begin(), tracefilename.end(), '/', '_');
		// std::replace(tracefilename.begin(), tracefilename.end(), ':', '_');
		// coutoutput << "write trace to " << tracefilename << BufferedWriter::Flush;
		// writeTrace(alignment.trace, tracefilename);
		// coutoutput << "trace write finished" << BufferedWriter::Flush;
	}
	assertSetRead("After all reads");
	coutoutput << "thread " << threadnum << " finished with " << results.size() << " alignments" << BufferedWriter::Flush;
}

AlignmentGraph getGraph(std::string graphFile)
{
	if (is_file_exist(graphFile)){
		std::cout << "load graph from " << graphFile << std::endl;
	}
	else{
		std::cerr << "No graph file exists" << std::endl;
		std::exit(0);
	}
	if (graphFile.substr(graphFile.size()-3) == ".vg")
	{
		return DirectedGraph::StreamVGGraphFromFile(graphFile);
	}
	else if (graphFile.substr(graphFile.size() - 4) == ".gfa")
	{
		return DirectedGraph::StreamGFAGraphFromFile(graphFile);
	}
	else
	{
		std::cerr << "Unknown graph type (" << graphFile << ")" << std::endl;
		std::exit(0);
	}
}

void alignReads(AlignerParams params)
{
	assertSetRead("Preprocessing");

	std::vector<FastQ> fastqs;
	if (is_file_exist(params.fastqFile)){
		fastqs = loadFastqFromFile(params.fastqFile);
		std::cout << fastqs.size() << " reads" << std::endl;
	}
	else{
		std::cerr << "No fastq file exists" << std::endl;
		std::exit(0);
	}

	const std::map<const FastQ*, std::vector<std::tuple<int, size_t, bool>>>* seedHitsToThreads = nullptr;
	std::map<const FastQ*, std::vector<std::tuple<int, size_t, bool>>> seedHits;

	if (params.seedFile != "")
	{
		std::map<std::string, std::vector<vg::Alignment>> seeds;
		{
			if (is_file_exist(params.seedFile)){
				std::ifstream seedfile { params.seedFile, std::ios::in | std::ios::binary };
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
	resultsPerThread.resize(params.numThreads);
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		readPointers.push_back(&(fastqs[i]));
	}

	auto alignmentGraph = getGraph(params.graphFile);

	std::vector<std::thread> threads;

	assertSetRead("Running alignments");
	std::mutex readMutex;

	for (int i = 0; i < params.numThreads; i++)
	{
		threads.emplace_back([&alignmentGraph, &readPointers, &readMutex, &resultsPerThread, i, seedHitsToThreads, params]() { runComponentMappings(alignmentGraph, readPointers, readMutex, resultsPerThread[i], i, seedHitsToThreads, params); });
	}

	for (int i = 0; i < params.numThreads; i++)
	{
		threads[i].join();
	}
	assertSetRead("Postprocessing");

	std::vector<vg::Alignment> alignments;

	for (int i = 0; i < params.numThreads; i++)
	{
		alignments.insert(alignments.end(), resultsPerThread[i].begin(), resultsPerThread[i].end());
	}

	std::cerr << "final result has " << alignments.size() << " alignments" << std::endl;

	if (params.alignmentFile != "")
	{
		std::ofstream alignmentOut { params.alignmentFile, std::ios::out | std::ios::binary };
		stream::write_buffered(alignmentOut, alignments, 0);
	}
	if (params.auggraphFile != "")
	{
		vg::Graph augmentedGraphAllReads;
		vg::Graph graph = CommonUtils::LoadVGGraph(params.graphFile);
		augmentedGraphAllReads = augmentGraphwithAlignment(graph, alignments);
		outputGraph(params.auggraphFile, augmentedGraphAllReads);
	}
}
