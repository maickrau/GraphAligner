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
#include "GraphAlignerBitvectorFull.h"
#include "GraphAlignerCellbycellFull.h"

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

void runComponentMappings(const AlignmentGraph& alignmentGraph, std::vector<const FastQ*>& fastQs, std::mutex& fastqMutex, int threadnum, const std::map<const FastQ*, std::vector<std::tuple<int, size_t, bool>>>* graphAlignerSeedHits, AlignerParams params, size_t& numAlignments, std::vector<vg::Alignment>& alignmentsOut, bool hasMergedAlignmentOut)
{
	assertSetRead("Before any read");
	BufferedWriter cerroutput;
	BufferedWriter coutoutput;
	if (!params.quietMode)
	{
		cerroutput = {std::cerr};
		coutoutput = {std::cout};
	}
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
				alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, params.quietMode, params.sloppyOptimizations, graphAlignerSeedHits->at(fastq));
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
			numAlignments += 1;
			replaceDigraphNodeIdsWithOriginalNodeIds(alignment.alignment);
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
		if (hasMergedAlignmentOut)
		{
			alignmentsOut.insert(alignmentsOut.end(), alignmentvec.begin(), alignmentvec.end());
		}
		else
		{
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
		}
	}
	assertSetRead("After all reads");
	coutoutput << "thread " << threadnum << " finished with " << numAlignments << " alignments" << BufferedWriter::Flush;
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
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		readPointers.push_back(&(fastqs[i]));
	}

	auto alignmentGraph = getGraph(params.graphFile);

	std::vector<std::thread> threads;

	assertSetRead("Running alignments");
	std::mutex readMutex;

	std::vector<size_t> numAlnsPerThread;
	numAlnsPerThread.resize(params.numThreads, 0);

	std::vector<std::vector<vg::Alignment>> resultsPerThread;
	bool hasMergedAlignmentOut = params.outputAlignmentFile != "";
	resultsPerThread.resize(params.numThreads);

	for (int i = 0; i < params.numThreads; i++)
	{
		threads.emplace_back([&alignmentGraph, &readPointers, &readMutex, i, seedHitsToThreads, params, &numAlnsPerThread, &resultsPerThread, hasMergedAlignmentOut]() { runComponentMappings(alignmentGraph, readPointers, readMutex, i, seedHitsToThreads, params, numAlnsPerThread[i], resultsPerThread[i], hasMergedAlignmentOut); });
	}

	for (int i = 0; i < params.numThreads; i++)
	{
		threads[i].join();
	}
	assertSetRead("Postprocessing");

	size_t numAlignments = 0;
	for (size_t i = 0; i < params.numThreads; i++)
	{
		numAlignments += numAlnsPerThread[i];
	}

	std::cout << "final result has " << numAlignments << " alignments" << std::endl;
	if (hasMergedAlignmentOut)
	{
		assert(params.outputAlignmentFile != "");
		std::cout << "merge alignments for writing" << std::endl;
		std::vector<vg::Alignment> finalResult;
		finalResult.reserve(numAlignments);
		for (size_t i = 0; i < resultsPerThread.size(); i++)
		{
			finalResult.insert(finalResult.end(), resultsPerThread[i].begin(), resultsPerThread[i].end());
		}
		std::cout << "write to " << params.outputAlignmentFile << std::endl;
		std::ofstream resultFile { params.outputAlignmentFile, std::ios::out | std::ios::binary };
		stream::write_buffered(resultFile, finalResult, 0);
		std::cout << "write finished" << std::endl;
	}
}

void wabiExperiments(AlignerParams params)
{
	std::cout << "load graph" << std::endl;
	auto alignmentGraph = getGraph(params.graphFile);
	std::cout << "load reads" << std::endl;
	auto fastqs = loadFastqFromFile(params.fastqFile);

	std::cout << "preprocess graph" << std::endl;
	auto componentOrder = alignmentGraph.TopologicalOrderOfComponents();
	bool isAcyclic = true;
	bool isTree = true;
	for (size_t i = 0; i < componentOrder.first.size(); i++)
	{
		if (componentOrder.first[i].size() > 1)
		{
			isAcyclic = false;
			break;
		}
		for (auto neighbor : alignmentGraph.inNeighbors[componentOrder.first[i][0]])
		{
			if (neighbor == componentOrder.first[i][0])
			{
				isAcyclic = false;
				break;
			}
		}
	}
	if (isAcyclic)
	{
		for (size_t i = 0; i < alignmentGraph.NodeSize(); i++)
		{
			if (alignmentGraph.inNeighbors[i].size() > 1)
			{
				isTree = false;
				break;
			}
		}
	}
	else
	{
		isTree = false;
	}

	if (isTree)
	{
		assert(isAcyclic);
		std::cout << "The graph is linear / a tree / a forest" << std::endl;
	}
	else if (isAcyclic)
	{
		std::cout << "The graph is a DAG" << std::endl;
	}
	else
	{
		std::cout << "The graph is cyclic" << std::endl;
	}
	size_t numEdges = 0;
	for (const auto& list : alignmentGraph.inNeighbors)
	{
		numEdges += list.size();
	}
	std::cout << "Collapsed nodes: " << alignmentGraph.NodeSize() << std::endl;
	std::cout << "Collapsed edges: " << numEdges << std::endl;
	std::cout << "Nodes: " << alignmentGraph.SizeInBp() << std::endl;
	std::cout << "Edges: " << numEdges + alignmentGraph.SizeInBp() - alignmentGraph.NodeSize() << std::endl;
	std::cout << "BPs: " << alignmentGraph.SizeInBp() << std::endl;

	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params alignerParams { 0, 0, alignmentGraph, 0, false, false };

	GraphAlignerBitvectorFull<size_t, int32_t, uint64_t> bv { alignerParams, componentOrder.first, componentOrder.second };
	GraphAlignerCellbycellFull<size_t, int32_t, uint64_t> cbc { alignerParams, componentOrder.first, componentOrder.second };

	std::vector<int32_t> bvScores;
	std::vector<int32_t> cbcScores;
	bvScores.reserve(fastqs.size());
	cbcScores.reserve(fastqs.size());
	std::cout << "start bitvector alignment" << std::endl;
	auto bvtimeStart = std::chrono::steady_clock::now();
	for (auto fastq : fastqs)
	{
		size_t score;
		if (isAcyclic)
		{
			score = bv.alignAndGetScoreAcyclic(fastq.sequence);
		}
		else
		{
			score = bv.alignAndGetScore(fastq.sequence);
		}
		bvScores.push_back(score);
	}
	auto bvtimeEnd = std::chrono::steady_clock::now();
	size_t bitvectorMicroseconds = std::chrono::duration_cast<std::chrono::microseconds>(bvtimeEnd - bvtimeStart).count();
	std::cout << "bitvector took " << bitvectorMicroseconds << "us" << std::endl;
	std::cout << "start cellbycell alignment" << std::endl;
	auto cbctimeStart = std::chrono::steady_clock::now();
	for (auto fastq : fastqs)
	{
		size_t score;
		if (isAcyclic)
		{
			score = cbc.alignAndGetScoreAcyclic(fastq.sequence);
		}
		else
		{
			score = cbc.alignAndGetScore(fastq.sequence);
		}
		cbcScores.push_back(score);
	}
	auto cbctimeEnd = std::chrono::steady_clock::now();
	size_t cellbycellMicroseconds = std::chrono::duration_cast<std::chrono::microseconds>(cbctimeEnd - cbctimeStart).count();
	std::cout << "cellbycell took " << cellbycellMicroseconds << "us" << std::endl;
	std::cout << "ratio: " << ((double)cellbycellMicroseconds / (double)bitvectorMicroseconds) << std::endl;
	bool scoresOK = true;
	for (size_t i = 0; i < bvScores.size(); i++)
	{
		if (bvScores[i] != cbcScores[i])
		{
			std::cout << "SCORES DON'T MATCH for read " << fastqs[i].seq_id << ": bitvector " << bvScores[i] << " cellbycell " << cbcScores[i] << std::endl;
			scoresOK = false;
		}
	}
	if (scoresOK) std::cout << "scores match" << std::endl;
}