#include <mutex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include "Aligner.h"
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerWrapper.h"
#include "STSeeder.h"

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

void consumeVGsAndWrite(const std::string& filename, moodycamel::ConcurrentQueue<std::string*>& writequeue, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, std::atomic<bool>& allThreadsDone, std::atomic<bool>& allWriteDone, bool quietMode)
{
	assertSetRead("Writer", "No seed");
	std::ofstream outfile { filename, std::ios::binary | std::ios::out };

	bool wroteAny = false;

	std::string* alns[100] {};

	BufferedWriter coutoutput;
	if (!quietMode)
	{
		coutoutput = {std::cout};
	}

	while (true)
	{
		size_t gotAlns = writequeue.try_dequeue_bulk(alns, 100);
		if (gotAlns == 0)
		{
			if (allThreadsDone) break;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			continue;
		}
		coutoutput << "write " << gotAlns << ", " << writequeue.size_approx() << " left" << BufferedWriter::Flush;
		for (size_t i = 0; i < gotAlns; i++)
		{
			outfile.write(alns[i]->data(), alns[i]->size());
		}
		deallocqueue.enqueue_bulk(alns, gotAlns);
		wroteAny = true;
	}

	if (!wroteAny)
	{
		::google::protobuf::io::ZeroCopyOutputStream *raw_out =
		      new ::google::protobuf::io::OstreamOutputStream(&outfile);
		::google::protobuf::io::GzipOutputStream *gzip_out =
		      new ::google::protobuf::io::GzipOutputStream(raw_out);
		::google::protobuf::io::CodedOutputStream *coded_out =
		      new ::google::protobuf::io::CodedOutputStream(gzip_out);
		coded_out->WriteVarint64(0);
		delete coded_out;
		delete gzip_out;
		delete raw_out;
	}

	allWriteDone = true;
}

void runComponentMappings(const AlignmentGraph& alignmentGraph, std::vector<const FastQ*>& fastQs, std::mutex& fastqMutex, int threadnum, const std::unordered_map<std::string, std::vector<SeedHit>>* graphAlignerSeedHits, const STSeeder* seeder, AlignerParams params, size_t& numAlignments, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, moodycamel::ProducerToken& token, moodycamel::ConcurrentQueue<std::string*>& deallocqueue)
{
	assertSetRead("Before any read", "No seed");
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, std::max(params.initialBandwidth, params.rampBandwidth), params.lowMemory, params.useSubgraph };
	BufferedWriter cerroutput;
	BufferedWriter coutoutput;
	if (!params.quietMode)
	{
		cerroutput = {std::cerr};
		coutoutput = {std::cout};
	}
	while (true)
	{
		std::string* dealloc;
		while (deallocqueue.try_dequeue(dealloc))
		{
			delete dealloc;
		}
		const FastQ* fastq;
		size_t fastqSize;
		{
			std::lock_guard<std::mutex> lock {fastqMutex};
			if (fastQs.size() == 0) break;
			fastq = fastQs.back();
			fastQs.pop_back();
			fastqSize = fastQs.size();
		}
		assertSetRead(fastq->seq_id, "No seed");
		coutoutput << "Thread " << threadnum << " " << fastqSize << " left\n";
		coutoutput << "Read " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" << BufferedWriter::Flush;

		AlignmentResult alignments;

		try
		{
			if (seeder != nullptr || graphAlignerSeedHits != nullptr)
			{
				std::vector<SeedHit> seeds;
				if (seeder != nullptr)
				{
					auto seedTimeStart = std::chrono::system_clock::now();
					seeds = seeder->getMumSeeds(fastq->sequence);
					auto seedTimeEnd = std::chrono::system_clock::now();
					auto seedTimems = std::chrono::duration_cast<std::chrono::milliseconds>(seedTimeEnd - seedTimeStart).count();
					coutoutput << "Read " << fastq->seq_id << " seeding took " << seedTimems << "ms" << BufferedWriter::Flush;
				}
				else
				{
					assert(graphAlignerSeedHits != nullptr);
					auto found = graphAlignerSeedHits->find(fastq->seq_id);
					if (found != graphAlignerSeedHits->end()) seeds = found->second;
				}
				if (seeds.size() == 0)
				{
					coutoutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					continue;
				}
				if (params.useSubgraph)
				{
					alignments = AlignOneWaySubgraph(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, params.quietMode, params.sloppyOptimizations, seeder->getMumSeeds(fastq->sequence), reusableState, params.lowMemory);
				}
				else
				{
					alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, params.quietMode, params.sloppyOptimizations, seeder->getMumSeeds(fastq->sequence), reusableState, params.lowMemory);
				}
			}
			else
			{
				alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.quietMode, reusableState, params.lowMemory);
			}
		}
		catch (const ThreadReadAssertion::AssertionFailure& a)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
			reusableState.clear();
			continue;
		}

		//failed alignment, don't output
		if (alignments.alignments.size() == 0)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			continue;
		}

		if (params.maxAlns != 0)
		{
			alignments.alignments = CommonUtils::SelectAlignments(alignments.alignments, params.maxAlns, [](const AlignmentResult::AlignmentItem& aln) { return aln.alignment.get(); });
		}
		
		std::sort(alignments.alignments.begin(), alignments.alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });

		std::string alignmentpositions;
		size_t timems = 0;
		size_t totalcells = 0;
		std::stringstream strstr;
		::google::protobuf::io::ZeroCopyOutputStream *raw_out =
		      new ::google::protobuf::io::OstreamOutputStream(&strstr);
		::google::protobuf::io::GzipOutputStream *gzip_out =
		      new ::google::protobuf::io::GzipOutputStream(raw_out);
		::google::protobuf::io::CodedOutputStream *coded_out =
		      new ::google::protobuf::io::CodedOutputStream(gzip_out);
		coded_out->WriteVarint64(alignments.alignments.size());
		for (size_t i = 0; i < alignments.alignments.size(); i++)
		{
			try
			{
				assert(!alignments.alignments[i].alignmentFailed());
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				reusableState.clear();
				continue;
			}
			numAlignments += 1;
			replaceDigraphNodeIdsWithOriginalNodeIds(*alignments.alignments[i].alignment);
			alignmentpositions += std::to_string(alignments.alignments[i].alignmentStart) + "-" + std::to_string(alignments.alignments[i].alignmentEnd) + ", ";
			timems += alignments.alignments[i].elapsedMilliseconds;
			totalcells += alignments.alignments[i].cellsProcessed;
			std::string s;
			alignments.alignments[i].alignment->SerializeToString(&s);
			coded_out->WriteVarint32(s.size());
			coded_out->WriteRaw(s.data(), s.size());
		}
		delete coded_out;
		delete gzip_out;
		delete raw_out;
		std::string* writeAlns = new std::string { strstr.str() };
		while (!alignmentsOut.try_enqueue(token, writeAlns))
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		alignmentpositions.pop_back();
		alignmentpositions.pop_back();

		coutoutput << "Read " << fastq->seq_id << " alignment took " << timems << "ms" << BufferedWriter::Flush;
		coutoutput << "Read " << fastq->seq_id << " aligned by thread " << threadnum << " with positions: " << alignmentpositions << " (read " << fastq->sequence.size() << "bp)" << BufferedWriter::Flush;
	}
	assertSetRead("After all reads", "No seed");
	coutoutput << "Thread " << threadnum << " finished with " << numAlignments << " alignments" << BufferedWriter::Flush;
}

AlignmentGraph getGraph(std::string graphFile, STSeeder** seeder, bool loadSeeder)
{
	if (is_file_exist(graphFile)){
		std::cout << "Load graph from " << graphFile << std::endl;
	}
	else{
		std::cerr << "No graph file exists" << std::endl;
		std::exit(0);
	}
	if (graphFile.substr(graphFile.size()-3) == ".vg")
	{
		if (loadSeeder)
		{
			auto graph = CommonUtils::LoadVGGraph(graphFile);
			*seeder = new STSeeder { graph };
			return DirectedGraph::BuildFromVG(graph);
		}
		else
		{
			return DirectedGraph::StreamVGGraphFromFile(graphFile);
		}
	}
	else if (graphFile.substr(graphFile.size() - 4) == ".gfa")
	{
		if (loadSeeder)
		{
			auto graph = GfaGraph::LoadFromFile(graphFile);
			*seeder = new STSeeder { graph };
			return DirectedGraph::BuildFromGFA(graph);
		}
		else
		{
			return DirectedGraph::StreamGFAGraphFromFile(graphFile);
		}
	}
	else
	{
		std::cerr << "Unknown graph type (" << graphFile << ")" << std::endl;
		std::exit(0);
	}
}

void alignReads(AlignerParams params)
{
	assertSetRead("Preprocessing", "No seed");

	std::vector<FastQ> fastqs;
	if (is_file_exist(params.fastqFile)){
		std::cout << "Load reads from " << params.fastqFile << std::endl;
		fastqs = loadFastqFromFile(params.fastqFile, false);
		std::cout << fastqs.size() << " reads" << std::endl;
	}
	else{
		std::cerr << "No fastq file exists" << std::endl;
		std::exit(0);
	}

	const std::unordered_map<std::string, std::vector<SeedHit>>* seedHitsToThreads = nullptr;
	std::unordered_map<std::string, std::vector<SeedHit>> seedHits;

	if (params.seedFile != "")
	{
		for (auto read : fastqs)
		{
			seedHits[read.seq_id] = {};
		}
		{
			if (is_file_exist(params.seedFile)){
				std::cout << "Load seeds from " << params.seedFile << std::endl;
				std::ifstream seedfile { params.seedFile, std::ios::in | std::ios::binary };
				size_t numSeeds = 0;
				std::function<void(vg::Alignment&)> alignmentLambda = [&seedHits, &numSeeds](vg::Alignment& seedhit) {
					auto pos = seedHits.find(seedhit.name());
					if (pos == seedHits.end()) return;
					pos->second.emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.path().mapping(0).position().offset(), seedhit.query_position(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).position().is_reverse());
					numSeeds += 1;
				};
				stream::for_each(seedfile, alignmentLambda);
				std::cout << numSeeds << " seeds" << std::endl;
			}
			else {
				std::cerr << "No seeds file exists" << std::endl;
				std::exit(0);
			}
		}
		seedHitsToThreads = &seedHits;
	}

	std::vector<const FastQ*> readPointers;
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		readPointers.push_back(&(fastqs[i]));
	}

	STSeeder* seeder = nullptr;

	auto alignmentGraph = getGraph(params.graphFile, &seeder, params.mums);

	std::vector<std::thread> threads;

	assertSetRead("Running alignments", "No seed");
	std::mutex readMutex;

	std::vector<size_t> numAlnsPerThread;
	numAlnsPerThread.resize(params.numThreads, 0);

	moodycamel::ConcurrentQueue<std::string*> outputAlns { ((params.numThreads * 50 + 31) / 32) * 32 };
	moodycamel::ConcurrentQueue<std::string*> deallocAlns;
	std::atomic<bool> allThreadsDone { false };
	std::atomic<bool> allWriteDone { false };
	std::vector<moodycamel::ProducerToken> tokens;
	tokens.reserve(params.numThreads);
	for (int i = 0; i < params.numThreads; i++)
	{
		tokens.emplace_back(outputAlns);
	}

	std::cout << "Align" << std::endl;
	std::thread writerThread { [file=params.outputAlignmentFile, &outputAlns, &deallocAlns, &allThreadsDone, &allWriteDone, quietMode=params.quietMode]() { consumeVGsAndWrite(file, outputAlns, deallocAlns, allThreadsDone, allWriteDone, quietMode); } };
	for (int i = 0; i < params.numThreads; i++)
	{
		threads.emplace_back([&alignmentGraph, &readPointers, &readMutex, i, seedHitsToThreads, params, &numAlnsPerThread, &outputAlns, &tokens, &deallocAlns, seeder]() { runComponentMappings(alignmentGraph, readPointers, readMutex, i, seedHitsToThreads, seeder, params, numAlnsPerThread[i], outputAlns, tokens[i], deallocAlns); });
	}

	for (int i = 0; i < params.numThreads; i++)
	{
		threads[i].join();
	}
	assertSetRead("Postprocessing", "No seed");

	allThreadsDone = true;

	writerThread.join();

	if (seeder != nullptr) delete seeder;

	std::string* dealloc;
	while (deallocAlns.try_dequeue(dealloc))
	{
		delete dealloc;
	}

	size_t numAlignments = 0;
	for (size_t i = 0; i < params.numThreads; i++)
	{
		numAlignments += numAlnsPerThread[i];
	}

	std::cout << "Final result has " << numAlignments << " alignments" << std::endl;
}
