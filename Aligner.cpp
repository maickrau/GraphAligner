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
#include "MummerSeeder.h"

struct Seeder
{
	enum Mode
	{
		File, Mum, Mem, None
	};
	Mode mode;
	size_t mumCount;
	size_t memCount;
	size_t mxmLength;
	const MummerSeeder* mummerSeeder;
	const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds;
	Seeder(const AlignerParams& params, const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds, const MummerSeeder* mummerSeeder) :
		mumCount(params.mumCount),
		memCount(params.memCount),
		mxmLength(params.mxmLength),
		mummerSeeder(mummerSeeder),
		fileSeeds(fileSeeds)
	{
		mode = Mode::None;
		if (fileSeeds != nullptr)
		{
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			mode = Mode::File;
		}
		if (mummerSeeder != nullptr)
		{
			assert(fileSeeds == nullptr);
			assert(mumCount != 0 || memCount != 0);
			if (mumCount != 0)
			{
				mode = Mode::Mum;
				assert(memCount == 0);
			}
			if (memCount != 0)
			{
				mode = Mode::Mem;
				assert(mumCount == 0);
			}
		}
	}
	std::vector<SeedHit> getSeeds(const std::string& seqName, const std::string& seq) const
	{
		switch(mode)
		{
			case Mode::File:
				assert(fileSeeds != nullptr);
				if (fileSeeds->count(seqName) == 0) return std::vector<SeedHit>{};
				return fileSeeds->at(seqName);
			case Mode::Mum:
				assert(mummerSeeder != nullptr);
				return mummerSeeder->getMumSeeds(seq, mumCount);
			case Mode::Mem:
				assert(mummerSeeder != nullptr);
				return mummerSeeder->getMemSeeds(seq, memCount);
			case Mode::None:
				assert(false);
		}
		return std::vector<SeedHit>{};
	}
};

struct AlignmentStats
{
	AlignmentStats() :
	reads(0),
	seeds(0),
	readsWithASeed(0),
	alignments(0),
	fullLengthAlignments(0),
	readsWithAnAlignment(0),
	bpInReads(0),
	bpInAlignments(0),
	bpInFullAlignments(0)
	{
	}
	std::atomic<size_t> reads;
	std::atomic<size_t> seeds;
	std::atomic<size_t> readsWithASeed;
	std::atomic<size_t> alignments;
	std::atomic<size_t> fullLengthAlignments;
	std::atomic<size_t> readsWithAnAlignment;
	std::atomic<size_t> bpInReads;
	std::atomic<size_t> bpInAlignments;
	std::atomic<size_t> bpInFullAlignments;
};

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
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

void readFastqs(const std::string& filename, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& writequeue, std::atomic<bool>& readStreamingFinished)
{
	assertSetRead("Read streamer", "No seed");
	FastQ::streamFastqFromFile(filename, false, [&writequeue](FastQ& read)
	{
		std::shared_ptr<FastQ> ptr = std::make_shared<FastQ>();
		std::swap(*ptr, read);
		writequeue.enqueue(ptr);
	});
	readStreamingFinished = true;
}

void consumeVGsAndWrite(const std::string& filename, moodycamel::ConcurrentQueue<std::string*>& writequeue, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, std::atomic<bool>& allThreadsDone, std::atomic<bool>& allWriteDone, bool verboseMode)
{
	assertSetRead("Writer", "No seed");
	std::ofstream outfile { filename, std::ios::binary | std::ios::out };

	bool wroteAny = false;

	std::string* alns[100] {};

	BufferedWriter coutoutput;
	if (verboseMode)
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

void runComponentMappings(const AlignmentGraph& alignmentGraph, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& readFastqsQueue, std::atomic<bool>& readStreamingFinished, int threadnum, const Seeder& seeder, AlignerParams params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, moodycamel::ProducerToken& token, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, AlignmentStats& stats)
{
	assertSetRead("Before any read", "No seed");
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, std::max(params.initialBandwidth, params.rampBandwidth), !params.highMemory, params.useSubgraph };
	BufferedWriter cerroutput;
	BufferedWriter coutoutput;
	if (params.verboseMode)
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
		std::shared_ptr<FastQ> fastq = nullptr;
		while (!readFastqsQueue.try_dequeue(fastq))
		{
			bool tryBreaking = readStreamingFinished;
			if (!readFastqsQueue.try_dequeue(fastq) && tryBreaking) break;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		if (fastq == nullptr) break;
		assertSetRead(fastq->seq_id, "No seed");
		coutoutput << "Read " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" << BufferedWriter::Flush;
		stats.reads += 1;
		stats.bpInReads += fastq->sequence.size();

		AlignmentResult alignments;

		try
		{
			if (seeder.mode != Seeder::Mode::None)
			{
				auto timeStart = std::chrono::system_clock::now();
				std::vector<SeedHit> seeds = seeder.getSeeds(fastq->seq_id, fastq->sequence);
				auto timeEnd = std::chrono::system_clock::now();
				size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
				coutoutput << "Read " << fastq->seq_id << " seeding took " << time << "ms" << BufferedWriter::Flush;
				stats.seeds += seeds.size();
				if (seeds.size() == 0)
				{
					coutoutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " has no seed hits" << BufferedWriter::Flush;
					coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
					continue;
				}
				stats.readsWithASeed += 1;
				if (params.useSubgraph)
				{
					alignments = AlignOneWaySubgraph(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, seeds, reusableState, !params.highMemory);
				}
				else
				{
					alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, seeds, reusableState, !params.highMemory);
				}
			}
			else
			{
				alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, !params.verboseMode, reusableState, !params.highMemory);
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

		stats.readsWithAnAlignment += 1;

		if (!params.outputAllAlns)
		{
			alignments.alignments = CommonUtils::SelectAlignments(alignments.alignments, std::numeric_limits<size_t>::max(), [](const AlignmentResult::AlignmentItem& aln) { return aln.alignment.get(); });
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
				assert(alignments.alignments[i].alignment != nullptr);
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				reusableState.clear();
				continue;
			}
			stats.alignments += 1;
			if (alignments.alignments[i].alignment->sequence().size() == fastq->sequence.size())
			{
				stats.fullLengthAlignments += 1;
				stats.bpInFullAlignments += alignments.alignments[i].alignment->sequence().size();
			}
			stats.bpInAlignments += alignments.alignments[i].alignment->sequence().size();
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
	coutoutput << "Thread " << threadnum << " finished" << BufferedWriter::Flush;
}

AlignmentGraph getGraph(std::string graphFile, MummerSeeder** seeder, bool loadSeeder, const std::string& seederCachePrefix, size_t mxmLength)
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
			std::cout << "Build seeder from the graph" << std::endl;
			*seeder = new MummerSeeder { graph, mxmLength, seederCachePrefix };
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
			std::cout << "Build seeder from the graph" << std::endl;
			*seeder = new MummerSeeder { graph, mxmLength, seederCachePrefix };
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

	const std::unordered_map<std::string, std::vector<SeedHit>>* seedHitsToThreads = nullptr;
	std::unordered_map<std::string, std::vector<SeedHit>> seedHits;

	if (params.seedFile != "")
	{
		{
			if (is_file_exist(params.seedFile)){
				std::cout << "Load seeds from " << params.seedFile << std::endl;
				std::ifstream seedfile { params.seedFile, std::ios::in | std::ios::binary };
				size_t numSeeds = 0;
				std::function<void(vg::Alignment&)> alignmentLambda = [&seedHits, &numSeeds](vg::Alignment& seedhit) {
					seedHits[seedhit.name()].emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.path().mapping(0).position().offset(), seedhit.query_position(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).position().is_reverse());
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

	MummerSeeder* mummerseeder = nullptr;

	auto alignmentGraph = getGraph(params.graphFile, &mummerseeder, params.mumCount != 0 || params.memCount != 0, params.seederCachePrefix, params.mxmLength);

	Seeder seeder { params, seedHitsToThreads, mummerseeder };

	switch(seeder.mode)
	{
		case Seeder::Mode::File:
			std::cout << "Seeds from file " << params.seedFile << std::endl;
			break;
		case Seeder::Mode::Mum:
			std::cout << "MUM seeds, min length " << seeder.mxmLength << ", max count " << seeder.mumCount << std::endl;
			break;
		case Seeder::Mode::Mem:
			std::cout << "MEM seeds, min length " << seeder.mxmLength << ", max count " << seeder.memCount << std::endl;
			break;
		case Seeder::Mode::None:
			std::cout << "No seeds, calculate the entire first row. VERY SLOW!" << std::endl;
			break;
	}

	std::cout << "Initial bandwidth " << params.initialBandwidth;
	if (params.rampBandwidth > 0) std::cout << ", ramp bandwidth " << params.rampBandwidth;
	if (params.maxCellsPerSlice != std::numeric_limits<size_t>::max()) std::cout << ", tangle effort " << params.maxCellsPerSlice;
	std::cout << std::endl;

	std::vector<std::thread> threads;

	assertSetRead("Running alignments", "No seed");

	moodycamel::ConcurrentQueue<std::string*> outputAlns { ((params.numThreads * 50 + 31) / 32) * 32 };
	moodycamel::ConcurrentQueue<std::string*> deallocAlns;
	moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>> readFastqsQueue;
	std::atomic<bool> readStreamingFinished { false };
	std::atomic<bool> allThreadsDone { false };
	std::atomic<bool> allWriteDone { false };
	std::vector<moodycamel::ProducerToken> tokens;
	tokens.reserve(params.numThreads);
	for (size_t i = 0; i < params.numThreads; i++)
	{
		tokens.emplace_back(outputAlns);
	}

	std::cout << "Align" << std::endl;
	AlignmentStats stats;
	std::thread fastqThread { [file=params.fastqFile, &readFastqsQueue, &readStreamingFinished]() { readFastqs(file, readFastqsQueue, readStreamingFinished); } };
	std::thread writerThread { [file=params.outputAlignmentFile, &outputAlns, &deallocAlns, &allThreadsDone, &allWriteDone, verboseMode=params.verboseMode]() { consumeVGsAndWrite(file, outputAlns, deallocAlns, allThreadsDone, allWriteDone, verboseMode); } };
	for (size_t i = 0; i < params.numThreads; i++)
	{
		threads.emplace_back([&alignmentGraph, &readFastqsQueue, &readStreamingFinished, i, seeder, params, &outputAlns, &tokens, &deallocAlns, &stats]() { runComponentMappings(alignmentGraph, readFastqsQueue, readStreamingFinished, i, seeder, params, outputAlns, tokens[i], deallocAlns, stats); });
	}

	for (size_t i = 0; i < params.numThreads; i++)
	{
		threads[i].join();
	}
	assertSetRead("Postprocessing", "No seed");

	allThreadsDone = true;

	writerThread.join();
	fastqThread.join();

	if (mummerseeder != nullptr) delete mummerseeder;

	std::string* dealloc;
	while (deallocAlns.try_dequeue(dealloc))
	{
		delete dealloc;
	}

	std::cout << "Alignment finished" << std::endl;
	std::cout << "Input reads: " << stats.reads << " (" << stats.bpInReads << "bp)" << std::endl;
	std::cout << "Seeds: " << stats.seeds << std::endl;
	std::cout << "Reads with a seed: " << stats.readsWithASeed << std::endl;
	std::cout << "Reads with an alignment: " << stats.readsWithAnAlignment << std::endl;
	std::cout << "Output alignments: " << stats.alignments << " (" << stats.bpInAlignments << "bp)" << std::endl;
	std::cout << "Output end-to-end alignments: " << stats.fullLengthAlignments << " (" << stats.bpInFullAlignments << "bp)" << std::endl;
}
