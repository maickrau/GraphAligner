#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include <google/protobuf/util/json_util.h>
#include "Aligner.h"
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerWrapper.h"
#include "MummerSeeder.h"
#include "MinimizerSeeder.h"

struct Seeder
{
	enum Mode
	{
		None, File, Mum, Mem, Minimizer
	};
	Mode mode;
	size_t mumCount;
	size_t memCount;
	size_t mxmLength;
	size_t minimizerLength;
	size_t minimizerWindowSize;
	size_t minimizerCount;
	const MummerSeeder* mummerSeeder;
	const MinimizerSeeder* minimizerSeeder;
	const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds;
	Seeder(const AlignerParams& params, const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds, const MummerSeeder* mummerSeeder, const MinimizerSeeder* minimizerSeeder) :
		mumCount(params.mumCount),
		memCount(params.memCount),
		mxmLength(params.mxmLength),
		minimizerLength(params.minimizerLength),
		minimizerWindowSize(params.minimizerWindowSize),
		minimizerCount(params.minimizerCount),
		mummerSeeder(mummerSeeder),
		minimizerSeeder(minimizerSeeder),
		fileSeeds(fileSeeds)
	{
		mode = Mode::None;
		if (fileSeeds != nullptr)
		{
			assert(minimizerSeeder == nullptr);
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerCount == 0);
			mode = Mode::File;
		}
		if (minimizerSeeder != nullptr)
		{
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerCount != 0);
			mode = Mode::Minimizer;
		}
		if (mummerSeeder != nullptr)
		{
			assert(minimizerSeeder == nullptr);
			assert(fileSeeds == nullptr);
			assert(mumCount != 0 || memCount != 0);
			assert(minimizerCount == 0);
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
				return mummerSeeder->getMumSeeds(seq, mumCount, mxmLength);
			case Mode::Mem:
				assert(mummerSeeder != nullptr);
				return mummerSeeder->getMemSeeds(seq, memCount, mxmLength);
			case Mode::Minimizer:
				assert(minimizerSeeder != nullptr);
				return minimizerSeeder->getSeeds(seq, minimizerCount);
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
	seedsFound(0),
	seedsExtended(0),
	readsWithASeed(0),
	alignments(0),
	fullLengthAlignments(0),
	readsWithAnAlignment(0),
	bpInReads(0),
	bpInReadsWithASeed(0),
	bpInAlignments(0),
	bpInFullAlignments(0),
	assertionBroke(false)
	{
	}
	std::atomic<size_t> reads;
	std::atomic<size_t> seeds;
	std::atomic<size_t> seedsFound;
	std::atomic<size_t> seedsExtended;
	std::atomic<size_t> readsWithASeed;
	std::atomic<size_t> alignments;
	std::atomic<size_t> fullLengthAlignments;
	std::atomic<size_t> readsWithAnAlignment;
	std::atomic<size_t> bpInReads;
	std::atomic<size_t> bpInReadsWithASeed;
	std::atomic<size_t> bpInAlignments;
	std::atomic<size_t> bpInFullAlignments;
	std::atomic<bool> assertionBroke;
};

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph)
{
	for (int i = 0; i < alignment.path().mapping_size(); i++)
	{
		int digraphNodeId = alignment.path().mapping(i).position().node_id();
		int originalNodeId = digraphNodeId / 2;
		alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_node_id(originalNodeId);
		std::string name = graph.OriginalNodeName(digraphNodeId);
		if (name.size() > 0)
		{
			alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_name(name);
		}
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

void readFastqs(const std::vector<std::string>& filenames, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& writequeue, std::atomic<bool>& readStreamingFinished)
{
	assertSetRead("Read streamer", "No seed");
	for (auto filename : filenames)
	{
		FastQ::streamFastqFromFile(filename, false, [&writequeue](FastQ& read)
		{
			std::shared_ptr<FastQ> ptr = std::make_shared<FastQ>();
			std::swap(*ptr, read);
			size_t slept = 0;
			while (writequeue.size_approx() > 200)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				slept++;
				if (slept > 100) break;
			}
			writequeue.enqueue(ptr);
		});
	}
	readStreamingFinished = true;
}

void consumeVGsAndWrite(const std::string& filename, moodycamel::ConcurrentQueue<std::string*>& writequeue, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, std::atomic<bool>& allThreadsDone, std::atomic<bool>& allWriteDone, bool verboseMode, bool outputJSON)
{
	assertSetRead("Writer", "No seed");
	auto openmode = std::ios::out;
	if (!outputJSON) openmode |= std::ios::binary;
	std::ofstream outfile { filename, openmode };

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
			if (!writequeue.try_dequeue(alns[0]))
			{
				if (allThreadsDone) break;
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				continue;
			}
			gotAlns = 1;
		}
		coutoutput << "write " << gotAlns << ", " << writequeue.size_approx() << " left" << BufferedWriter::Flush;
		for (size_t i = 0; i < gotAlns; i++)
		{
			outfile.write(alns[i]->data(), alns[i]->size());
		}
		deallocqueue.enqueue_bulk(alns, gotAlns);
		wroteAny = true;
	}

	if (!outputJSON && !wroteAny)
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
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, std::max(params.initialBandwidth, params.rampBandwidth), !params.highMemory };
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
				stats.seedsFound += seeds.size();
				stats.readsWithASeed += 1;
				stats.bpInReadsWithASeed += fastq->sequence.size();
				alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, seeds, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping);
			}
			else
			{
				alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, !params.verboseMode, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping);
			}
		}
		catch (const ThreadReadAssertion::AssertionFailure& a)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
			reusableState.clear();
			stats.assertionBroke = true;
			continue;
		}

		//failed alignment, don't output
		if (alignments.alignments.size() == 0)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			continue;
		}

		stats.seedsExtended += alignments.seedsExtended;
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
		::google::protobuf::io::ZeroCopyOutputStream *raw_out;
		::google::protobuf::io::GzipOutputStream *gzip_out;
		::google::protobuf::io::CodedOutputStream *coded_out;
		if (!params.outputJSON)
		{
			raw_out = new ::google::protobuf::io::OstreamOutputStream(&strstr);
		    gzip_out = new ::google::protobuf::io::GzipOutputStream(raw_out);
		    coded_out = new ::google::protobuf::io::CodedOutputStream(gzip_out);
			coded_out->WriteVarint64(alignments.alignments.size());
		}
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
				stats.assertionBroke = true;
				continue;
			}
			stats.alignments += 1;
			if (alignments.alignments[i].alignment->sequence().size() == fastq->sequence.size())
			{
				stats.fullLengthAlignments += 1;
				stats.bpInFullAlignments += alignments.alignments[i].alignment->sequence().size();
			}
			stats.bpInAlignments += alignments.alignments[i].alignment->sequence().size();
			replaceDigraphNodeIdsWithOriginalNodeIds(*alignments.alignments[i].alignment, alignmentGraph);
			alignmentpositions += std::to_string(alignments.alignments[i].alignmentStart) + "-" + std::to_string(alignments.alignments[i].alignmentEnd) + ", ";
			timems += alignments.alignments[i].elapsedMilliseconds;
			totalcells += alignments.alignments[i].cellsProcessed;
			if (params.outputJSON)
			{
				google::protobuf::util::JsonPrintOptions options;
				options.preserve_proto_field_names = true;
				std::string s;
				google::protobuf::util::MessageToJsonString(*alignments.alignments[i].alignment, &s, options);
				strstr << s;
				strstr << '\n';
			}
			else
			{
				std::string s;
				alignments.alignments[i].alignment->SerializeToString(&s);
				coded_out->WriteVarint32(s.size());
				coded_out->WriteRaw(s.data(), s.size());
			}
		}
		if (!params.outputJSON)
		{
			delete coded_out;
			delete gzip_out;
			delete raw_out;
		}
		std::string* writeAlns = new std::string { strstr.str() };
		size_t waited = 0;
		while (!alignmentsOut.try_enqueue(token, writeAlns) && !alignmentsOut.try_enqueue(writeAlns))
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			waited++;
			if (waited >= 1000)
			{
				if (alignmentsOut.size_approx() < 100 && alignmentsOut.enqueue(writeAlns)) break;
			}
		}
		alignmentpositions.pop_back();
		alignmentpositions.pop_back();

		coutoutput << "Read " << fastq->seq_id << " alignment took " << timems << "ms" << BufferedWriter::Flush;
		coutoutput << "Read " << fastq->seq_id << " aligned by thread " << threadnum << " with positions: " << alignmentpositions << " (read " << fastq->sequence.size() << "bp)" << BufferedWriter::Flush;
	}
	assertSetRead("After all reads", "No seed");
	coutoutput << "Thread " << threadnum << " finished" << BufferedWriter::Flush;
}

AlignmentGraph getGraph(std::string graphFile, MummerSeeder** mxmSeeder, MinimizerSeeder** minimizerSeeder, const AlignerParams& params)
{
	bool loadMxmSeeder = params.mumCount > 0 || params.memCount > 0;
	bool loadMinimizerSeeder = params.minimizerCount > 0;
	bool tryDAG = params.maxCellsPerSlice == std::numeric_limits<size_t>::max();
	if (is_file_exist(graphFile)){
		std::cout << "Load graph from " << graphFile << std::endl;
	}
	else{
		std::cerr << "No graph file exists" << std::endl;
		std::exit(0);
	}
	try
	{
		if (graphFile.substr(graphFile.size()-3) == ".vg")
		{
			if (loadMxmSeeder || loadMinimizerSeeder)
			{
				auto graph = CommonUtils::LoadVGGraph(graphFile);
				if (loadMxmSeeder)
				{
					std::cout << "Build MUM/MEM seeder from the graph" << std::endl;
					*mxmSeeder = new MummerSeeder { graph, params.seederCachePrefix };
				}
				if (loadMinimizerSeeder)
				{
					std::cout << "Build minimizer seeder from the graph" << std::endl;
					*minimizerSeeder = new MinimizerSeeder(graph, params.minimizerLength, params.minimizerWindowSize);
				}
				return DirectedGraph::BuildFromVG(graph, tryDAG);
			}
			else
			{
				return DirectedGraph::StreamVGGraphFromFile(graphFile, tryDAG);
			}
		}
		else if (graphFile.substr(graphFile.size() - 4) == ".gfa")
		{
			auto graph = GfaGraph::LoadFromFile(graphFile, true);
			if (loadMxmSeeder)
			{
				std::cout << "Build MUM/MEM seeder from the graph" << std::endl;
				*mxmSeeder = new MummerSeeder { graph, params.seederCachePrefix };
			}
			if (loadMinimizerSeeder)
			{
				std::cout << "Build minimizer seeder from the graph" << std::endl;
				*minimizerSeeder = new MinimizerSeeder(graph, params.minimizerLength, params.minimizerWindowSize);
			}
			return DirectedGraph::BuildFromGFA(graph, tryDAG);
		}
		else
		{
			std::cerr << "Unknown graph type (" << graphFile << ")" << std::endl;
			std::exit(0);
		}
	}
	catch (const CommonUtils::InvalidGraphException& e)
	{
		std::cout << "Error in the graph: " << e.what() << std::endl;
		std::cerr << "Error in the graph: " << e.what() << std::endl;
		std::exit(1);
	}
}

void alignReads(AlignerParams params)
{
	assertSetRead("Preprocessing", "No seed");

	const std::unordered_map<std::string, std::vector<SeedHit>>* seedHitsToThreads = nullptr;
	std::unordered_map<std::string, std::vector<SeedHit>> seedHits;
	MummerSeeder* mummerseeder = nullptr;
	MinimizerSeeder* minimizerseeder = nullptr;
	auto alignmentGraph = getGraph(params.graphFile, &mummerseeder, &minimizerseeder, params);

	if (params.seedFiles.size() > 0)
	{
		for (auto file : params.seedFiles)
		{
			if (is_file_exist(file)){
				std::cout << "Load seeds from " << file << std::endl;
				std::ifstream seedfile { file, std::ios::in | std::ios::binary };
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

	Seeder seeder { params, seedHitsToThreads, mummerseeder, minimizerseeder };

	switch(seeder.mode)
	{
		case Seeder::Mode::File:
			std::cout << "Seeds from file" << std::endl;
			break;
		case Seeder::Mode::Mum:
			std::cout << "MUM seeds, min length " << seeder.mxmLength << ", max count " << seeder.mumCount << std::endl;
			break;
		case Seeder::Mode::Mem:
			std::cout << "MEM seeds, min length " << seeder.mxmLength << ", max count " << seeder.memCount << std::endl;
			break;
		case Seeder::Mode::Minimizer:
			std::cout << "Minimizer seeds, length " << seeder.minimizerLength << ", window size " << seeder.minimizerWindowSize << ", max count " << seeder.minimizerCount << std::endl;
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

	moodycamel::ConcurrentQueue<std::string*> outputAlns;
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
	std::thread fastqThread { [files=params.fastqFiles, &readFastqsQueue, &readStreamingFinished]() { readFastqs(files, readFastqsQueue, readStreamingFinished); } };
	std::thread writerThread { [file=params.outputAlignmentFile, &outputAlns, &deallocAlns, &allThreadsDone, &allWriteDone, verboseMode=params.verboseMode, outputJSON=params.outputJSON]() { consumeVGsAndWrite(file, outputAlns, deallocAlns, allThreadsDone, allWriteDone, verboseMode, outputJSON); } };
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
	if (minimizerseeder != nullptr) delete minimizerseeder;

	std::string* dealloc;
	while (deallocAlns.try_dequeue(dealloc))
	{
		delete dealloc;
	}

	std::cout << "Alignment finished" << std::endl;
	std::cout << "Input reads: " << stats.reads << " (" << stats.bpInReads << "bp)" << std::endl;
	std::cout << "Seeds found: " << stats.seedsFound << std::endl;
	std::cout << "Seeds extended: " << stats.seedsExtended << std::endl;
	std::cout << "Reads with a seed: " << stats.readsWithASeed << " (" << stats.bpInReadsWithASeed << "bp)" << std::endl;
	std::cout << "Reads with an alignment: " << stats.readsWithAnAlignment << std::endl;
	std::cout << "Output alignments: " << stats.alignments << " (" << stats.bpInAlignments << "bp)" << std::endl;
	std::cout << "Output end-to-end alignments: " << stats.fullLengthAlignments << " (" << stats.bpInFullAlignments << "bp)" << std::endl;
	if (stats.assertionBroke)
	{
		std::cout << "Alignment broke with some reads. Look at stderr output." << std::endl;
	}
}
