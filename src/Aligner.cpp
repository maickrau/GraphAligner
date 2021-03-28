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
#include "ReadCorrection.h"
#include "MinimizerSeeder.h"
#include "AlignmentSelection.h"

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
	double minimizerSeedDensity;
	const MummerSeeder* mummerSeeder;
	const MinimizerSeeder* minimizerSeeder;
	const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds;
	Seeder(const AlignerParams& params, const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds, const MummerSeeder* mummerSeeder, const MinimizerSeeder* minimizerSeeder) :
		mumCount(params.mumCount),
		memCount(params.memCount),
		mxmLength(params.mxmLength),
		minimizerLength(params.minimizerLength),
		minimizerWindowSize(params.minimizerWindowSize),
		minimizerSeedDensity(params.minimizerSeedDensity),
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
			assert(minimizerSeedDensity == 0);
			mode = Mode::File;
		}
		if (minimizerSeeder != nullptr)
		{
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerSeedDensity != 0);
			mode = Mode::Minimizer;
		}
		if (mummerSeeder != nullptr)
		{
			assert(minimizerSeeder == nullptr);
			assert(fileSeeds == nullptr);
			assert(mumCount != 0 || memCount != 0);
			assert(minimizerSeedDensity == 0);
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
				return minimizerSeeder->getSeeds(seq, minimizerSeedDensity);
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
	allAlignmentsCount(0),
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
	std::atomic<size_t> allAlignmentsCount;
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

void readFastqs(const std::vector<std::string>& filenames, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& writequeue, std::atomic<bool>& readStreamingFinished)
{
	assertSetNoRead("Read streamer");
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

void consumeBytesAndWrite(const std::string& filename, moodycamel::ConcurrentQueue<std::string*>& writequeue, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, std::atomic<bool>& allThreadsDone, std::atomic<bool>& allWriteDone, bool verboseMode, bool textMode)
{
	assertSetNoRead("Writer");
	auto openmode = std::ios::out;
	if (!textMode) openmode |= std::ios::binary;
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

	if (!textMode && !wroteAny)
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

void QueueInsertSlowly(moodycamel::ProducerToken& token, moodycamel::ConcurrentQueue<std::string*>& queue, std::string&& str)
{
	std::string* write = new std::string { std::move(str) };
	size_t waited = 0;
	while (!queue.try_enqueue(token, write) && !queue.try_enqueue(token, write))
	{
		if (queue.size_approx() < 100 && queue.enqueue(token, write)) break;
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
		waited++;
		if (waited >= 10)
		{
			if (queue.enqueue(token, write)) break;
		}
	}
}

void writeGAMToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	::google::protobuf::io::ZeroCopyOutputStream *raw_out = new ::google::protobuf::io::OstreamOutputStream(&strstr);
	::google::protobuf::io::GzipOutputStream *gzip_out = new ::google::protobuf::io::GzipOutputStream(raw_out);
	::google::protobuf::io::CodedOutputStream *coded_out = new ::google::protobuf::io::CodedOutputStream(gzip_out);
	coded_out->WriteVarint64(alignments.alignments.size());
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].alignment != nullptr);
		std::string s;
		alignments.alignments[i].alignment->SerializeToString(&s);
		coded_out->WriteVarint32(s.size());
		coded_out->WriteRaw(s.data(), s.size());
	}
	delete coded_out;
	delete gzip_out;
	delete raw_out;
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeJSONToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	google::protobuf::util::JsonPrintOptions options;
	options.preserve_proto_field_names = true;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].alignment != nullptr);
		std::string s;
		google::protobuf::util::MessageToJsonString(*alignments.alignments[i].alignment, &s, options);
		strstr << s;
		strstr << '\n';
	}
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeGAFToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].GAFline.size() > 0);
		strstr << alignments.alignments[i].GAFline;
		strstr << '\n';
	}
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeCorrectedToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, const std::string& readName, const std::string& original, size_t maxOverlap, moodycamel::ConcurrentQueue<std::string*>& correctedOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	zstr::ostream *compressed = nullptr;
	if (params.compressCorrected)
	{
		compressed = new zstr::ostream(strstr);
	}
	std::vector<Correction> corrections;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].corrected.size() > 0);
		corrections.emplace_back();
		corrections.back().startIndex = alignments.alignments[i].alignmentStart;
		corrections.back().endIndex = alignments.alignments[i].alignmentEnd;
		corrections.back().corrected = alignments.alignments[i].corrected;
	}
	std::string corrected = getCorrected(original, corrections, maxOverlap);
	if (compressed != nullptr)
	{
		(*compressed) << ">" << readName << std::endl;
		(*compressed) << corrected << std::endl;
		delete compressed;
	}
	else
	{
		strstr << ">" << readName << std::endl;
		strstr << corrected << std::endl;
	}
	QueueInsertSlowly(token, correctedOut, strstr.str());
}

void writeCorrectedClippedToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& correctedClippedOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	zstr::ostream *compressed = nullptr;
	if (params.compressClipped)
	{
		compressed = new zstr::ostream(strstr);
	}
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].corrected.size() > 0);
		if (compressed != nullptr)
		{
			(*compressed) << ">" << alignments.readName << "_" << i << "_" << alignments.alignments[i].alignmentStart << "_" << alignments.alignments[i].alignmentEnd << std::endl;
			(*compressed) << alignments.alignments[i].corrected << std::endl;
		}
		else
		{
			strstr << ">" << alignments.readName << "_" << i << "_" << alignments.alignments[i].alignmentStart << "_" << alignments.alignments[i].alignmentEnd << std::endl;
			strstr << alignments.alignments[i].corrected << std::endl;
		}
	}
	if (compressed != nullptr)
	{
		delete compressed;
	}
	QueueInsertSlowly(token, correctedClippedOut, strstr.str());
}

void runComponentMappings(const AlignmentGraph& alignmentGraph, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& readFastqsQueue, std::atomic<bool>& readStreamingFinished, int threadnum, const Seeder& seeder, AlignerParams params, moodycamel::ConcurrentQueue<std::string*>& GAMOut, moodycamel::ConcurrentQueue<std::string*>& JSONOut, moodycamel::ConcurrentQueue<std::string*>& GAFOut, moodycamel::ConcurrentQueue<std::string*>& correctedOut, moodycamel::ConcurrentQueue<std::string*>& correctedClippedOut, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, AlignmentStats& stats)
{
	moodycamel::ProducerToken GAMToken { GAMOut };
	moodycamel::ProducerToken JSONToken { JSONOut };
	moodycamel::ProducerToken GAFToken { GAFOut };
	moodycamel::ProducerToken correctedToken { correctedOut };
	moodycamel::ProducerToken clippedToken { correctedClippedOut };
	assertSetNoRead("Before any read");
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, std::max(params.initialBandwidth, params.rampBandwidth), !params.highMemory };
	AlignmentSelection::SelectionOptions selectionOptions;
	selectionOptions.method = params.alignmentSelectionMethod;
	selectionOptions.graphSize = alignmentGraph.SizeInBP();
	selectionOptions.ECutoff = params.selectionECutoff;
	selectionOptions.minAlignmentScore = params.minAlignmentScore;
	selectionOptions.EValueCalc = EValueCalculator { params.preciseClippingIdentityCutoff };
	selectionOptions.AlignmentScoreFractionCutoff = params.multimapScoreFraction;
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
		while (fastq == nullptr && !readFastqsQueue.try_dequeue(fastq))
		{
			bool tryBreaking = readStreamingFinished;
			if (!readFastqsQueue.try_dequeue(fastq) && tryBreaking) break;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		if (fastq == nullptr) break;
		assertSetNoRead(fastq->seq_id);
		coutoutput << "Read " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" << BufferedWriter::Flush;
		selectionOptions.readSize = fastq->sequence.size();
		stats.reads += 1;
		stats.bpInReads += fastq->sequence.size();

		AlignmentResult alignments;

		size_t alntimems = 0;
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
					if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, alignments);
					continue;
				}
				stats.seedsFound += seeds.size();
				stats.readsWithASeed += 1;
				stats.bpInReadsWithASeed += fastq->sequence.size();
				auto clusterTimeStart = std::chrono::system_clock::now();
				if (params.multiseedDP)
				{
					PrepareMultiseeds(alignmentGraph, seeds, fastq->sequence.size());
				}
				else
				{
					OrderSeeds(alignmentGraph, seeds);
				}
				auto clusterTimeEnd = std::chrono::system_clock::now();
				size_t clusterTime = std::chrono::duration_cast<std::chrono::milliseconds>(clusterTimeEnd - clusterTimeStart).count();
				coutoutput << "Read " << fastq->seq_id << " clustering took " << clusterTime << "ms" << BufferedWriter::Flush;
				auto alntimeStart = std::chrono::system_clock::now();
				if (params.multiseedDP)
				{
					alignments = AlignMultiseed(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, seeds, reusableState, !params.highMemory, params.forceGlobal, params.seedClusterMinSize, params.seedExtendDensity, params.nondeterministicOptimizations, params.preciseClippingIdentityCutoff, params.Xdropcutoff, params.multimapScoreFraction);
					AlignmentSelection::AddMappingQualities(alignments.alignments);
				}
				else
				{
					alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, seeds, reusableState, !params.highMemory, params.forceGlobal, params.seedClusterMinSize, params.seedExtendDensity, params.nondeterministicOptimizations, params.preciseClippingIdentityCutoff, params.Xdropcutoff);
				}
				auto alntimeEnd = std::chrono::system_clock::now();
				alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
			}
			else
			{
				auto alntimeStart = std::chrono::system_clock::now();
				alignments = AlignOneWay(alignmentGraph, fastq->seq_id, fastq->sequence, params.initialBandwidth, params.rampBandwidth, !params.verboseMode, reusableState, !params.highMemory, params.forceGlobal, params.nondeterministicOptimizations, params.preciseClippingIdentityCutoff, params.Xdropcutoff, params.DPRestartStride);
				auto alntimeEnd = std::chrono::system_clock::now();
				alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
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

		stats.allAlignmentsCount += alignments.alignments.size();

		coutoutput << "Read " << fastq->seq_id << " alignment took " << alntimems << "ms" << BufferedWriter::Flush;
		if (alignments.alignments.size() > 0) alignments.alignments = AlignmentSelection::SelectAlignments(alignments.alignments, selectionOptions);

		//failed alignment, don't output
		if (alignments.alignments.size() == 0)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			try
			{
				if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, alignments);
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				reusableState.clear();
				stats.assertionBroke = true;
				continue;
			}
			continue;
		}

		stats.seedsExtended += alignments.seedsExtended;
		stats.readsWithAnAlignment += 1;

		size_t totalcells = 0;
		for (size_t i = 0; i < alignments.alignments.size(); i++)
		{
			totalcells += alignments.alignments[i].cellsProcessed;
		}
		
		std::sort(alignments.alignments.begin(), alignments.alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });

		if (params.outputGAMFile != "" || params.outputJSONFile != "")
		{
			for (size_t i = 0; i < alignments.alignments.size(); i++)
			{
				AddAlignment(fastq->seq_id, fastq->sequence, alignments.alignments[i]);
				replaceDigraphNodeIdsWithOriginalNodeIds(*alignments.alignments[i].alignment, alignmentGraph);
			}
		}

		if (params.outputGAFFile != "")
		{
			for (size_t i = 0; i < alignments.alignments.size(); i++)
			{
				AddGAFLine(alignmentGraph, fastq->seq_id, fastq->sequence, alignments.alignments[i], params.cigarMatchMismatchMerge);
			}
		}
		
		std::sort(alignments.alignments.begin(), alignments.alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });

		std::string alignmentpositions;

		for (size_t i = 0; i < alignments.alignments.size(); i++)
		{
			stats.alignments += 1;
			size_t alignmentSize = alignments.alignments[i].alignmentEnd - alignments.alignments[i].alignmentStart;
			if (alignmentSize == fastq->sequence.size())
			{
				stats.fullLengthAlignments += 1;
				stats.bpInFullAlignments += alignmentSize;
			}
			stats.bpInAlignments += alignmentSize;
			if (params.outputCorrectedFile != "" || params.outputCorrectedClippedFile != "") AddCorrected(alignments.alignments[i]);
			alignmentpositions += std::to_string(alignments.alignments[i].alignmentStart) + "-" + std::to_string(alignments.alignments[i].alignmentEnd) + ", ";
		}

		alignmentpositions.pop_back();
		alignmentpositions.pop_back();
		coutoutput << "Read " << fastq->seq_id << " aligned by thread " << threadnum << " with positions: " << alignmentpositions << " (read " << fastq->sequence.size() << "bp)" << BufferedWriter::Flush;

		try
		{
			if (params.outputGAMFile != "") writeGAMToQueue(GAMToken, params, GAMOut, alignments);
			if (params.outputJSONFile != "") writeJSONToQueue(JSONToken, params, JSONOut, alignments);
			if (params.outputGAFFile != "") writeGAFToQueue(GAFToken, params, GAFOut, alignments);
			if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, alignments);
			if (params.outputCorrectedClippedFile != "") writeCorrectedClippedToQueue(clippedToken, params, correctedClippedOut, alignments);
		}
		catch (const ThreadReadAssertion::AssertionFailure& a)
		{
			reusableState.clear();
			stats.assertionBroke = true;
			continue;
		}

	}
	assertSetNoRead("After all reads");
	coutoutput << "Thread " << threadnum << " finished" << BufferedWriter::Flush;
}

AlignmentGraph getGraph(std::string graphFile, MummerSeeder** mxmSeeder, const AlignerParams& params)
{
	bool loadMxmSeeder = params.mumCount > 0 || params.memCount > 0;
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
			if (loadMxmSeeder)
			{
				auto graph = CommonUtils::LoadVGGraph(graphFile);
				if (loadMxmSeeder)
				{
					std::cout << "Build MUM/MEM seeder from the graph" << std::endl;
					*mxmSeeder = new MummerSeeder { graph, params.seederCachePrefix };
				}
				std::cout << "Build alignment graph" << std::endl;
				auto result = DirectedGraph::BuildFromVG(graph);
				return result;
			}
			else
			{
				return DirectedGraph::StreamVGGraphFromFile(graphFile);
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
			std::cout << "Build alignment graph" << std::endl;
			auto result = DirectedGraph::BuildFromGFA(graph);
			return result;
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
	assertSetNoRead("Preprocessing");

	const std::unordered_map<std::string, std::vector<SeedHit>>* seedHitsToThreads = nullptr;
	std::unordered_map<std::string, std::vector<SeedHit>> seedHits;
	MummerSeeder* mummerseeder = nullptr;
	auto alignmentGraph = getGraph(params.graphFile, &mummerseeder, params);
	bool loadMinimizerSeeder = params.minimizerSeedDensity != 0;
	MinimizerSeeder* minimizerseeder = nullptr;
	if (loadMinimizerSeeder)
	{
		std::cout << "Build minimizer seeder from the graph" << std::endl;
		minimizerseeder = new MinimizerSeeder(alignmentGraph, params.minimizerLength, params.minimizerWindowSize, params.numThreads, 1.0 - params.minimizerDiscardMostNumerousFraction);
		if (!minimizerseeder->canSeed())
		{
			std::cout << "Warning: Minimizer seeder has no seed hits. Reads cannot be aligned. Try unchopping the graph with vg or a different seeding mode" << std::endl;
		}
	}

	if (params.seedFiles.size() > 0)
	{
		for (auto file : params.seedFiles)
		{
			if (is_file_exist(file)){
				std::cout << "Load seeds from " << file << std::endl;
				std::ifstream seedfile { file, std::ios::in | std::ios::binary };
				size_t numSeeds = 0;
				std::function<void(vg::Alignment&)> alignmentLambda = [&seedHits, &numSeeds](vg::Alignment& seedhit) {
					seedHits[seedhit.name()].emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.path().mapping(0).position().offset(), seedhit.query_position(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).position().is_reverse());
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
			std::cout << "MUM seeds, min length " << seeder.mxmLength;
			if (seeder.mumCount != std::numeric_limits<size_t>::max()) std::cout << ", max count " << seeder.mumCount;
			std::cout << std::endl;
			break;
		case Seeder::Mode::Mem:
			std::cout << "MEM seeds, min length " << seeder.mxmLength;
			if (seeder.memCount != std::numeric_limits<size_t>::max()) std::cout << ", max count " << seeder.memCount;
			std::cout << std::endl;
			break;
		case Seeder::Mode::Minimizer:
			std::cout << "Minimizer seeds, length " << seeder.minimizerLength << ", window size " << seeder.minimizerWindowSize << ", density " << seeder.minimizerSeedDensity << std::endl;
			break;
		case Seeder::Mode::None:
			std::cout << "No seeds, calculate the entire first row. VERY SLOW!" << std::endl;
			break;
	}
	if (seeder.mode != Seeder::Mode::None) std::cout << "Seed cluster size " << params.seedClusterMinSize << std::endl;
	if (seeder.mode != Seeder::Mode::None && params.seedExtendDensity != -1) std::cout << "Extend up to best " << params.seedExtendDensity << " fraction of seeds" << std::endl;

	std::cout << "Initial bandwidth " << params.initialBandwidth;
	if (params.rampBandwidth > 0) std::cout << ", ramp bandwidth " << params.rampBandwidth;
	if (params.maxCellsPerSlice != std::numeric_limits<size_t>::max()) std::cout << ", tangle effort " << params.maxCellsPerSlice;
	std::cout << std::endl;

	if (params.selectionECutoff != -1) std::cout << "Discard alignments with an E-value > " << params.selectionECutoff << std::endl;
	std::cout << "Clip alignment ends with identity < " << params.preciseClippingIdentityCutoff * 100 << "%" << std::endl;
	std::cout << "X-drop DP score cutoff " << params.Xdropcutoff << std::endl;

	if (params.outputGAMFile != "") std::cout << "write alignments to " << params.outputGAMFile << std::endl;
	if (params.outputJSONFile != "") std::cout << "write alignments to " << params.outputJSONFile << std::endl;
	if (params.outputGAFFile != "") std::cout << "write alignments to " << params.outputGAFFile << std::endl;
	if (params.outputCorrectedFile != "") std::cout << "write corrected reads to " << params.outputCorrectedFile << std::endl;
	if (params.outputCorrectedClippedFile != "") std::cout << "write corrected & clipped reads to " << params.outputCorrectedClippedFile << std::endl;

	std::vector<std::thread> threads;

	assertSetNoRead("Running alignments");

	moodycamel::ConcurrentQueue<std::string*> outputGAM { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputGAF { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputJSON { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> deallocAlns;
	moodycamel::ConcurrentQueue<std::string*> outputCorrected { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputCorrectedClipped { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>> readFastqsQueue;
	std::atomic<bool> readStreamingFinished { false };
	std::atomic<bool> allThreadsDone { false };
	std::atomic<bool> GAMWriteDone { false };
	std::atomic<bool> GAFWriteDone { false };
	std::atomic<bool> JSONWriteDone { false };
	std::atomic<bool> correctedWriteDone { false };
	std::atomic<bool> correctedClippedWriteDone { false };

	std::cout << "Align" << std::endl;
	AlignmentStats stats;
	std::thread fastqThread { [files=params.fastqFiles, &readFastqsQueue, &readStreamingFinished]() { readFastqs(files, readFastqsQueue, readStreamingFinished); } };
	std::thread GAMwriterThread { [file=params.outputGAMFile, &outputGAM, &deallocAlns, &allThreadsDone, &GAMWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputGAM, deallocAlns, allThreadsDone, GAMWriteDone, verboseMode, false); else GAMWriteDone = true; } };
	std::thread GAFwriterThread { [file=params.outputGAFFile, &outputGAF, &deallocAlns, &allThreadsDone, &GAFWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputGAF, deallocAlns, allThreadsDone, GAFWriteDone, verboseMode, true); else GAFWriteDone = true; } };
	std::thread JSONwriterThread { [file=params.outputJSONFile, &outputJSON, &deallocAlns, &allThreadsDone, &JSONWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputJSON, deallocAlns, allThreadsDone, JSONWriteDone, verboseMode, true); else JSONWriteDone = true; } };
	std::thread correctedWriterThread { [file=params.outputCorrectedFile, &outputCorrected, &deallocAlns, &allThreadsDone, &correctedWriteDone, verboseMode=params.verboseMode, uncompressed=!params.compressCorrected]() { if (file != "") consumeBytesAndWrite(file, outputCorrected, deallocAlns, allThreadsDone, correctedWriteDone, verboseMode, uncompressed); else correctedWriteDone = true; } };
	std::thread correctedClippedWriterThread { [file=params.outputCorrectedClippedFile, &outputCorrectedClipped, &deallocAlns, &allThreadsDone, &correctedClippedWriteDone, verboseMode=params.verboseMode, uncompressed=!params.compressClipped]() { if (file != "") consumeBytesAndWrite(file, outputCorrectedClipped, deallocAlns, allThreadsDone, correctedClippedWriteDone, verboseMode, uncompressed); else correctedClippedWriteDone = true; } };

	for (size_t i = 0; i < params.numThreads; i++)
	{
		threads.emplace_back([&alignmentGraph, &readFastqsQueue, &readStreamingFinished, i, seeder, params, &outputGAM, &outputJSON, &outputGAF, &outputCorrected, &outputCorrectedClipped, &deallocAlns, &stats]() { runComponentMappings(alignmentGraph, readFastqsQueue, readStreamingFinished, i, seeder, params, outputGAM, outputJSON, outputGAF, outputCorrected, outputCorrectedClipped, deallocAlns, stats); });
	}

	for (size_t i = 0; i < params.numThreads; i++)
	{
		threads[i].join();
	}
	assertSetNoRead("Postprocessing");

	allThreadsDone = true;

	GAMwriterThread.join();
	GAFwriterThread.join();
	JSONwriterThread.join();
	correctedWriterThread.join();
	correctedClippedWriterThread.join();
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
	std::cout << "Alignments: " << stats.alignments << " (" << stats.bpInAlignments << "bp)";
	if (stats.allAlignmentsCount > stats.alignments) std::cout << " (" << (stats.allAlignmentsCount - stats.alignments) << " additional alignments discarded)";
	std::cout << std::endl;
	std::cout << "End-to-end alignments: " << stats.fullLengthAlignments << " (" << stats.bpInFullAlignments << "bp)" << std::endl;
	if (stats.assertionBroke)
	{
		std::cout << "Alignment broke with some reads. Look at stderr output." << std::endl;
	}
}
