#include <unordered_set>
#include <vector>
#include <atomic>
#include <unordered_map>
#include <fstream>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "CommonUtils.h"
#include "AlignmentSelection.h"

std::atomic<bool> readingDone;
std::atomic<bool> splittingDone;
std::vector<vg::Alignment*> cleanup;

size_t allAlnsCount = 0;
size_t selectedAlnCount = 0;
size_t fullLengthAlnCount = 0;
size_t readsWithAnAlnCount = 0;
size_t bpInReads = 0;
size_t bpInSelected = 0;
size_t bpInFull = 0;

std::unordered_map<std::string, size_t> getReadLengths(std::string readFile)
{
	std::unordered_map<std::string, size_t> result;
	auto reads = loadFastqFromFile(readFile);
	for (auto read : reads)
	{
		result[read.seq_id] = read.sequence.size();
	}
	return result;
}

void loadAlignments(std::string filename, moodycamel::ConcurrentQueue<vg::Alignment*>& output)
{
	vg::Alignment* current[100];
	size_t countCurrent = 0;
	std::ifstream alnFile { filename, std::ios::in | std::ios::binary };
	std::function<void(vg::Alignment&)> lambda = [&cleanup, &output, &current, &countCurrent](vg::Alignment& g) {
		vg::Alignment* ptr = new vg::Alignment;
		*ptr = g;
		cleanup.push_back(ptr);
		current[countCurrent] = ptr;
		countCurrent++;
		if (countCurrent == 100)
		{
			output.enqueue_bulk(current, 100);
			countCurrent = 0;
		}
	};
	stream::for_each(alnFile, lambda);
	if (countCurrent > 0)
	{
		output.enqueue_bulk(current, countCurrent);
	}

	readingDone = true;
}

void splitAlignmentsIntoSelectedAndFullLength(const std::unordered_map<std::string, size_t>& readLengths, moodycamel::ConcurrentQueue<vg::Alignment*>& inputAlns, moodycamel::ConcurrentQueue<vg::Alignment*>& outputSelected, moodycamel::ConcurrentQueue<vg::Alignment*>& outputFullLength)
{
	vg::Alignment* alns[100] {};

	std::unordered_map<std::string, std::vector<vg::Alignment*>> alnsPerRead;
	while (true)
	{
		size_t gotAlns = inputAlns.try_dequeue_bulk(alns, 100);
		if (gotAlns == 0)
		{
			if (readingDone) break;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			continue;
		}
		for (size_t i = 0; i < gotAlns; i++)
		{
			alnsPerRead[alns[i]->name()].push_back(alns[i]);
		}
	}

	AlignmentSelection::SelectionOptions options;
	options.method = AlignmentSelection::SelectionMethod::GreedyLength;
	for (auto pair : alnsPerRead)
	{
		auto selected = AlignmentSelection::SelectAlignments(pair.second, options);
		outputSelected.enqueue_bulk(selected.data(), selected.size());
		allAlnsCount += pair.second.size();
		selectedAlnCount += selected.size();
		for (auto ptr : selected)
		{
			bpInSelected += ptr->sequence().size();
		}
		if (selected[0]->sequence().size() >= readLengths.at(pair.first) - 1)
		{
			outputFullLength.enqueue(selected[0]);
			bpInFull += selected[0]->sequence().size();
			fullLengthAlnCount += 1;
		}
	}

	readsWithAnAlnCount = alnsPerRead.size();

	splittingDone = true;
}

void writeAlignments(std::string filename, moodycamel::ConcurrentQueue<vg::Alignment*>& inputAlns)
{
	std::ofstream outfile { filename,  std::ios::out | std::ios::binary };

	std::vector<vg::Alignment*> alns;
	alns.resize(1000, nullptr);

	while (true)
	{
		alns.resize(1000);
		size_t gotAlns = inputAlns.try_dequeue_bulk(alns.data(), 1000);
		if (gotAlns == 0)
		{
			if (splittingDone) break;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			continue;
		}
		alns.resize(gotAlns);
		stream::write_buffered_ptr(outfile, alns, 0);
	}
}

int main(int argc, char** argv)
{
	std::string rawAlnFile { argv[1] };
	std::string readsFile { argv[2] };
	std::string outputSelectedAlnFile { argv[3] };
	std::string outputFullLengthAlnFile { argv[4] };
	std::string outputSummaryFile { argv[5] };

	readingDone = false;
	splittingDone = false;

	auto readLengths = getReadLengths(readsFile);

	moodycamel::ConcurrentQueue<vg::Alignment*> readToSplitting;
	moodycamel::ConcurrentQueue<vg::Alignment*> splitToSelected;
	moodycamel::ConcurrentQueue<vg::Alignment*> splitToFullLength;

	std::thread readThread {[&rawAlnFile, &readToSplitting](){loadAlignments(rawAlnFile, readToSplitting);}};
	std::thread splitter {[&readToSplitting, &splitToSelected, &splitToFullLength, &readLengths](){splitAlignmentsIntoSelectedAndFullLength(readLengths, readToSplitting, splitToSelected, splitToFullLength);}};
	std::thread selectedWriter {[&splitToSelected, &outputSelectedAlnFile](){writeAlignments(outputSelectedAlnFile, splitToSelected);}};
	std::thread fullLengthWriter {[&splitToFullLength, &outputFullLengthAlnFile](){writeAlignments(outputFullLengthAlnFile, splitToFullLength);}};

	readThread.join();
	splitter.join();
	selectedWriter.join();
	fullLengthWriter.join();

	for (auto aln : cleanup)
	{
		delete aln;
	}

	for (auto pair : readLengths)
	{
		bpInReads += pair.second;
	}

	std::ofstream summary {outputSummaryFile};
	summary << readLengths.size() << "\tnumber of reads" << std::endl;
	summary << selectedAlnCount << "\tnumber of selected alignments" << std::endl;
	summary << fullLengthAlnCount << "\tnumber of full length alignments" << std::endl;
	summary << readsWithAnAlnCount << "\treads with an alignment" << std::endl;
	summary << bpInReads << "\tbp in reads" << std::endl;
	summary << bpInSelected << "\tbp in selected alignments" << std::endl;
	summary << bpInFull << "\tbp in full length alignments" << std::endl;
}