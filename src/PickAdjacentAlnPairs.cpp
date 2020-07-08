#include <unordered_map>
#include <fstream>
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"
#include "fastqloader.h"

void outputPairs(std::ofstream& alignmentOut, const std::string& readname, const std::unordered_map<std::string, size_t>& readLens, const std::vector<vg::Alignment>& starts, const std::vector<vg::Alignment>& ends, const int maxSplitDist, const int minPartialLen)
{
	std::vector<vg::Alignment> pairs;
	size_t currentPairNum = 0;
	for (auto& start : starts)
	{
		assert(start.query_position() == 0);
		int startEnd = 0;
		for (int i = 0; i < start.path().mapping_size(); i++)
		{
			startEnd += start.path().mapping(i).edit(0).to_length();
		}
		assert(startEnd >= minPartialLen);
		for (auto& end : ends)
		{
			int endStart = end.query_position();
			if (abs(startEnd-endStart) > maxSplitDist) continue;
			vg::Alignment left { start };
			vg::Alignment right { end };
			left.set_name(readname + "_pair" + std::to_string(currentPairNum) + "_1");
			right.set_name(readname + "_pair" + std::to_string(currentPairNum) + "_2");
			pairs.push_back(std::move(left));
			pairs.push_back(std::move(right));
			currentPairNum++;
		}
	}
	if (pairs.size() > 0) stream::write_buffered(alignmentOut, pairs, 0);
}

void pickAndWritePairs(std::string inputFile, std::string outputFile, const std::unordered_map<std::string, size_t>& readLens, const int maxSplitDist, const int minPartialLen)
{
	std::ofstream alignmentOut { outputFile, std::ios::out | std::ios::binary };
	std::string currentRead;
	std::vector<vg::Alignment> starts;
	std::vector<vg::Alignment> ends;

	std::ifstream alignmentIn { inputFile, std::ios::in | std::ios::binary };
	std::function<void(vg::Alignment&)> lambda = [&alignmentOut, &readLens, &currentRead, &starts, &ends, maxSplitDist, minPartialLen](vg::Alignment& aln) {
		if (aln.name() != currentRead)
		{
			outputPairs(alignmentOut, currentRead, readLens, starts, ends, maxSplitDist, minPartialLen);
			starts.clear();
			ends.clear();
			currentRead = aln.name();
		}
		assert(readLens.count(aln.name()) == 1);
		size_t alnlen = 0;
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			alnlen += aln.path().mapping(i).edit(0).to_length();
		}
		if (alnlen < minPartialLen) return;
		if (aln.query_position() == 0)
		{
			starts.push_back(aln);
		}
		if (aln.query_position() + alnlen == readLens.at(aln.name()))
		{
			ends.push_back(aln);
		}
	};
	stream::for_each(alignmentIn, lambda);
	outputPairs(alignmentOut, currentRead, readLens, starts, ends, maxSplitDist, minPartialLen);
}

std::unordered_map<std::string, size_t> getReadLens(std::string filename)
{
	std::unordered_map<std::string, size_t> result;
	FastQ::streamFastqFromFile(filename, false, [&result](const FastQ& read)
	{
		result[read.seq_id] = read.sequence.size();
	});
	return result;
}

int main(int argc, char** argv)
{
	std::string inputAlns { argv[1] };
	int maxSplitDist = std::stoi(argv[2]);
	std::string readFile { argv[3] };
	std::string outputAlns { argv[4] };
	int minPartialLen = std::stoi(argv[5]);

	auto readLens = getReadLens(readFile);
	pickAndWritePairs(inputAlns, outputAlns, readLens, maxSplitDist, minPartialLen);
}