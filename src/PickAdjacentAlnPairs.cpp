#include <unordered_map>
#include <fstream>
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"
#include "fastqloader.h"

std::vector<vg::Alignment> pickPairs(const std::vector<vg::Alignment>& alns, const std::unordered_map<std::string, size_t>& readLens, int maxSplitDist, int minPartialLen)
{
	std::unordered_map<std::string, std::vector<const vg::Alignment*>> startsPerRead;
	std::unordered_map<std::string, std::vector<const vg::Alignment*>> endsPerRead;
	for (auto& aln : alns)
	{
		assert(readLens.count(aln.name()) == 1);
		if (aln.sequence().size() < minPartialLen) continue;
		if (aln.query_position() == 0)
		{
			startsPerRead[aln.name()].push_back(&aln);
		}
		if (aln.query_position() + aln.sequence().size() == readLens.at(aln.name()))
		{
			endsPerRead[aln.name()].push_back(&aln);
		}
	}
	std::vector<vg::Alignment> result;
	for (auto pair : startsPerRead)
	{
		size_t currentPairNum = 0;
		for (auto start : pair.second)
		{
			for (auto end : endsPerRead[pair.first])
			{
				assert(start->query_position() == 0);
				int startEnd = start->sequence().size();
				int endStart = end->query_position();
				if (abs(startEnd-endStart) > maxSplitDist) continue;
				vg::Alignment left { *start };
				vg::Alignment right { *end };
				left.set_name(pair.first + "_pair" + std::to_string(currentPairNum) + "_1");
				right.set_name(pair.first + "_pair" + std::to_string(currentPairNum) + "_2");
				result.push_back(std::move(left));
				result.push_back(std::move(right));
				currentPairNum++;
			}
		}
	}
	return result;
}

std::unordered_map<std::string, size_t> getReadLens(std::string filename)
{
	auto reads = loadFastqFromFile(filename);
	std::unordered_map<std::string, size_t> result;
	for (auto read : reads)
	{
		result[read.seq_id] = read.sequence.size();
	}
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
	auto alns = CommonUtils::LoadVGAlignments(inputAlns);
	auto pairs = pickPairs(alns, readLens, maxSplitDist, minPartialLen);

	std::ofstream alignmentOut { outputAlns, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, pairs, 0);
}