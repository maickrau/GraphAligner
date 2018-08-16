#include <queue>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"
#include "fastqloader.h"

struct MummerSeed
{
	size_t readpos;
	size_t len;
	int nodeId;
	size_t nodepos;
	bool reverse;
};

class AlignmentLengthCompare
{
public:
	bool operator()(const MummerSeed& left, const MummerSeed& right) const
	{
		return left.len > right.len;
	}
};

vg::Alignment createAlignment(const std::string& readname, const MummerSeed& seed)
{
	vg::Alignment result;
	result.set_name(readname);
	result.set_query_position(seed.readpos);
	auto path = new vg::Path;
	result.set_allocated_path(path);
	auto vgmapping = path->add_mapping();
	auto position = new vg::Position;
	vgmapping->set_allocated_position(position);
	position->set_node_id(seed.nodeId);
	position->set_is_reverse(seed.reverse);
	position->set_offset(seed.nodepos);
	auto edit = vgmapping->add_edit();
	edit->set_from_length(seed.len);
	edit->set_to_length(seed.len);
	return result;
}

int getNodeIndex(size_t pos, const std::vector<size_t>& nodeMappingPositions)
{
	auto iter = std::upper_bound(nodeMappingPositions.begin(), nodeMappingPositions.end(), pos);
	int index = iter - nodeMappingPositions.begin();
	assert(index > 0);
	assert(index <= nodeMappingPositions.size());
	return index-1;
}

int main(int argc, char** argv)
{
	std::string outputFileName { argv[1] };
	std::string gfaReferenceFilename { argv[2] };
	int maxSeeds = std::stoi(argv[3]);
	std::string readFile { argv[4] };
	std::unordered_map<std::string, size_t> readLengths;
	std::unordered_map<int, size_t> nodeLengths;

	{
		auto reads = loadFastqFromFile(readFile);
		for (auto read : reads)
		{
			readLengths[read.seq_id] = read.sequence.size();
		}
	}
	{
		auto reads = loadFastqFromFile(gfaReferenceFilename);
		for (size_t i = 0; i < reads.size(); i++)
		{
			nodeLengths[std::stoi(reads[i].seq_id)] = reads[i].sequence.size();
		}
	}

	std::unordered_map<std::string, std::priority_queue<MummerSeed, std::vector<MummerSeed>, AlignmentLengthCompare>> alignments;
	size_t numElems = 0;
	std::string currentRead;
	std::string line;
	bool currentReverse = false;
	size_t currentReadLength;
	std::priority_queue<MummerSeed, std::vector<MummerSeed>, AlignmentLengthCompare>* currentQueue;
	while (std::getline(std::cin, line))
	{
		if (line[0] == '>')
		{
			if (line.size() > 8 && (std::string{line.end()-8, line.end()} == " Reverse" || std::string{line.end()-8, line.end()} == "_Reverse"))
			{
				currentReverse = true;
				currentRead = std::string { line.begin()+2, line.end()-8 };
			}
			else
			{
				currentReverse = false;
				currentRead = std::string { line.begin()+2, line.end() };
			}
			currentReadLength = readLengths[currentRead];
			currentQueue = &alignments[currentRead];
		}
		else
		{
			std::stringstream str { line };
			MummerSeed newSeed;
			str >> newSeed.nodeId >> newSeed.nodepos >> newSeed.readpos >> newSeed.len;
			newSeed.reverse = currentReverse;
			assert(newSeed.nodepos >= 1);
			assert(newSeed.readpos >= 1);
			newSeed.nodepos -= 1;
			newSeed.readpos -= 1;
			if (currentReverse)
			{
				newSeed.readpos = currentReadLength - (newSeed.readpos + newSeed.len);
			}
			//there's some weird bug, possibly even in mummer
			//ignore it until we figure out what's going on
			if (newSeed.readpos >= currentReadLength) continue;
			if (newSeed.nodepos >= nodeLengths[newSeed.nodeId]) continue;
			assert(newSeed.readpos < currentReadLength);
			assert(newSeed.nodepos < nodeLengths[newSeed.nodeId]);
			assert(newSeed.readpos >= 0);
			assert(newSeed.nodepos >= 0);
			if (currentQueue->size() < maxSeeds)
			{
				currentQueue->emplace(newSeed);
				numElems++;
			}
			else if (AlignmentLengthCompare{}(newSeed, currentQueue->top()))
			{
				currentQueue->pop();
				currentQueue->emplace(newSeed);
			}
		}
	}
	std::vector<vg::Alignment> writeAlignments;
	writeAlignments.reserve(numElems);
	std::vector<vg::Alignment> insertAlns;
	insertAlns.reserve(maxSeeds);
	for (auto& pair : alignments)
	{
		insertAlns.clear();
		while (pair.second.size() > 0)
		{
			auto aln = createAlignment(pair.first, pair.second.top());
			insertAlns.push_back(aln);
			pair.second.pop();
		}
		std::reverse(insertAlns.begin(), insertAlns.end());
		writeAlignments.insert(writeAlignments.end(), insertAlns.begin(), insertAlns.end());
	}
	assert(writeAlignments.size() == numElems);
	std::ofstream alignmentOut { outputFileName, std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, writeAlignments, 0);
}