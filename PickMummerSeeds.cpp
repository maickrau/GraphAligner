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
	std::string refNodesFileName { argv[3] };
	int maxSeeds = std::stoi(argv[4]);
	std::string readFile { argv[5] };
	std::unordered_map<std::string, size_t> readLengths;

	{
		auto reads = loadFastqFromFile(readFile);
		for (auto read : reads)
		{
			readLengths[read.seq_id] = read.sequence.size();
		}
	}

	std::unordered_map<std::string, std::priority_queue<MummerSeed, std::vector<MummerSeed>, AlignmentLengthCompare>> alignments;
	size_t numElems = 0;
	std::vector<size_t> refNodes;
	{
		std::ifstream refNodesFile { refNodesFileName };
		while (refNodesFile.good())
		{
			size_t found;
			refNodesFile >> found;
			if (!refNodesFile.good()) break;
			refNodes.push_back(found);
		}
	}
	std::vector<size_t> nodeMappingPositions;
	{
		std::ifstream mappingfile { gfaReferenceFilename };
		std::string line;
		std::getline(mappingfile, line);
		char lastChar = '-';
		size_t currentPos = 0;
		while (mappingfile.good())
		{
			if (lastChar == 'N' && mappingfile.peek() == 'N')
			{
				nodeMappingPositions.push_back(currentPos);
			}
			lastChar = mappingfile.get();
			currentPos++;
		}
		nodeMappingPositions.push_back(currentPos+1);
	}
	std::string currentRead;
	std::string line;
	bool currentReverse = false;
	while (std::getline(std::cin, line))
	{
		if (line[0] == '>')
		{
			if (line.size() > 8 && std::string{line.end()-8, line.end()} == " Reverse")
			{
				currentReverse = true;
				currentRead = std::string { line.begin()+2, line.end()-8 };
			}
			else
			{
				currentReverse = false;
				currentRead = std::string { line.begin()+2, line.end() };
			}
		}
		else
		{
			std::stringstream str { line };
			MummerSeed newSeed;
			size_t seqpos;
			str >> seqpos >> newSeed.readpos >> newSeed.len;
			newSeed.reverse = currentReverse;
			size_t index = getNodeIndex(seqpos, nodeMappingPositions);
			assert(index < refNodes.size());
			newSeed.nodeId = refNodes[index];
			assert(seqpos >= nodeMappingPositions[index]);
			assert(seqpos < nodeMappingPositions[index+1]);
			newSeed.nodepos = seqpos - nodeMappingPositions[index];
			assert(newSeed.nodepos >= 2);
			newSeed.nodepos -= 2;
			if (currentReverse)
			{
				assert(newSeed.readpos <= readLengths[currentRead] + 2 - newSeed.len);
				newSeed.readpos = readLengths[currentRead] + 2 - newSeed.readpos - newSeed.len;
				size_t nodeLen = nodeMappingPositions[index+1] - nodeMappingPositions[index];
				assert(newSeed.nodepos <= nodeLen - 2 - newSeed.len);
				newSeed.nodepos = nodeLen - 2 - newSeed.nodepos - newSeed.len;
			}
			assert(newSeed.readpos > 0);
			newSeed.readpos -= 1;
			assert(newSeed.readpos >= 0);
			if (alignments[currentRead].size() < maxSeeds)
			{
				alignments[currentRead].emplace(newSeed);
				numElems++;
			}
			else if (AlignmentLengthCompare{}(newSeed, alignments[currentRead].top()))
			{
				alignments[currentRead].pop();
				alignments[currentRead].emplace(newSeed);
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