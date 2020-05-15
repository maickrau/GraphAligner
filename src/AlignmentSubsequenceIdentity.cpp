#include <vector>
#include <algorithm>
#include <iostream>
#include <tuple>
#include "CommonUtils.h"
#include "GfaGraph.h"
#include "fastqloader.h"

bool fakeLengths = false;

namespace std 
{
	template <> 
	struct hash<std::pair<size_t, size_t>>
	{
		size_t operator()(const std::pair<size_t, size_t>& x) const
		{
			return hash<size_t>()(x.first) ^ hash<size_t>()(x.second);
		}
	};
}

struct Node
{
	int nodeId;
	bool reverse;
	bool operator==(const Node& other) const
	{
		return nodeId == other.nodeId && reverse == other.reverse;
	}
};

struct Alignment
{
	std::vector<Node> path;
	std::vector<size_t> length;
	std::string name;
};

Alignment convertVGtoAlignment(const vg::Alignment& vgAln)
{
	Alignment result;
	result.name = vgAln.name();
	for (int i = 0; i < vgAln.path().mapping_size(); i++)
	{
		result.path.emplace_back();
		result.path.back().nodeId = vgAln.path().mapping(i).position().node_id();
		result.path.back().reverse = vgAln.path().mapping(i).position().is_reverse();
		result.length.emplace_back(vgAln.path().mapping(i).edit(0).to_length());
	}
	return result;
}

Alignment reverse(const Alignment& old)
{
	Alignment result;
	result.name = old.name;
	for (size_t i = 0; i < old.path.size(); i++)
	{
		result.path.emplace_back();
		result.path.back().nodeId = old.path[i].nodeId;
		result.path.back().reverse = !old.path[i].reverse;
	}
	result.length = old.length;
	std::reverse(result.path.begin(), result.path.end());
	std::reverse(result.length.begin(), result.length.end());
	return result;
}

std::pair<double, size_t> getAlignmentIdentity(const Alignment& read, const Alignment& transcript, const std::unordered_map<std::string, size_t>& readLengths)
{
	std::vector<std::vector<size_t>> matchLen;
	matchLen.resize(read.path.size()+1);
	for (size_t i = 0; i < read.path.size()+1; i++)
	{
		matchLen[i].resize(transcript.path.size()+1, 0);
	}
	size_t maxMatch = 0;
	size_t maxmatchIndex = 0;
	size_t maxmatchMissing = 0;
	for (size_t i = 0; i < read.path.size(); i++)
	{
		for (size_t j = 0; j < transcript.path.size(); j++)
		{
			matchLen[i+1][j+1] = std::max(matchLen[i+1][j], matchLen[i][j+1]);
			if (read.path[i] == transcript.path[j])
			{
				matchLen[i+1][j+1] = std::max(matchLen[i+1][j+1], matchLen[i][j] + std::min(read.length[i], transcript.length[j]));
			}
			else
			{
				matchLen[i+1][j+1] = std::max(matchLen[i+1][j+1], matchLen[i][j]);
			}
			if (matchLen[i+1][j+1] > maxMatch)
			{
				maxMatch = matchLen[i+1][j+1];
				maxmatchIndex = j;
				maxmatchMissing = 0;
				if (read.length[i] < transcript.length[j]) maxmatchMissing = transcript.length[j] - read.length[i];
			}
		}
	}
	size_t threeprimeDistance = maxmatchMissing;
	for (size_t i = maxmatchIndex+1; i < transcript.length.size(); i++)
	{
		threeprimeDistance += transcript.length[i];
	}
	assert(maxMatch >= 0);
	double length;
	if (fakeLengths)
	{
		length = 1;
	}
	else
	{
		assert(readLengths.count(read.name) == 1);
		assert(maxMatch <= readLengths.at(read.name));
		length = (double)maxMatch / (double)readLengths.at(read.name);
	}
	return std::make_pair(length, threeprimeDistance);
}

int main(int argc, char** argv)
{
	std::string transcriptFile { argv[1] };
	std::string readAlignmentFile { argv[2] };
	std::string readFastaFile { argv[3] };

	if (argc >= 5 && argv[4][0] == '1') fakeLengths = true;

	std::unordered_map<std::string, size_t> readLengths;
	{
		auto reads = loadFastqFromFile(readFastaFile);
		for (auto read : reads)
		{
			readLengths[read.seq_id] = read.sequence.size();
		}
	}

	std::vector<Alignment> transcripts;
	std::vector<Alignment> reads;
	{
		auto vgtranscripts = CommonUtils::LoadVGAlignments(transcriptFile);
		for (auto vg : vgtranscripts)
		{
			transcripts.push_back(convertVGtoAlignment(vg));
		}
	}
	{
		auto vgreads = CommonUtils::LoadVGAlignments(readAlignmentFile);
		for (auto vg : vgreads)
		{
			reads.push_back(convertVGtoAlignment(vg));
		}
	}

	std::unordered_map<int, std::vector<size_t>> transcriptsCrossingNode;
	for (size_t i = 0; i < transcripts.size(); i++)
	{
		for (int j = 0; j < transcripts[i].path.size(); j++)
		{
			transcriptsCrossingNode[transcripts[i].path[j].nodeId].push_back(i);
		}
	}

	std::unordered_map<std::pair<size_t, size_t>, std::pair<double, size_t>> readTranscriptBestPair;

	for (size_t readi = 0; readi < reads.size(); readi++)
	{
		auto read = reads[readi];
		std::set<size_t> possibleTranscripts;
		for (size_t i = 0; i < read.path.size(); i++)
		{
			possibleTranscripts.insert(transcriptsCrossingNode[read.path[i].nodeId].begin(), transcriptsCrossingNode[read.path[i].nodeId].end());
		}
		auto reverseread = reverse(read);
		for (auto i : possibleTranscripts)
		{
			auto identityFw = getAlignmentIdentity(read, transcripts[i], readLengths);
			auto identityBw = getAlignmentIdentity(reverseread, transcripts[i], readLengths);
			auto bigger = identityFw;
			if (identityBw.first > identityFw.first) bigger = identityBw;
			if (bigger.first > 0 && (readTranscriptBestPair.count(std::make_pair(readi, i)) == 0 || readTranscriptBestPair[std::make_pair(readi, i)].first < bigger.first))
			{
				readTranscriptBestPair[std::make_pair(readi, i)] = bigger;
			}
		}
	}
	for (auto mapping : readTranscriptBestPair)
	{
		std::cout << reads[mapping.first.first].name << "\t" << transcripts[mapping.first.second].name << "\t" << mapping.second.first << "\t" << mapping.second.second << std::endl;
	}
}
