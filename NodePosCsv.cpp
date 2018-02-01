#include <algorithm>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "CommonUtils.h"

int main(int argc, char** argv)
{
	std::string infilename {argv[1]};
	std::string outfilename {argv[2]};

	std::unordered_map<int, std::unordered_map<std::string, std::vector<std::pair<int, int>>>> positions;
	std::unordered_map<int, std::unordered_map<std::string, int>> minRepeatCounts;
	auto alignments = CommonUtils::LoadVGAlignments(infilename);
	std::unordered_set<std::string> alignmentNames;
	for (auto aln : alignments)
	{
		alignmentNames.insert(aln.name());
		int pos = aln.query_position();
		for (size_t i = 0; i < aln.path().mapping_size(); i++)
		{
			auto mapping = aln.path().mapping(i);
			positions[mapping.position().node_id()][aln.name()].emplace_back(pos, pos+mapping.edit(0).to_length());
			pos += mapping.edit(0).to_length();
			minRepeatCounts[mapping.position().node_id()][aln.name()] += 1;
		}
	}
	std::vector<std::string> readnames { alignmentNames.begin(), alignmentNames.end() };
	std::sort(readnames.begin(), readnames.end());
	std::ofstream out {outfilename};
	out << "node,_numreads,_minalnrepeatcount,_traversingreads";
	for (auto read : readnames)
	{
		out << "," << read;
	}
	out << std::endl;
	std::vector<int> nodevec;
	for (auto node : positions)
	{
		nodevec.push_back(node.first);
	}
	std::sort(nodevec.begin(), nodevec.end());
	for (auto node : nodevec)
	{
		out << node;
		out << "," << positions[node].size();
		int minRepeatCount = 0;
		for (auto pair : minRepeatCounts[node])
		{
			minRepeatCount = std::max(minRepeatCount, pair.second);
		}
		out << "," << minRepeatCount;
		out << ",";
		bool first = true;
		for (auto read : positions[node])
		{
			if (read.second.size() > 0)
			{
				if (!first) out << ";";
				out << read.first;
				first = false;
			}
		}
		for (auto read : readnames)
		{
			out << ",";
			if (positions[node].count(read) == 1)
			{
				for (size_t i = 0; i < positions[node][read].size(); i++)
				{
					if (i > 0) out << ";";
					out << positions[node][read][i].first << "-" << positions[node][read][i].second;
				}
			}
		}
		out << std::endl;
	}
}
