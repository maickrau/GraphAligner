#include <algorithm>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <iostream>
#include <functional>
#include "GfaGraph.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"
#include "fastqloader.h"

std::string toupper(std::string seq)
{
	for (auto& c : seq)
	{
		c = toupper(c);
	}
	return seq;
}

std::string tolower(std::string seq)
{
	for (auto& c : seq)
	{
		c = tolower(c);
	}
	return seq;
}

struct PartialAlignment
{
	size_t start;
	size_t end;
	std::string seq;
};

void addPartial(const std::unordered_map<int, int>& ids, std::unordered_map<std::string, std::vector<PartialAlignment>>& partials, std::function<std::string(int)> seqGetter, const vg::Alignment& v)
{
	PartialAlignment result;
	result.start = v.query_position();
	result.end = v.query_position() + v.sequence().size();
	result.seq = "";
	for (int i = 0; i < v.path().mapping_size(); i++)
	{
		auto nodeid = v.path().mapping(i).position().node_id();
		auto sequence = seqGetter(ids.at(nodeid));
		int len = 0;
		for (int j = 0; j < v.path().mapping(i).edit_size(); j++)
		{
			len += v.path().mapping(i).edit(j).from_length();
		}
		if (v.path().mapping(i).position().is_reverse())
		{
			sequence = CommonUtils::ReverseComplement(sequence);
		}
		if (v.path().mapping(i).position().offset() > 0)
		{
			sequence = sequence.substr(v.path().mapping(i).position().offset());
		}
		sequence = sequence.substr(0, len);
		result.seq += sequence;
	}
	partials[v.name()].push_back(result);
}

void addPartial(const vg::Graph& g, const std::unordered_map<int, int>& ids, const vg::Alignment& v, std::unordered_map<std::string, std::vector<PartialAlignment>>& partials)
{
	addPartial(ids, partials, [&g](int id) {return g.node(id).sequence();}, v);
}

void addPartial(const GfaGraph& g, const std::unordered_map<int, int>& ids, const vg::Alignment& v, std::unordered_map<std::string, std::vector<PartialAlignment>>& partials)
{
	addPartial(ids, partials, [&g](int id) {return g.nodes.at(id);}, v);
}

size_t getLongestOverlap(const std::string& left, const std::string& right, size_t maxOverlap)
{
	if (left.size() < maxOverlap) maxOverlap = left.size();
	if (right.size() < maxOverlap) maxOverlap = right.size();
	for (size_t i = maxOverlap; i > 0; i--)
	{
		bool match = true;
		for (size_t a = 0; a < i && match; a++)
		{
			if (left[left.size() - maxOverlap + a] != right[a]) match = false;
		}
		if (match) return i;
	}
	return 0;
}

void mergePartials(const std::unordered_map<std::string, std::vector<PartialAlignment>>& partials, const std::vector<FastQ>& reads, size_t maxOverlap)
{
	for (auto read : reads)
	{
		if (partials.count(read.seq_id) == 0)
		{
			std::cout << ">" << read.seq_id << std::endl << tolower(read.sequence) << std::endl;
			continue;
		}
		auto p = partials.at(read.seq_id);
		std::sort(p.begin(), p.end(), [](const PartialAlignment& left, const PartialAlignment& right) { return left.start < right.start; });
		std::string correctedSequence;
		if (p[0].start > 0)
		{
			correctedSequence = tolower(read.sequence.substr(0, p[0].start));
		}
		assert(p.size() > 0);
		for (size_t i = 0; i < p.size(); i++)
		{
			assert(i == 0 || p[i].start > p[i-1].start);
			if (i > 0 && p[i].start < p[i-1].end)
			{
				size_t overlap = getLongestOverlap(correctedSequence, p[i].seq, maxOverlap);
				correctedSequence += toupper(p[i].seq.substr(overlap));
			}
			else
			{
				if (i > 0 && p[i].start > p[i-1].end) correctedSequence += tolower(read.sequence.substr(p[i-1].end, p[i].start - p[i-1].end));
				correctedSequence += toupper(p[i].seq);
			}
		}
		std::cout << ">" << read.seq_id << std::endl << correctedSequence << std::endl;
	}
}

int main(int argc, char** argv)
{
	std::string graphfilename {argv[1]};
	std::string alnfilename { argv[2] };
	std::string readfilename { argv[3] };
	//output in stdout
	auto reads = loadFastqFromFile(readfilename);
	for (int i = 4; i < argc; i++)
	{
		auto extrareads = loadFastqFromFile(argv[i]);
		reads.insert(reads.end(), extrareads.begin(), extrareads.end());
	}

	size_t maxOverlap = 0;

	std::unordered_map<std::string, std::vector<PartialAlignment>> partials;
	if (graphfilename.substr(graphfilename.size()-3) == ".vg")
	{
		vg::Graph graph = CommonUtils::LoadVGGraph(argv[1]);
		std::unordered_map<int, int> ids;
		for (int i = 0; i < graph.node_size(); i++)
		{
			ids[graph.node(i).id()] = i;
		}
		{
			std::ifstream alnfile { argv[2], std::ios::in | std::ios::binary };
			std::function<void(vg::Alignment&)> lambda = [&graph, &ids, &partials](vg::Alignment& aln) {
				addPartial(graph, ids, aln, partials);
			};
			stream::for_each(alnfile, lambda);
		}
	}
	else if (graphfilename.substr(graphfilename.size() - 4) == ".gfa")
	{
		GfaGraph graph = GfaGraph::LoadFromFile(argv[1]);
		maxOverlap = graph.edgeOverlap;
		std::unordered_map<int, int> ids;
		for (auto node : graph.nodes)
		{
			ids[node.first] = node.first;
		}
		{
			std::ifstream alnfile { argv[2], std::ios::in | std::ios::binary };
			std::function<void(vg::Alignment&)> lambda = [&graph, &ids, &partials](vg::Alignment& aln) {
				addPartial(graph, ids, aln, partials);
			};
			stream::for_each(alnfile, lambda);
		}
	}


	mergePartials(partials, reads, maxOverlap);
}