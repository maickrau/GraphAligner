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

void mergePartials(const std::unordered_map<std::string, std::vector<PartialAlignment>>& partials, const std::vector<FastQ>& reads)
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
		for (size_t i = 0; i < p.size()-1; i++)
		{
			assert(p[i+1].start > p[i].start);
			if (p[i+1].start < p[i].end)
			{
				size_t grab = (double)p[i].seq.size() * ((double)(p[i+1].start - p[i].start) / (double)(p[i].end - p[i].start));
				correctedSequence += toupper(p[i].seq.substr(0, grab));
			}
			else
			{
				correctedSequence += toupper(p[i].seq);
				correctedSequence += tolower(read.sequence.substr(p[i].end, p[i+1].start - p[i].end));
			}
		}
		correctedSequence += toupper(p.back().seq);
		if (p.back().end < read.sequence.size())
		{
			correctedSequence += tolower(read.sequence.substr(p.back().end));
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


	mergePartials(partials, reads);
}