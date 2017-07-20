#include <iostream>
#include <vector>
#include <algorithm>
#include "vg.pb.h"
#include "stream.hpp"

struct MafEntry {
	std::string readname;
	int startpos;
	int length;
	bool forward;
};

vg::Graph mergeGraphs(const std::vector<vg::Graph>& parts)
{
	vg::Graph newGraph;
	std::vector<const vg::Node*> allNodes;
	std::vector<const vg::Edge*> allEdges;
	for (size_t i = 0; i < parts.size(); i++)
	{
		for (int j = 0; j < parts[i].node_size(); j++)
		{
			allNodes.push_back(&parts[i].node(j));
		}
		for (int j = 0; j < parts[i].edge_size(); j++)
		{
			allEdges.push_back(&parts[i].edge(j));
		}
	}
	for (size_t i = 0; i < allNodes.size(); i++)
	{
		auto node = newGraph.add_node();
		node->set_id(allNodes[i]->id());
		node->set_sequence(allNodes[i]->sequence());
		node->set_name(allNodes[i]->name());
	}
	for (size_t i = 0; i < allEdges.size(); i++)
	{
		auto edge = newGraph.add_edge();
		edge->set_from(allEdges[i]->from());
		edge->set_to(allEdges[i]->to());
		edge->set_from_start(allEdges[i]->from_start());
		edge->set_to_end(allEdges[i]->to_end());
		edge->set_overlap(allEdges[i]->overlap());
	}
	return newGraph;
}

std::vector<vg::Alignment> mafsToAlignments(const std::vector<MafEntry>& mafs, const std::vector<int>& posToNode, const std::map<int, bool>& nodeIsReverse)
{
	std::vector<vg::Alignment> result;
	for (size_t i = 0; i < mafs.size(); i++)
	{
		std::vector<int> nodeIds;
		nodeIds.push_back(posToNode[mafs[i].startpos]);
		for (int j = 1; j < mafs[i].length; j++)
		{
			if (posToNode[mafs[i].startpos+j] != nodeIds.back())
			{
				nodeIds.push_back(posToNode[mafs[i].startpos+j]);
			}
		}
		if (!mafs[i].forward)
		{
			std::reverse(nodeIds.begin(), nodeIds.end());
		}
		vg::Alignment mafResult;
		mafResult.set_name(mafs[i].readname);
		auto path = new vg::Path;
		mafResult.set_allocated_path(path);
		for (size_t j = 0; j < nodeIds.size(); j++)
		{
			auto vgmapping = path->add_mapping();
			auto position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(j);
			position->set_node_id(nodeIds[j]);
			position->set_is_reverse(nodeIsReverse.at(nodeIds[j]) ^ !mafs[i].forward);
		}
		result.push_back(mafResult);
	}
	return result;
}

std::vector<MafEntry> getMafEntries(std::string filename)
{
	std::vector<MafEntry> result;

	std::ifstream mafFile { filename };
	while (mafFile.good())
	{
		std::string line;
		std::getline(mafFile, line);
		if (line.size() == 0 || line[0] != 'a') continue;
		MafEntry maf;
		std::string checks, checkref;
		mafFile >> checks >> checkref;
		assert(checkref == "ref");
		assert(checks == "s");
		mafFile >> maf.startpos >> maf.length;
		std::getline(mafFile, line);
		mafFile >> checks >> maf.readname;
		assert(checks == "s");
		result.push_back(maf);
	}

	return result;
}

int main(int argc, char** argv)
{
	vg::Graph graph;
	{
		std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
		std::vector<vg::Graph> parts;
		std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
			parts.push_back(g);
		};
		stream::for_each(graphfile, lambda);
		graph = mergeGraphs(parts);
	}

	vg::Alignment referenceAlignment;
	{
		std::ifstream referenceFile { argv[2], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda = [&referenceAlignment](vg::Alignment& g) {
			referenceAlignment = g;
		};
		stream::for_each(referenceFile, lambda);
	}

	std::vector<int> posToNode;
	std::map<int, bool> nodeIsReverse;

	for (int i = 0; i < referenceAlignment.path().mapping_size(); i++)
	{
		auto mapping = referenceAlignment.path().mapping(i);
		int currentNodeSize = mapping.edit(0).to_length();
		for (int j = 0; j < currentNodeSize; j++)
		{
			posToNode.push_back(mapping.position().node_id());
		}
		nodeIsReverse[mapping.position().node_id()] = mapping.position().is_reverse();
	}

	auto mafs = getMafEntries(argv[3]);
	auto alignments = mafsToAlignments(mafs, posToNode, nodeIsReverse);

	std::ofstream alignmentOut { argv[4], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, alignments, 0);

}