#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include "GfaGraph.h"
#include "ThreadReadAssertion.h"
#include "CommonUtils.h"

NodePos::NodePos() :
id(0),
end(false)
{
}

NodePos::NodePos(int id, bool end) :
id(id),
end(end)
{
}

bool NodePos::operator==(const NodePos& other) const
{
	return id == other.id && end == other.end;
}

bool NodePos::operator!=(const NodePos& other) const
{
	return !(*this == other);
}

NodePos NodePos::Reverse() const
{
	return NodePos { id, !end };
}

GfaGraph::GfaGraph() :
nodes(),
edges(),
varyingOverlaps(),
edgeOverlap(std::numeric_limits<size_t>::max())
{
}

GfaGraph GfaGraph::GetSubgraph(const std::unordered_set<int>& ids) const
{
	GfaGraph result;
	result.edgeOverlap = edgeOverlap;
	for (auto node : ids)
	{
		if (nodes.count(node) == 0) continue;
		result.nodes[node] = nodes.at(node);
		if (originalNodeName.count(node) == 1) result.originalNodeName[node] = originalNodeName.at(node);
		if (tags.count(node) == 1) result.tags[node] = tags.at(node);
		NodePos end {node, true};
		if (edges.count(end) == 1)
		{
			for (auto target : edges.at(end))
			{
				if (ids.count(target.id) == 0) continue;
				if (varyingOverlaps.count(std::make_pair(end, target)) == 1) result.varyingOverlaps[std::make_pair(end, target)] = varyingOverlaps.at(std::make_pair(end, target));
				result.edges[end].push_back(target);
			}
		}
		NodePos start {node, false};
		if (edges.count(start) == 1)
		{
			for (auto target : edges.at(start))
			{
				if (ids.count(target.id) == 0) continue;
				if (varyingOverlaps.count(std::make_pair(start, target)) == 1) result.varyingOverlaps[std::make_pair(start, target)] = varyingOverlaps.at(std::make_pair(start, target));
				result.edges[start].push_back(target);
			}
		}
	}
	return result;
}

GfaGraph GfaGraph::GetSubgraph(const std::unordered_set<int>& nodeids, const std::unordered_set<std::pair<NodePos, NodePos>>& selectedEdges) const
{
	GfaGraph result;
	result.edgeOverlap = edgeOverlap;
	for (auto node : nodeids)
	{
		if (nodes.count(node) == 0) continue;
		result.nodes[node] = nodes.at(node);
		if (tags.count(node) == 1) result.tags[node] = tags.at(node);
		NodePos end {node, true};
		if (edges.count(end) == 1)
		{
			for (auto target : edges.at(end))
			{
				if (nodeids.count(target.id) == 0) continue;
				if (selectedEdges.count({end, target}) == 1 || selectedEdges.count({target, end}) == 1) result.edges[end].push_back(target);
			}
		}
		NodePos start {node, false};
		if (edges.count(start) == 1)
		{
			for (auto target : edges.at(start))
			{
				if (nodeids.count(target.id) == 0) continue;
				if (selectedEdges.count({start, target}) == 1 || selectedEdges.count({target, start}) == 1) result.edges[start].push_back(target);
			}
		}
	}
	return result;
}

void GfaGraph::SaveToFile(std::string filename) const
{
	std::ofstream file {filename};
	SaveToStream(file);
}

void GfaGraph::SaveToStream(std::ostream& file) const
{
	for (auto node : nodes)
	{
		file << "S\t" << node.first << "\t" << node.second;
		if (tags.count(node.first) == 1) file << "\t" << tags.at(node.first);
		file << std::endl;
	}
	for (auto edge : edges)
	{
		for (auto target : edge.second)
		{
			auto overlap = edgeOverlap;
			if (varyingOverlaps.count(std::make_pair(edge.first, target)) == 1)
			{
				overlap = varyingOverlaps.at(std::make_pair(edge.first, target));
			}
			file << "L\t" << edge.first.id << "\t" << (edge.first.end ? "+" : "-") << "\t" << target.id << "\t" << (target.end ? "+" : "-") << "\t" << overlap << "M" << std::endl;
		}
	}
}

void GfaGraph::AddSubgraph(const GfaGraph& other)
{
	for (auto node : other.nodes)
	{
		assert(nodes.count(node.first) == 0 || nodes.at(node.first) == node.second);
		nodes[node.first] = node.second;
		if (other.tags.count(node.first) == 1) tags[node.first] = other.tags.at(node.first);
	}
	for (auto edge : other.edges)
	{
		for (auto target : edge.second)
		{
			edges[edge.first].push_back(target);
			if (other.varyingOverlaps.count(std::make_pair(edge.first, target)) == 1)
			{
				varyingOverlaps[std::make_pair(edge.first, target)] = other.varyingOverlaps.at(std::make_pair(edge.first, target));
			}
		}
	}
}

GfaGraph GfaGraph::LoadFromFile(std::string filename, bool allowVaryingOverlaps)
{
	std::ifstream file {filename};
	return LoadFromStream(file, allowVaryingOverlaps);
}

int getNameId(std::unordered_map<std::string, int>& assigned, const std::string& name)
{
	auto found = assigned.find(name);
	if (found == assigned.end())
	{
		int result = assigned.size();
		assigned[name] = result;
		return result;
	}
	return found->second;
}

void GfaGraph::numberBackToIntegers()
{
	std::unordered_map<int, std::string> newNodes;
	std::unordered_map<NodePos, std::vector<NodePos>> newEdges;
	std::unordered_map<std::pair<NodePos, NodePos>, size_t> newVaryingOverlaps;
	std::unordered_map<int, std::string> newTags;
	for (auto pair : varyingOverlaps)
	{
		auto key = pair.first;
		key.first.id = std::stoi(originalNodeName[key.first.id]);
		key.second.id = std::stoi(originalNodeName[key.second.id]);
		newVaryingOverlaps[key] = pair.second;
	}
	for (auto pair : nodes)
	{
		assert(originalNodeName.count(pair.first) == 1);
		newNodes[std::stoi(originalNodeName[pair.first])] = pair.second;
	}
	for (auto edge : edges)
	{
		for (auto target : edge.second)
		{
			newEdges[NodePos { std::stoi(originalNodeName[edge.first.id]), edge.first.end }].push_back(NodePos { std::stoi(originalNodeName[target.id]), target.end });
		}
	}
	for (auto tag : tags)
	{
		newTags[std::stoi(originalNodeName[tag.first])] = tag.second;
	}
	varyingOverlaps = std::move(newVaryingOverlaps);
	nodes = std::move(newNodes);
	edges = std::move(newEdges);
	tags = std::move(newTags);
	originalNodeName.clear();
}

GfaGraph GfaGraph::LoadFromStream(std::istream& file, bool allowVaryingOverlaps)
{
	std::unordered_map<std::string, int> nameMapping;
	GfaGraph result;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		if (line.size() == 0) continue;
		if (line[0] != 'S' && line[0] != 'L') continue;
		if (line[0] == 'S')
		{
			std::stringstream sstr {line};
			std::string idstr;
			std::string dummy;
			std::string seq;
			sstr >> dummy;
			assert(dummy == "S");
			sstr >> idstr;
			int id = getNameId(nameMapping, idstr);
			sstr >> seq;
			assert(seq.size() >= 1);
			std::string tags;
			while (sstr.good())
			{
				char c = sstr.get();
				if (sstr.good() && c != '\r' && c != '\n' && (c != '\t' || tags.size() > 0))
				{
					tags += c;
				}
			}
			result.nodes[id] = seq;
			if (tags.size() > 0) result.tags[id] = tags;
		}
		if (line[0] == 'L')
		{
			std::stringstream sstr {line};
			std::string fromstr;
			std::string tostr;
			std::string fromstart;
			std::string toend;
			std::string dummy;
			int overlap;
			sstr >> dummy;
			assert(dummy == "L");
			sstr >> fromstr;
			int from = getNameId(nameMapping, fromstr);
			sstr >> fromstart;
			sstr >> tostr;
			int to = getNameId(nameMapping, tostr);
			sstr >> toend;
			assert(fromstart == "+" || fromstart == "-");
			assert(toend == "+" || toend == "-");
			sstr >> overlap;
			char dummyc;
			sstr >> dummyc;
			assert(dummyc == 'M');
			if (overlap < 0) throw CommonUtils::InvalidGraphException { "Edge overlap cannot be negative. Fix the graph" };
			assert(overlap >= 0);
			if (!allowVaryingOverlaps && result.edgeOverlap != std::numeric_limits<size_t>::max() && (size_t)overlap != result.edgeOverlap)
			{
				throw CommonUtils::InvalidGraphException { "Varying edge overlaps are not allowed" };
			}
			result.edgeOverlap = overlap;
			NodePos frompos {from, fromstart == "+"};
			NodePos topos {to, toend == "+"};
			result.edges[frompos].push_back(topos);
			if (allowVaryingOverlaps)
			{
				result.varyingOverlaps[std::make_pair(frompos, topos)] = overlap;
			}
		}
	}
	if (result.edges.size() == 0) result.edgeOverlap = 0;
	bool allIdsIntegers = true;
	for (auto pair : nameMapping)
	{
		assert(result.originalNodeName.count(pair.second) == 0);
		result.originalNodeName[pair.second] = pair.first;
		if (allIdsIntegers)
		{
			char* p;
			strtol(pair.first.c_str(), &p, 10);
			if (*p) {
				allIdsIntegers = false;
			}
		}
	}
	if (allIdsIntegers)
	{
		result.numberBackToIntegers();
	}
	std::vector<NodePos> nonexistantEdges;
	bool hasNonexistant = false;
	for (auto& edge : result.edges)
	{
		if (result.nodes.count(edge.first.id) == 0)
		{
			nonexistantEdges.push_back(edge.first);
			for (auto target : edge.second)
			{
				std::cerr << "WARNING: The graph has an edge between non-existant node(s) " << (result.originalNodeName.count(edge.first.id) == 1 ? result.originalNodeName.at(edge.first.id) : std::to_string(edge.first.id)) << (edge.first.end ? "+" : "-") << " and " << (result.originalNodeName.count(target.id) == 1 ? result.originalNodeName.at(target.id) : std::to_string(target.id)) << (target.end ? "+" : "-") << std::endl;
				hasNonexistant = true;
			}
			continue;
		}
		for (size_t i = edge.second.size()-1; i < edge.second.size()+1; i--)
		{
			if (result.nodes.count(edge.second[i].id) == 0)
			{
				std::cerr << "WARNING: The graph has an edge between non-existant node(s) " << (result.originalNodeName.count(edge.first.id) == 1 ? result.originalNodeName.at(edge.first.id) : std::to_string(edge.first.id)) << (edge.first.end ? "+" : "-") << " and " << (result.originalNodeName.count(edge.second[i].id) == 1 ? result.originalNodeName.at(edge.second[i].id) : std::to_string(edge.second[i].id)) << (edge.second[i].end ? "+" : "-") << std::endl;
				hasNonexistant = true;
				edge.second.erase(edge.second.begin()+i);
			}
		}
	}
	if (hasNonexistant)
	{
		std::cerr << "WARNING: Edges between non-existant nodes have been removed." << std::endl;
		std::cout << "WARNING: The graph has edges between non-existant nodes. Check the stderr output." << std::endl;
	}
	for (auto nonexistant : nonexistantEdges)
	{
		assert(result.edges.find(nonexistant) != result.edges.end());
		result.edges.erase(result.edges.find(nonexistant));
	}
	return result;
}

std::string GfaGraph::OriginalNodeName(int nodeId) const
{
	auto found = originalNodeName.find(nodeId);
	if (found == originalNodeName.end()) return "";
	return found->second;
}

void GfaGraph::confirmDoublesidedEdges()
{
	for (auto node : nodes)
	{
		NodePos source;
		source.id = node.first;
		source.end = true;
		NodePos revSource = source;
		revSource.end = false;
		if (edges.count(source) == 1)
		{
			auto targets = edges[source];
			for (auto target : targets)
			{
				bool found = false;
				auto revTarget = target;
				revTarget.end = !revTarget.end;
				for (auto check : edges[revTarget])
				{
					if (check == revSource)
					{
						found = true;
						break;
					}
				}
				if (!found)
				{
					edges[revTarget].emplace_back(revSource);
				}
			}
		}
		if (edges.count(revSource) == 1)
		{
			auto targets = edges[revSource];
			for (auto target : targets)
			{
				bool found = false;
				auto revTarget = target;
				revTarget.end = !revTarget.end;
				for (auto check : edges[revTarget])
				{
					if (check == source)
					{
						found = true;
						break;
					}
				}
				if (!found)
				{
					edges[revTarget].emplace_back(source);
				}
			}
		}
	}
}