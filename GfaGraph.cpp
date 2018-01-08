#include <fstream>
#include <sstream>
#include "GfaGraph.h"
#include "ThreadReadAssertion.h"

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
edgeOverlap(-1)
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
		NodePos end {node, true};
		if (edges.count(end) == 1)
		{
			for (auto target : edges.at(end))
			{
				if (ids.count(target.id) == 0) continue;
				result.edges[end].push_back(target);
			}
		}
		NodePos start {node, false};
		if (edges.count(start) == 1)
		{
			for (auto target : edges.at(start))
			{
				result.edges[start].push_back(target);
			}
		}
	}
	return result;
}

void GfaGraph::SaveToFile(std::string filename) const
{
	std::ofstream file {filename};
	for (auto node : nodes)
	{
		file << "S\t" << node.first << "\t" << node.second << std::endl;
	}
	for (auto edge : edges)
	{
		for (auto target : edge.second)
		{
			file << "L\t" << edge.first.id << "\t" << (edge.first.end ? "+" : "-") << "\t" << target.id << "\t" << (target.end ? "+" : "-") << "\t" << edgeOverlap << "M" << std::endl;
		}
	}
}

void GfaGraph::AddSubgraph(const GfaGraph& other)
{
	for (auto node : other.nodes)
	{
		assert(nodes.count(node.first) == 0 || nodes.at(node.first) == node.second);
		nodes[node.first] = node.second;
	}
	for (auto edge : other.edges)
	{
		for (auto target : edge.second)
		{
			edges[edge.first].push_back(target);
		}
	}
}

GfaGraph GfaGraph::LoadFromFile(std::string filename)
{
	GfaGraph result;
	std::ifstream file {filename};
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
			int id;
			std::string dummy;
			std::string seq;
			sstr >> dummy;
			assert(dummy == "S");
			sstr >> id;
			sstr >> seq;
			result.nodes[id] = seq;
		}
		if (line[0] == 'L')
		{
			std::stringstream sstr {line};
			int from;
			int to;
			std::string fromstart;
			std::string toend;
			std::string dummy;
			int overlap;
			sstr >> dummy;
			assert(dummy == "L");
			sstr >> from;
			sstr >> fromstart;
			sstr >> to;
			sstr >> toend;
			sstr >> overlap;
			assert(overlap >= 0);
			assert(result.edgeOverlap == -1 || overlap == result.edgeOverlap);
			result.edgeOverlap = overlap;
			NodePos frompos {from, fromstart == "+"};
			NodePos topos {to, toend == "+"};
			result.edges[frompos].push_back(topos);
		}
	}
	return result;
}

