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
		if (tags.count(node) == 1) result.tags[node] = tags.at(node);
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
				if (ids.count(target.id) == 0) continue;
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
		if (other.tags.count(node.first) == 1) tags[node.first] = other.tags.at(node.first);
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
	std::ifstream file {filename};
	return LoadFromStream(file);
}

GfaGraph GfaGraph::LoadFromStream(std::istream& file)
{
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
			int id;
			std::string dummy;
			std::string seq;
			sstr >> dummy;
			assert(dummy == "S");
			sstr >> id;
			sstr >> seq;
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
	std::vector<NodePos> nonexistantEdges;
	for (auto& edge : result.edges)
	{
		if (result.nodes.count(edge.first.id) == 0)
		{
			nonexistantEdges.push_back(edge.first);
			continue;
		}
		for (size_t i = edge.second.size()-1; i < edge.second.size()+1; i--)
		{
			if (result.nodes.count(edge.second[i].id) == 0) edge.second.erase(edge.second.begin()+i);
		}
	}
	for (auto nonexistant : nonexistantEdges)
	{
		assert(result.edges.find(nonexistant) != result.edges.end());
		result.edges.erase(result.edges.find(nonexistant));
	}
	return result;
}

