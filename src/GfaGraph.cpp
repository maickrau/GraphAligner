#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <zstr.hpp> //https://github.com/mateidavid/zstr
#include "GfaGraph.h"
#include "ThreadReadAssertion.h"
#include "CommonUtils.h"

bool operator<(const NodePos& left, const NodePos& right)
{
	if (left.id < right.id) return true;
	if (right.id < left.id) return false;
	if (!left.end && right.end) return true;
	if (left.end && !right.end) return false;
	assert(left == right);
	return false;
}

bool operator>(const NodePos& left, const NodePos& right)
{
	return right < left;
}

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
edges()
{
}

GfaGraph GfaGraph::LoadFromFile(std::string filename)
{
	if (filename.substr(filename.size()-3) == ".gz")
	{
		assert(filename.substr(filename.size()-7) == ".gfa.gz");
		zstr::ifstream file { filename };
		return LoadFromStream(file);
	}
	else
	{
		assert(filename.substr(filename.size()-4) == ".gfa");
		std::ifstream file { filename };
		return LoadFromStream(file);
	}
	zstr::ifstream file { filename };
}

size_t getNameId(std::unordered_map<std::string, size_t>& assigned, const std::string& name, std::vector<DNAString>& nodeSeqs, std::vector<std::string>& originalNodeName)
{
	auto found = assigned.find(name);
	if (found == assigned.end())
	{
		assert(assigned.size() == nodeSeqs.size());
		assert(assigned.size() == originalNodeName.size());
		int result = assigned.size();
		assigned[name] = result;
		originalNodeName.emplace_back(name);
		nodeSeqs.emplace_back();
		return result;
	}
	assert(found->second < originalNodeName.size());
	assert(name == originalNodeName[found->second]);
	return found->second;
}

GfaGraph GfaGraph::LoadFromStream(std::istream& file)
{
	std::unordered_map<std::string, size_t> nameMapping;
	GfaGraph result;
	while (file.good())
	{
		std::string line = "";
		std::getline(file, line);
		if (line.size() == 0 && !file.good()) break;
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
			size_t id = getNameId(nameMapping, idstr, result.nodes, result.originalNodeName);
			sstr >> seq;
			if (seq == "*") throw CommonUtils::InvalidGraphException { std::string { "Nodes without sequence (*) are not currently supported (nodeid " + idstr + ")" } };
			assert(seq.size() >= 1);
			result.nodes[id] = seq;
		}
		if (line[0] == 'L')
		{
			std::stringstream sstr {line};
			std::string fromstr;
			std::string tostr;
			std::string fromstart;
			std::string toend;
			std::string dummy;
			std::string overlapstr;
			int overlap = 0;
			sstr >> dummy;
			assert(dummy == "L");
			sstr >> fromstr;
			size_t from = getNameId(nameMapping, fromstr, result.nodes, result.originalNodeName);
			sstr >> fromstart;
			sstr >> tostr;
			size_t to = getNameId(nameMapping, tostr, result.nodes, result.originalNodeName);
			sstr >> toend;
			assert(fromstart == "+" || fromstart == "-");
			assert(toend == "+" || toend == "-");
			sstr >> overlapstr;
			if (overlapstr == "*")
			{
				throw CommonUtils::InvalidGraphException { "Unspecified edge overlaps (*) are not supported" };
			}
			if (overlapstr == "")
			{
				throw CommonUtils::InvalidGraphException { "Edge overlap missing between edges " + fromstr + " and " + tostr };
			}
			assert(overlapstr.size() >= 1);
			size_t charAfterIndex = 0;
			overlap = std::stol(overlapstr, &charAfterIndex, 10);
			if (charAfterIndex != overlapstr.size() - 1 || overlapstr.back() != 'M')
			{
				throw CommonUtils::InvalidGraphException { "Edge overlaps other than exact match are not supported (non supported overlap: " + overlapstr + ")" };
			}
			if (overlap < 0) throw CommonUtils::InvalidGraphException { std::string { "Edge overlap between nodes " + std::to_string(from) + " and " + std::to_string(to) + " is negative" } };
			assert(overlap >= 0);
			NodePos frompos {(int)from, fromstart == "+"};
			NodePos topos {(int)to, toend == "+"};
			result.edges.emplace_back(frompos, topos, overlap);
			result.edges.emplace_back(topos.Reverse(), frompos.Reverse(), overlap);
		}
	}
	for (size_t i = 0; i < result.nodes.size(); i++)
	{
		if (result.nodes[i].size() > 0) continue;
		throw CommonUtils::InvalidGraphException { std::string { "Node " + result.originalNodeName[i] + " is present in edges but missing in nodes" } };
	}
	for (auto t : result.edges)
	{
		auto from = std::get<0>(t);
		auto to = std::get<1>(t);
		auto overlap = std::get<2>(t);
		assert(from.id < result.nodes.size());
		assert(to.id < result.nodes.size());
		if (result.nodes.at(from.id).size() <= overlap || result.nodes.at(to.id).size() <= overlap)
		{
			throw CommonUtils::InvalidGraphException { std::string { "Overlap between nodes " + result.originalNodeName.at(from.id) + " and " + result.originalNodeName.at(to.id) + " fully contains one of the nodes. Fix the overlap to be strictly smaller than both nodes" } };
		}
	}
	for (const auto& t : result.edges)
	{
		if (std::get<0>(t).id >= result.nodes.size() || result.nodes[std::get<0>(t).id].size() == 0)
		{
			throw CommonUtils::InvalidGraphException { std::string { "The graph has an edge between non-existant node(s) " + result.originalNodeName.at(std::get<0>(t).id) + (std::get<0>(t).end ? "+" : "-") + " and " + result.originalNodeName.at(std::get<1>(t).id) + (std::get<1>(t).end ? "+" : "-") } };
		}
		if (std::get<1>(t).id >= result.nodes.size() || result.nodes[std::get<1>(t).id].size() == 0)
		{
			throw CommonUtils::InvalidGraphException { std::string { "The graph has an edge between non-existant node(s) " + result.originalNodeName.at(std::get<0>(t).id) + (std::get<0>(t).end ? "+" : "-") + " and " + result.originalNodeName.at(std::get<1>(t).id) + (std::get<1>(t).end ? "+" : "-") } };
		}
	}
	std::sort(result.edges.begin(), result.edges.end(), [](std::tuple<NodePos, NodePos, size_t> left, std::tuple<NodePos, NodePos, size_t> right){
		if (std::get<0>(left) < std::get<0>(right)) return true;
		if (std::get<0>(left) > std::get<0>(right)) return false;
		if (std::get<1>(left) < std::get<1>(right)) return true;
		if (std::get<1>(left) > std::get<1>(right)) return false;
		if (std::get<2>(left) < std::get<2>(right)) return true;
		if (std::get<2>(left) > std::get<2>(right)) return false;
		return false;
	});
	if (result.edges.size() > 0)
	{
		for (size_t i = result.edges.size()-1; i > 0; i--)
		{
			if (result.edges[i] != result.edges[i-1]) continue;
			std::swap(result.edges[i], result.edges.back());
			result.edges.pop_back();
		}
	}
	// edges are not sorted anymore but doesn't matter as long as the order is deterministic
	return result;
}

std::string GfaGraph::OriginalNodeName(int nodeId) const
{
	assert(nodeId < originalNodeName.size());
	return originalNodeName[nodeId];
}

size_t GfaGraph::totalBp() const
{
	size_t result = 0;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		result += nodes[i].size();
	}
	return result;
}
