#include <limits>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include "vg.pb.h"
#include "CommonUtils.h"
#include "stream.hpp"

enum EndType
{
	LeftEdge,
	LeftOff,
	RightEdge,
	RightOff
};

enum NodeKeepingType
{
	KeepLeft,
	KeepRight,
	KeepAll
};

class PreGraph
{
public:
	class Edge
	{
	public:
		Edge(size_t from, bool fromStart, size_t to, bool toEnd) :
		from(from), to(to), fromStart(fromStart), toEnd(toEnd)
		{}
		size_t from;
		size_t to;
		bool fromStart;
		bool toEnd;
	};
	std::vector<Edge> edges;
	std::vector<std::string> nodeSequences;
};

void setKeepingType(const std::vector<std::set<size_t>>& goodEdges, const std::vector<std::set<size_t>>& badEdges, std::vector<bool>& hasKeepingType, std::vector<NodeKeepingType>& result, size_t node, NodeKeepingType type)
{
	std::vector<std::tuple<size_t, NodeKeepingType>> stack;
	stack.emplace_back(node, type);
	while (stack.size() > 0)
	{
		auto top = stack.back();
		stack.pop_back();
		node = std::get<0>(top);
		type = std::get<1>(top);
		bool madeKeepAll = false;
		if (hasKeepingType[node])
		{
			if (result[node] != type) result[node] = KeepAll;
			continue;
		}
		assert(type == KeepLeft || type == KeepRight);
		hasKeepingType[node] = true;
		result[node] = type;
		for (auto neighbor : goodEdges[node])
		{
			if (!hasKeepingType[neighbor]) continue;
			if (result[neighbor] == KeepAll) continue;
			if (result[neighbor] != result[node])
			{
				result[node] = KeepAll;
				madeKeepAll = true;
				break;
			}
		}
		if (madeKeepAll) continue;
		for (auto neighbor : badEdges[node])
		{
			if (!hasKeepingType[neighbor]) continue;
			if (result[neighbor] == KeepAll) continue;
			if (result[neighbor] == result[node])
			{
				result[node] = KeepAll;
				madeKeepAll = true;
				break;
			}
		}
		if (madeKeepAll) continue;
		for (auto neighbor : goodEdges[node])
		{
			if (hasKeepingType[neighbor]) continue;
			stack.emplace_back(neighbor, type);
		}
		for (auto neighbor : badEdges[node])
		{
			if (hasKeepingType[neighbor]) continue;
			assert(type == KeepLeft || type == KeepRight);
			stack.emplace_back(neighbor, type == KeepLeft ? KeepRight : KeepLeft);
		}
	}
}

std::vector<NodeKeepingType> getNodeKeepingTypes(const PreGraph& graph)
{
	std::vector<bool> hasKeepingType;
	std::vector<NodeKeepingType> result;
	hasKeepingType.resize(graph.nodeSequences.size(), false);
	result.resize(graph.nodeSequences.size());
	// result.resize(graph.nodeSequences.size(), KeepAll);
	// return result;
	{
		std::vector<bool> hasLeftEdge;
		std::vector<bool> hasRightEdge;
		hasLeftEdge.resize(graph.nodeSequences.size(), false);
		hasRightEdge.resize(graph.nodeSequences.size(), false);
		for (auto edge : graph.edges)
		{
			if (edge.fromStart)
			{
				hasLeftEdge[edge.from] = true;
			}
			else
			{
				hasRightEdge[edge.from] = true;
			}
			if (edge.toEnd)
			{
				hasRightEdge[edge.to] = true;
			}
			else
			{
				hasLeftEdge[edge.to] = true;
			}
		}
		for (size_t i = 0; i < result.size(); i++)
		{
			if (!hasLeftEdge[i] || !hasRightEdge[i])
			{
				result[i] = KeepAll;
				hasKeepingType[i] = true;
			}
		}
	}
	std::vector<std::set<size_t>> goodEdges;
	std::vector<std::set<size_t>> badEdges;
	goodEdges.resize(graph.nodeSequences.size());
	badEdges.resize(graph.nodeSequences.size());
	for (auto edge : graph.edges)
	{
		if (edge.fromStart == edge.toEnd)
		{
			goodEdges[edge.from].insert(edge.to);
			goodEdges[edge.to].insert(edge.from);
		}
		else
		{
			badEdges[edge.from].insert(edge.to);
			badEdges[edge.to].insert(edge.from);
		}
	}
	result.resize(graph.nodeSequences.size());
	hasKeepingType.resize(graph.nodeSequences.size(), false);
	for (size_t i = 0; i < graph.nodeSequences.size(); i++)
	{
		if (!hasKeepingType[i]) setKeepingType(goodEdges, badEdges, hasKeepingType, result, i, KeepLeft);
	}
	return result;
}

PreGraph loadGraphFromGfa(std::string filename)
{
	std::unordered_map<int, std::string> nodeSequences;
	std::ifstream file {filename};
	PreGraph result;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		if (line.size() == 0) continue;
		if (line[0] == 'S')
		{
			std::stringstream sstr {line};
			int id;
			std::string sequence;
			std::string empty;
			sstr >> empty >> id >> sequence;
			assert(sequence.size() > 0);
#ifndef NDEBUG
			for (size_t i = 0; i < sequence.size(); i++)
			{
				assert(sequence[i] == 'a' || sequence[i] == 'A' || sequence[i] == 't' || sequence[i] == 'T' || sequence[i] == 'c' || sequence[i] == 'C' || sequence[i] == 'g' || sequence[i] == 'G');
			}
#endif
			assert(nodeSequences.count(id) == 0);
			nodeSequences[id] = sequence;
		}
		else if (line[0] == 'L')
		{
			std::stringstream sstr {line};
			int fromid, toid;
			std::string fromstart, toend;
			std::string empty;
			sstr >> empty >> fromid >> fromstart >> toid >> toend;
			assert(fromstart == "+" || fromstart == "-");
			assert(toend == "+" || toend == "-");
			result.edges.emplace_back(fromid, fromstart == "-", toid, toend == "-");
		}
	}
	result.nodeSequences.resize(nodeSequences.size());
	for (auto iter : nodeSequences)
	{
		assert(iter.first < result.nodeSequences.size());
		assert(iter.second.size() > 0);
		result.nodeSequences[iter.first] = iter.second;
		result.nodeSequences[iter.first].shrink_to_fit();
	}
	result.edges.shrink_to_fit();
	return result;
}

std::pair<size_t, bool> getNewIndexAndDirection(size_t oldNodeSize, int kmin1, size_t oldNodeId, bool oldNodeEnd, bool oldNodeOff)
{
	if (oldNodeEnd && !oldNodeOff)
	{
		return std::make_pair(oldNodeId * 3 + 1, true);
	}
	if (!oldNodeEnd && !oldNodeOff)
	{
		return std::make_pair(oldNodeId * 3, false);
	}
	if (oldNodeEnd && oldNodeOff)
	{
		if (oldNodeSize > 2 * kmin1)
		{
			return std::make_pair(oldNodeId * 3 + 2, true);
		}
		else
		{
			return std::make_pair(oldNodeId * 3, true);
		}
	}
	if (!oldNodeEnd & oldNodeOff)
	{
		if (oldNodeSize > 2 * kmin1)
		{
			return std::make_pair(oldNodeId * 3 + 2, false);
		}
		else
		{
			return std::make_pair(oldNodeId * 3 + 1, false);
		}
	}
	assert(false);
	return std::make_pair(0, false);
}

PreGraph bluntifyViaKeepingTypes(const PreGraph& graph, const std::vector<NodeKeepingType>& keepingType, int k)
{
	assert(k > 1);
	const int kmin1 = k - 1;
	PreGraph result;
	std::vector<bool> hasLeftPart;
	std::vector<bool> hasRightPart;
	std::vector<bool> hasMiddlePart;
	hasLeftPart.resize(graph.nodeSequences.size(), false);
	hasRightPart.resize(graph.nodeSequences.size(), false);
	hasMiddlePart.resize(graph.nodeSequences.size(), false);
	for (size_t i = 0; i < graph.nodeSequences.size(); i++)
	{
		result.nodeSequences.push_back("");
		result.nodeSequences.push_back("");
		result.nodeSequences.push_back("");
		const size_t seqsize = graph.nodeSequences[i].size();
		if (seqsize < 2 * kmin1)
		{
			if (keepingType[i] == KeepLeft || keepingType[i] == KeepAll)
			{
				result.nodeSequences[i * 3] = graph.nodeSequences[i].substr(0, seqsize - kmin1);
				hasLeftPart[i] = true;
			}
			if (keepingType[i] == KeepRight || keepingType[i] == KeepAll)
			{
				result.nodeSequences[i * 3 + 1] = graph.nodeSequences[i].substr(kmin1);
				hasRightPart[i] = true;
			}
			if (keepingType[i] == KeepAll)
			{
				result.nodeSequences[i * 3 + 2] = graph.nodeSequences[i].substr(seqsize - kmin1, 2 * kmin1 - seqsize);
				hasMiddlePart[i] = true;
			}
		}
		else if (seqsize == 2 * kmin1)
		{
			if (keepingType[i] == KeepLeft || keepingType[i] == KeepAll)
			{
				result.nodeSequences[i * 3] = graph.nodeSequences[i].substr(0, kmin1);
				hasLeftPart[i] = true;
			}
			if (keepingType[i] == KeepRight || keepingType[i] == KeepAll)
			{
				result.nodeSequences[i * 3 + 1] = graph.nodeSequences[i].substr(seqsize - kmin1);
				hasRightPart[i] = true;
			}
		}
		else
		{
			assert(seqsize > 2 * kmin1);
			if (keepingType[i] == KeepLeft || keepingType[i] == KeepAll)
			{
				result.nodeSequences[i * 3] = graph.nodeSequences[i].substr(0, kmin1);
				hasLeftPart[i] = true;
			}
			if (keepingType[i] == KeepRight || keepingType[i] == KeepAll)
			{
				result.nodeSequences[i * 3 + 1] = graph.nodeSequences[i].substr(seqsize - kmin1);
				hasRightPart[i] = true;
			}
			result.nodeSequences[i * 3 + 2] = graph.nodeSequences[i].substr(kmin1, seqsize - 2 * kmin1);
			hasMiddlePart[i] = true;
		}
		if (hasLeftPart[i] && hasRightPart[i])
		{
			assert(result.nodeSequences[i * 3].size() == result.nodeSequences[i * 3 + 1].size());
		}
		if (hasLeftPart[i] && hasMiddlePart[i])
		{
			result.edges.emplace_back(i * 3, false, i * 3 + 2, false);
		}
		if (hasMiddlePart[i] && hasRightPart[i])
		{
			result.edges.emplace_back(i * 3 + 2, false, i * 3 + 1, false);
		}
		if (seqsize == 2 * kmin1 && hasLeftPart[i] && hasRightPart[i])
		{
			result.edges.emplace_back(i * 3, false, i * 3 + 1, false);
		}
	}
	for (auto edge : graph.edges)
	{
		auto newfrom = getNewIndexAndDirection(graph.nodeSequences[edge.from].size(), kmin1, edge.from, !edge.fromStart, false);
		auto newto = getNewIndexAndDirection(graph.nodeSequences[edge.to].size(), kmin1, edge.to, edge.toEnd, true);
		if (newfrom.first % 3 == 0 && !hasLeftPart[edge.from]) continue;
		if (newfrom.first % 3 == 1 && !hasRightPart[edge.from]) continue;
		if (newfrom.first % 3 == 2 && !hasMiddlePart[edge.from]) continue;
		if (newto.first % 3 == 0 && !hasLeftPart[edge.to]) continue;
		if (newto.first % 3 == 1 && !hasRightPart[edge.to]) continue;
		if (newto.first % 3 == 2 && !hasMiddlePart[edge.to]) continue;
		result.edges.emplace_back(newfrom.first, !newfrom.second, newto.first, newto.second);
	}
	for (auto edge : graph.edges)
	{
		auto newfrom = getNewIndexAndDirection(graph.nodeSequences[edge.from].size(), kmin1, edge.from, !edge.fromStart, true);
		auto newto = getNewIndexAndDirection(graph.nodeSequences[edge.to].size(), kmin1, edge.to, edge.toEnd, false);
		if (newfrom.first % 3 == 0 && !hasLeftPart[edge.from]) continue;
		if (newfrom.first % 3 == 1 && !hasRightPart[edge.from]) continue;
		if (newfrom.first % 3 == 2 && !hasMiddlePart[edge.from]) continue;
		if (newto.first % 3 == 0 && !hasLeftPart[edge.to]) continue;
		if (newto.first % 3 == 1 && !hasRightPart[edge.to]) continue;
		if (newto.first % 3 == 2 && !hasMiddlePart[edge.to]) continue;
		result.edges.emplace_back(newfrom.first, !newfrom.second, newto.first, newto.second);
	}
	return result;
}

void writeGFA(const PreGraph& graph, std::string filename)
{
	std::ofstream file {filename};
	//start at 1 because 0 is not a valid node id in vg
	const int off = 1;
	for (int i = 0; i < graph.nodeSequences.size(); i++)
	{
		if (graph.nodeSequences[i].size() == 0) continue;
		file << "S\t" << (i + off) << "\t" << graph.nodeSequences[i] << std::endl;
	}
	for (int i = 0; i < graph.edges.size(); i++)
	{
		assert(graph.nodeSequences[graph.edges[i].from].size() > 0);
		assert(graph.nodeSequences[graph.edges[i].to].size() > 0);
		file << "L\t" << (graph.edges[i].from + off) << "\t" << (graph.edges[i].fromStart ? "-" : "+") << "\t" << (graph.edges[i].to + off) << "\t" << (graph.edges[i].toEnd ? "-" : "+") << "\t0M" << std::endl;
	}
}

int main(int argc, char** argv)
{
	int k = std::stoi(argv[1]);
	std::string inFile = argv[2];
	std::string outFile = argv[3];
	auto graph = loadGraphFromGfa(inFile);
	auto keepingTypes = getNodeKeepingTypes(graph);
	size_t left, right, all;
	left = 0;
	right = 0;
	all = 0;
	for (size_t i = 0; i < keepingTypes.size(); i++)
	{
		switch(keepingTypes[i])
		{
		case KeepLeft:
			left++;
			break;
		case KeepRight:
			right++;
			break;
		case KeepAll:
			all++;
			break;
		default:
			assert(false);
		}
	}
	std::cerr << "left: " << left << " right: " << right << " all: " << all << std::endl;
	auto result = bluntifyViaKeepingTypes(graph, keepingTypes, k);
	writeGFA(result, outFile);
}
