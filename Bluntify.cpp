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

class BluntifyGraph
{
public:
	class Edge
	{
	public:
		Edge(size_t from, size_t to, EndType fromPos, EndType toPos) :
		from(from),
		to(to),
		fromPos(fromPos),
		toPos(toPos)
		{
		}
		size_t from;
		size_t to;
		EndType fromPos;
		EndType toPos;
	};
	std::vector<std::string> nodeSequences;
	std::vector<Edge> edges;
};

BluntifyGraph loadGraphFromGfa(std::string filename, int k)
{
	std::unordered_map<int, std::string> nodeSequences;
	std::ifstream file {filename};
	BluntifyGraph result;
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
			EndType frompos, topos;
			frompos = fromstart == "-" ? LeftEdge : RightEdge;
			topos = toend == "-" ? RightOff : LeftOff;
			result.edges.emplace_back(fromid, toid, frompos, topos);
			frompos = fromstart == "-" ? LeftOff : RightOff;
			topos = toend == "-" ? RightEdge : LeftEdge;
			result.edges.emplace_back(fromid, toid, frompos, topos);
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

std::pair<size_t, bool> getNewIndexAndDirection(const std::vector<bool>& smallerThan2k, const std::vector<bool>& biggerThan2k, size_t node, EndType pos)
{
	assert(node < smallerThan2k.size());
	assert(node < biggerThan2k.size());
	assert(!(smallerThan2k[node] && biggerThan2k[node]));
	if (pos == LeftEdge)
	{
		return std::make_pair(node * 3 + 1, false);
	}
	else if (pos == RightEdge)
	{
		return std::make_pair(node * 3 + 2, true);
	}
	else if (pos == LeftOff)
	{
		if (smallerThan2k[node])
		{
			return std::make_pair(node * 3 + 2, false);
		}
		else if (biggerThan2k[node])
		{
			return std::make_pair(node * 3 + 3, false);
		}
		else
		{
			return std::make_pair(node * 3 + 2, false);
		}
	}
	else if (pos == RightOff)
	{
		if (smallerThan2k[node])
		{
			return std::make_pair(node * 3 + 1, true);
		}
		else if (biggerThan2k[node])
		{
			return std::make_pair(node * 3 + 3, true);
		}
		else
		{
			return std::make_pair(node * 3 + 1, true);
		}
	}
	assert(false);
	return std::make_pair(0, true);
}

vg::Graph convertBluntifyToVg(const BluntifyGraph& graph, int k)
{
	vg::Graph result;
	std::vector<bool> smallerThan2k;
	std::vector<bool> biggerThan2k;
	smallerThan2k.resize(graph.nodeSequences.size(), false);
	biggerThan2k.resize(graph.nodeSequences.size(), false);
	for (size_t i = 0; i < graph.nodeSequences.size(); i++)
	{
		auto& seq = graph.nodeSequences[i];
		assert(seq.size() > 0);
		if (seq.size() < 2 * (k - 1)) smallerThan2k[i] = true;
		if (seq.size() > 2 * (k - 1)) biggerThan2k[i] = true;
		auto leftnode = result.add_node();
		assert(leftnode != nullptr);
		leftnode->set_id(i * 3 + 1);
		leftnode->set_name(std::to_string(i * 3 + 1));
		if (smallerThan2k[i])
		{
			leftnode->set_sequence(seq.substr(0, seq.size() - (k - 1)));
		}
		else
		{
			leftnode->set_sequence(seq.substr(0, k-1));
		}
		assert(leftnode->sequence().size() > 0);
		auto rightnode = result.add_node();
		assert(rightnode != nullptr);
		rightnode->set_id(i * 3 + 2);
		rightnode->set_name(std::to_string(i * 3 + 2));
		if (smallerThan2k[i])
		{
			rightnode->set_sequence(seq.substr(k - 1));
		}
		else
		{
			rightnode->set_sequence(seq.substr(seq.size() - (k - 1)));
		}
		assert(rightnode->sequence().size() > 0);
		assert(leftnode->sequence().size() == rightnode->sequence().size());
		if (smallerThan2k[i])
		{
			auto middlenode = result.add_node();
			assert(middlenode != nullptr);
			middlenode->set_id(i * 3 + 3);
			middlenode->set_name(std::to_string(i * 3 + 3));
			assert(seq.size() >= k - 1);
			assert(seq.size() - (k - 1) + 2 * (k - 1) - seq.size() < seq.size());
			middlenode->set_sequence(seq.substr(seq.size() - (k - 1), 2 * (k - 1) - seq.size()));
			assert(middlenode->sequence().size() > 0);
			assert(leftnode->sequence().size() + rightnode->sequence().size() + middlenode->sequence().size() == seq.size());
			auto edge1 = result.add_edge();
			assert(edge1 != nullptr);
			edge1->set_from(i * 3 + 1);
			edge1->set_to(i * 3 + 3);
			edge1->set_from_start(false);
			edge1->set_to_end(false);
			edge1->set_overlap(0);
			auto edge2 = result.add_edge();
			assert(edge2 != nullptr);
			edge2->set_from(i * 3 + 3);
			edge2->set_to(i * 3 + 2);
			edge2->set_from_start(false);
			edge2->set_to_end(false);
			edge2->set_overlap(0);
		}
		else if (biggerThan2k[i])
		{
			auto middlenode = result.add_node();
			assert(middlenode != nullptr);
			middlenode->set_id(i * 3 + 3);
			middlenode->set_name(std::to_string(i * 3 + 3));
			assert(k - 1 + seq.size() - 2 * (k - 1) < seq.size());
			middlenode->set_sequence(seq.substr(k - 1, seq.size() - 2 * (k - 1)));
			assert(middlenode->sequence().size() > 0);
			assert(leftnode->sequence().size() + rightnode->sequence().size() + middlenode->sequence().size() == seq.size());
			auto edge1 = result.add_edge();
			assert(edge1 != nullptr);
			edge1->set_from(i * 3 + 1);
			edge1->set_to(i * 3 + 3);
			edge1->set_from_start(false);
			edge1->set_to_end(false);
			edge1->set_overlap(0);
			auto edge2 = result.add_edge();
			assert(edge2 != nullptr);
			edge2->set_from(i * 3 + 3);
			edge2->set_to(i * 3 + 2);
			edge2->set_from_start(false);
			edge2->set_to_end(false);
			edge2->set_overlap(0);
		}
		else
		{
			auto edge1 = result.add_edge();
			assert(edge1 != nullptr);
			edge1->set_from(i * 3 + 1);
			edge1->set_to(i * 3 + 2);
			edge1->set_from_start(false);
			edge1->set_to_end(false);
			edge1->set_overlap(0);
		}
	}
	for (size_t i = 0; i < graph.edges.size(); i++)
	{
		auto left = getNewIndexAndDirection(smallerThan2k, biggerThan2k, graph.edges[i].from, graph.edges[i].fromPos);
		auto right = getNewIndexAndDirection(smallerThan2k, biggerThan2k, graph.edges[i].to, graph.edges[i].toPos);
		assert(left.first > 0);
		assert(right.first > 0);
		assert(left.first <= graph.nodeSequences.size() * 3 + 3);
		assert(right.first <= graph.nodeSequences.size() * 3 + 3);
		left.second = !left.second;
		auto edge = result.add_edge();
		edge->set_from(left.first);
		edge->set_from_start(left.second);
		edge->set_to(right.first);
		edge->set_to_end(right.second);
		edge->set_overlap(0);
	}
	assert(result.node_size() >= graph.nodeSequences.size() * 2);
	assert(result.node_size() <= graph.nodeSequences.size() * 3);
	assert(result.edge_size() >= graph.edges.size());
	assert(result.edge_size() <= graph.edges.size() + graph.nodeSequences.size() * 2);
	return result;
}

void writeGraph(const vg::Graph& graph, std::string filename)
{
	std::vector<vg::Graph> writeVector;
	writeVector.push_back(graph);
	std::ofstream graphOut { filename, std::ios::out | std::ios::binary };
	stream::write_buffered(graphOut, writeVector, 0);
}

void writeGFA(const vg::Graph& graph, std::string filename)
{
	std::ofstream file {filename};
	for (int i = 0; i < graph.node_size(); i++)
	{
		file << "S\t" << graph.node(i).id() << "\t" << graph.node(i).sequence() << std::endl;
	}
	for (int i = 0; i < graph.edge_size(); i++)
	{
		file << "L\t" << graph.edge(i).from() << "\t" << (graph.edge(i).from_start() ? "-" : "+") << "\t" << graph.edge(i).to() << "\t" << (graph.edge(i).to_end() ? "-" : "+") << "\t0M" << std::endl;
	}
}

int main(int argc, char** argv)
{
	int k = std::stoi(argv[1]);
	std::string inFile = argv[2];
	std::string outFile = argv[3];
	auto graph = convertBluntifyToVg(loadGraphFromGfa(inFile, k), k);
	writeGFA(graph, outFile);
}
