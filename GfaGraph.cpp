#include <fstream>
#include <sstream>
#include "GfaGraph.h"
#include "ThreadReadAssertion.h"

GfaGraph::GfaGraph() :
nodes(),
edges(),
edgeOverlap(-1)
{
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
			Node newnode;
			newnode.id = id;
			newnode.sequence = seq;
			result.nodes.push_back(newnode);
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
			Edge newedge;
			newedge.from = from;
			newedge.to = to;
			newedge.fromStart = fromstart == "-";
			newedge.toEnd = toend == "-";
			result.edges.push_back(newedge);
		}
	}
	return result;
}

