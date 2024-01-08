#include "CommonUtils.h"
#include "stream.hpp"

namespace CommonUtils
{
	InvalidGraphException::InvalidGraphException(const char* c) : std::runtime_error(c) 
	{
	}

	InvalidGraphException::InvalidGraphException(std::string c) : std::runtime_error(c) 
	{
	}

	void mergeGraphs(vg::Graph& graph, const vg::Graph& part)
	{
		for (int i = 0; i < part.node_size(); i++)
		{
			auto node = graph.add_node();
			node->set_id(part.node(i).id());
			node->set_sequence(part.node(i).sequence());
			node->set_name(part.node(i).name());
		}
		for (int i = 0; i < part.edge_size(); i++)
		{
			auto edge = graph.add_edge();
			edge->set_from(part.edge(i).from());
			edge->set_to(part.edge(i).to());
			edge->set_from_start(part.edge(i).from_start());
			edge->set_to_end(part.edge(i).to_end());
			edge->set_overlap(part.edge(i).overlap());
		}
	}

	vg::Graph LoadVGGraph(std::string filename)
	{
		vg::Graph result;
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Graph&)> lambda = [&result](vg::Graph& g) {
			mergeGraphs(result, g);
		};
		stream::for_each(graphfile, lambda);
		if (result.node_size() == 0)
		{
			std::cerr << "Failed to load vg graph." << std::endl;
			std::cerr << "This graph might be in a format other than protobuf. GraphAligner only supports vg graphs in protobuf format." << std::endl;
			std::cerr << "Check the graph format with:" << std::endl;
			std::cerr << "vg stats -F " << filename << std::endl;
			std::cerr << "or try converting the graph to GFA" << std::endl;
			std::abort();
		}
		return result;
	}

	std::vector<vg::Alignment> LoadVGAlignments(std::string filename)
	{
		std::vector<vg::Alignment> result;
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda = [&result](vg::Alignment& g) {
			result.push_back(g);
		};
		stream::for_each(graphfile, lambda);
		return result;
	}

	vg::Alignment LoadVGAlignment(std::string filename)
	{
		vg::Alignment result;
		std::ifstream graphfile { filename, std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda = [&result](vg::Alignment& g) {
			result = g;
		};
		stream::for_each(graphfile, lambda);
		return result;
	}

	std::string ReverseComplement(std::string str)
	{
		std::string result;
		result.reserve(str.size());
		for (int i = str.size()-1; i >= 0; i--)
		{
			result += Complement(str[i]);
		}
		return result;
	}

	char Complement(char c)
	{
		switch (c)
		{
			case 'A':
			case 'a':
				return 'T';
			case 'C':
			case 'c':
				return 'G';
			case 'T':
			case 't':
				return 'A';
			case 'G':
			case 'g':
				return 'C';
			case 'N':
			case 'n':
				return 'N';
			case 'U':
			case 'u':
				return 'A';
			case 'R':
			case 'r':
				return 'Y';
			case 'Y':
			case 'y':
				return 'R';
			case 'K':
			case 'k':
				return 'M';
			case 'M':
			case 'm':
				return 'K';
			case 'S':
			case 's':
				return 'S';
			case 'W':
			case 'w':
				return 'W';
			case 'B':
			case 'b':
				return 'V';
			case 'V':
			case 'v':
				return 'B';
			case 'D':
			case 'd':
				return 'H';
			case 'H':
			case 'h':
				return 'D';
			case '-':
				return '-';
			default:
				assert(false);
				return 'N';
	}
	}

}

BufferedWriter::BufferedWriter() : stream(nullptr) {};
BufferedWriter::BufferedWriter(std::ostream& stream) : stream(&stream) {};
BufferedWriter& BufferedWriter::operator<<(FlushClass)
{
	if (stream == nullptr) return *this;
	flush();
	return *this;
}
void BufferedWriter::flush()
{
	if (stream == nullptr) return;
	stringstream << std::endl;
	(*stream) << stringstream.str();
	stringstream.str("");
}
bool BufferedWriter::inputDiscarded() const
{
	return stream == nullptr;
}