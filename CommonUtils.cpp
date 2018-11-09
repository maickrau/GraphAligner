#include "CommonUtils.h"
#include "stream.hpp"

namespace CommonUtils
{
	namespace inner
	{
		//an overlap which is larger than the fraction cutoff of the smaller alignment means the alignments are incompatible
		//eg alignments 12000bp and 15000bp, overlap of 12000*0.05 = 600bp means they are incompatible
		const float OverlapIncompatibleFractionCutoff = 0.05;

		//longer alignments are better
		bool alignmentLengthCompare(const vg::Alignment* const left, const vg::Alignment* const right)
		{
			return left->sequence().size() > right->sequence().size();
		}

		//lower scores are better
		bool alignmentScoreCompare(const vg::Alignment* const left, const vg::Alignment* const right)
		{
			return left->score() < right->score();
		}

		bool alignmentIncompatible(const vg::Alignment* const left, const vg::Alignment* const right)
		{
			auto minOverlapLen = std::min(left->sequence().size(), right->sequence().size()) * OverlapIncompatibleFractionCutoff;
			assert(left->query_position() >= 0);
			assert(right->query_position() >= 0);
			size_t leftStart = left->query_position();
			size_t leftEnd = leftStart + left->sequence().size();
			size_t rightStart = right->query_position();
			size_t rightEnd = rightStart + right->sequence().size();
			if (leftStart > rightStart)
			{
				std::swap(leftStart, rightStart);
				std::swap(leftEnd, rightEnd);
			}
			int overlap = 0;
			assert(leftStart <= rightStart);
			if (leftEnd > rightStart) overlap = leftEnd - rightStart;
			return overlap > minOverlapLen;
		}
	}

	std::vector<vg::Alignment> SelectAlignments(std::vector<vg::Alignment> alns, size_t maxnum)
	{
		return SelectAlignments(alns, maxnum, [](const vg::Alignment& aln) { return &aln; });
	}

	std::vector<vg::Alignment*> SelectAlignments(std::vector<vg::Alignment*> alns, size_t maxnum)
	{
		return SelectAlignments(alns, maxnum, [](vg::Alignment* aln) { return aln; });
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
			switch (str[i])
			{
				case 'A':
				case 'a':
				result += 'T';
				break;
				case 'C':
				case 'c':
				result += 'G';
				break;
				case 'T':
				case 't':
				result += 'A';
				break;
				case 'G':
				case 'g':
				result += 'C';
				break;
				case 'N':
				case 'n':
				result += 'N';
				break;
				case 'U':
				case 'u':
				result += 'A';
				break;
				case 'R':
				case 'r':
				result += 'Y';
				break;
				case 'Y':
				case 'y':
				result += 'R';
				break;
				case 'K':
				case 'k':
				result += 'M';
				break;
				case 'M':
				case 'm':
				result += 'K';
				break;
				case 'S':
				case 's':
				result += 'S';
				break;
				case 'W':
				case 'w':
				result += 'W';
				break;
				case 'B':
				case 'b':
				result += 'V';
				break;
				case 'V':
				case 'v':
				result += 'B';
				break;
				case 'D':
				case 'd':
				result += 'H';
				break;
				case 'H':
				case 'h':
				result += 'D';
				default:
				assert(false);
			}
		}
		return result;
	}

}

BufferedWriter::BufferedWriter() : stream(nullptr) {};
BufferedWriter::BufferedWriter(std::ostream& stream) : stream(&stream) {};
BufferedWriter& BufferedWriter::operator<<(FlushClass f)
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
