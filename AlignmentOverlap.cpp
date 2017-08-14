#include <algorithm>
#include <vector>
#include <set>
#include "vg.pb.h"
#include "stream.hpp"
#include "CommonUtils.h"
#include "fastqloader.h"

const int LEFT = 1;
const int UP = 2;
const int DIAGONAL = 3;

class Overlap
{
public:
	std::string readname1;
	std::string readname2;
	int length1;
	int length2;
	bool backward1;
	bool backward2;
};

class NodeMovement
{
public:
	int nodeId;
	bool backwards;
	int offset;
	int length;
	bool operator==(const NodeMovement& other) const
	{
		return nodeId == other.nodeId && backwards == other.backwards && ((offset <= other.offset && offset + length >= other.offset + other.length) || (other.offset <= offset && other.offset + other.length >= offset + length));
	}
	bool operator<(const NodeMovement& other) const
	{
		return nodeId < other.nodeId || (nodeId == other.nodeId && backwards && !other.backwards);
	}
};

std::vector<NodeMovement> reverse(const std::vector<NodeMovement>& vec)
{
	std::vector<NodeMovement> result { vec };
	std::reverse(result.begin(), result.end());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i].backwards = !result[i].backwards;
	}
	return result;
}

std::vector<NodeMovement> getNodeMovements(const vg::Alignment& alignment, const std::map<int, int>& nodeSizes)
{
	std::vector<NodeMovement> result;
	for (int i = 0; i < alignment.path().mapping_size(); i++)
	{
		NodeMovement node;
		node.nodeId = alignment.path().mapping(i).position().node_id();
		node.backwards = alignment.path().mapping(i).position().is_reverse();
		node.length = alignment.path().mapping(i).edit(0).from_length();
		node.offset = alignment.path().mapping(i).position().offset();
		result.push_back(node);
	}
	return result;
}

Overlap backtrace(const std::vector<std::vector<int>>& direction, std::string readname1, std::string readname2, const std::vector<NodeMovement>& read1, const std::vector<NodeMovement>& read2, int starti, int startj, bool backward2)
{
	Overlap result;
	result.readname1 = readname1;
	result.readname2 = readname2;
	result.length1 = 0;
	result.length2 = 0;
	result.backward1 = false;
	result.backward2 = backward2;
	auto i = starti;
	auto j = startj;
	while (i != 0 && j != 0)
	{
		switch(direction[i][j])
		{
			case LEFT:
			j--;
			break;
			case UP:
			i--;
			break;
			case DIAGONAL:
			i--;
			j--;
			break;
			default:
			assert(false);
		}
	}
	auto endi = i;
	auto endj = j;
	if (starti == read1.size() && endi == 0) return result;
	if (startj == read2.size() && endj == 0) return result;
	for (i = endi; i <= starti; i++)
	{
		result.length1 += read1[i-1].length;
	}
	for (j = endj; j <= startj; j++)
	{
		result.length2 += read2[j-1].length;
	}
	// assert(result.length1 == result.length2);
	assert((endi == 0 && startj == read2.size()) || (endj == 0 && starti == read1.size()));
	if (endi == 0 && startj == read2.size())
	{
		std::swap(result.readname1, result.readname2);
		std::swap(result.length1, result.length2);
		std::swap(result.backward1, result.backward2);
	}
	return result;
}

std::vector<Overlap> getExactOverlaps(const std::string& readname1, const std::vector<NodeMovement>& read1, const std::string& readname2, const std::vector<NodeMovement>& read2, double minMatchFraction, double minSizeFraction, bool backward2)
{
	int read1size = 0;
	int read2size = 0;
	for (size_t i = 0; i < read1.size(); i++)
	{
		read1size += read1[i].length;
	}
	for (size_t i = 0; i < read2.size(); i++)
	{
		read2size += read2[i].length;
	}
	int minMatchSize = std::min(read1size * minSizeFraction, read2size * minSizeFraction);
	std::vector<Overlap> result;
	for (size_t i = 0; i < read1.size(); i++)
	{
		if (read1.size() >= read2.size() && i <= read1.size() - read2.size())
		{
			i = read1.size()-read2.size()+1;
		}
		bool match = true;
		int length = 0;
		for (size_t k = 0; k < read1.size()-i; k++)
		{
			assert(k < read2.size());
			if (read1[i+k] == read2[k])
			{
				length += std::min(read1[i+k].length, read2[k].length);
			}
			else
			{
				match = false;
				break;
			}
		}
		assert(length <= read1size);
		assert(length <= read2size);
		if (match && length > minMatchSize)
		{
			Overlap readmatch;
			readmatch.readname1 = readname1;
			readmatch.readname2 = readname2;
			readmatch.length1 = length;
			readmatch.length2 = length;
			readmatch.backward1 = false;
			readmatch.backward2 = backward2;
			result.push_back(readmatch);
			break;
		}
	}
	return result;
}

std::vector<Overlap> getOverlaps(const std::string& readname1, const std::vector<NodeMovement>& read1, const std::string& readname2, const std::vector<NodeMovement>& read2, double minMatchFraction, double minSizeFraction, bool backward2)
{
	int read1size = 0;
	int read2size = 0;
	for (size_t i = 0; i < read1.size(); i++)
	{
		read1size += read1[i].length;
	}
	for (size_t i = 0; i < read2.size(); i++)
	{
		read2size += read2[i].length;
	}
	int minMatchSize = std::min(read1size * minSizeFraction, read2size * minSizeFraction);
	std::vector<std::vector<int>> mismatches;
	std::vector<std::vector<int>> direction;
	std::vector<std::vector<int>> length;
	mismatches.resize(read1.size()+1);
	direction.resize(read1.size()+1);
	length.resize(read1.size()+1);
	for (size_t i = 0; i < read1.size()+1; i++)
	{
		mismatches[i].resize(read2.size()+1, std::numeric_limits<int>::min());
		direction[i].resize(read2.size()+1, 0);
		length[i].resize(read2.size()+1, 0);
	}
	for (size_t i = 0; i <= read1.size(); i++)
	{
		mismatches[i][0] = 0;
	}
	for (size_t i = 0; i <= read2.size(); i++)
	{
		mismatches[0][i] = 0;
	}
	for (size_t i = 1; i <= read1.size(); i++)
	{
		for (size_t j = 1; j <= read2.size(); j++)
		{
			mismatches[i][j] = mismatches[i-1][j]+read1[i].length;
			length[i][j] = length[i-1][j]+read1[i].length;
			direction[i][j] = UP;
			if (mismatches[i][j-1]+read2[j].length < mismatches[i][j])
			{
				mismatches[i][j] = mismatches[i][j-1]+read2[j].length;
				length[i][j] = length[i][j-1]+read2[j].length;
				direction[i][j] = LEFT;
			}
			if (read1[i-1] == read2[j-1])
			{
				if (mismatches[i-1][j-1] < mismatches[i][j])
				{
					mismatches[i][j] = mismatches[i-1][j-1];
					length[i][j] = length[i-1][j-1] + std::max(read1[i-1].length, read2[j-1].length);
					direction[i][j] = DIAGONAL;
				}
			}
			else
			{
				if (mismatches[i-1][j-1] + std::max(read1[i-1].length, read2[j-1].length) < mismatches[i][j])
				{
					mismatches[i][j] = mismatches[i-1][j-1] + std::max(read1[i-1].length, read2[j-1].length);
					length[i][j] = length[i-1][j-1] + std::max(read1[i-1].length, read2[j-1].length);
					direction[i][j] = DIAGONAL;
				}
			}
		}
	}
	std::vector<Overlap> result;
	int row = read2.size();
	for (size_t i = read1.size(); i > 0; i--)
	{
		assert(mismatches[i][row] <= length[i][row]);
		if (length[i][row] >= minMatchSize && 1.0 - (double)mismatches[i][row] / (double)length[i][row] >= minMatchFraction)
		{
			auto overlap = backtrace(direction, readname1, readname2, read1, read2, i, row, backward2);
			if (overlap.length1 == 0 || overlap.length2 == 0) continue;
			result.push_back(overlap);
			break;
		}
	}
	int col = read1.size();
	for (size_t j = read2.size(); j > 0; j--)
	{
		assert(mismatches[col][j] <= length[col][j]);
		if (length[col][j] >= minMatchSize && 1.0 - (double)mismatches[col][j] / (double)length[col][j] >= minMatchFraction)
		{
			auto overlap = backtrace(direction, readname1, readname2, read1, read2, col, j, backward2);
			if (overlap.length1 == 0 || overlap.length2 == 0) continue;
			result.push_back(overlap);
			break;
		}
	}
	return result;
}

bool alignmentPossible(const std::vector<int>& alignmentSizes, const std::vector<std::vector<NodeMovement>>& alignmentNodes, size_t first, size_t second, double minSizeFraction)
{
	int size = 0;
	int i = 0;
	int j = 0;
	auto minSize = std::min(alignmentSizes[first], alignmentSizes[second]) * minSizeFraction;
	while (i < alignmentNodes[first].size() && j < alignmentNodes[second].size())
	{
		if (alignmentNodes[first][i].nodeId < alignmentNodes[second][j].nodeId)
		{
			i++;
		}
		else if (alignmentNodes[second][j].nodeId < alignmentNodes[first][i].nodeId)
		{
			j++;
		}
		else
		{
			size += alignmentNodes[first][i].length;
			i++;
			j++;
			if (size >= minSize) return true;
		}
	}
	return false;
}

int main(int argc, char** argv)
{
	auto graph = CommonUtils::LoadVGGraph(argv[1]);
	auto reads = loadFastqFromFile(argv[5]);
	std::vector<vg::Alignment> alignments;
	std::map<int, int> nodeSizes;
	for (int i = 0; i < graph.node_size(); i++)
	{
		nodeSizes[graph.node(i).id()] = graph.node(i).sequence().size();
	}
	double minMatchFraction = std::stod(argv[3]);
	double minSizeFraction = std::stod(argv[4]);
	{
		std::ifstream graphfile { argv[2], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> lambda = [&alignments](vg::Alignment& g) {
			alignments.push_back(g);
		};
		stream::for_each(graphfile, lambda);
	}
	std::vector<int> alignmentSizes;
	std::vector<std::vector<NodeMovement>> alignmentNodeMovements;
	std::vector<std::vector<NodeMovement>> alignmentNodeComparison;
	alignmentSizes.resize(alignments.size(), 0);
	alignmentNodeMovements.resize(alignments.size());
	alignmentNodeComparison.resize(alignments.size());
	for (size_t i = 0; i < alignments.size(); i++)
	{
		alignmentNodeMovements[i] = getNodeMovements(alignments[i], nodeSizes);
		alignmentSizes[i] = 0;
		for (size_t j = 0; j < alignmentNodeMovements[i].size(); j++)
		{
			alignmentSizes[i] += alignmentNodeMovements[i][j].length;
			NodeMovement node = alignmentNodeMovements[i][j];
			node.backwards = false;
			alignmentNodeComparison[i].push_back(node);
		}
		std::sort(alignmentNodeComparison[i].begin(), alignmentNodeComparison[i].end());
	}
	std::vector<Overlap> validOverlaps;
	for (size_t i = 0; i < alignments.size(); i++)
	{
		for (size_t j = 0; j < alignments.size(); j++)
		{
			if (!alignmentPossible(alignmentSizes, alignmentNodeComparison, i, j, minSizeFraction)) continue;

			auto fw = getExactOverlaps(alignments[i].name(), alignmentNodeMovements[i], alignments[j].name(), alignmentNodeMovements[j], minMatchFraction, minSizeFraction, false);
			auto bw = getExactOverlaps(alignments[i].name(), alignmentNodeMovements[i], alignments[j].name(), reverse(alignmentNodeMovements[j]), minMatchFraction, minSizeFraction, true);
			validOverlaps.insert(validOverlaps.begin(), fw.begin(), fw.end());
			validOverlaps.insert(validOverlaps.begin(), bw.begin(), bw.end());
		}
	}

	std::ofstream outfile { argv[6], std::ios::out };
	for (size_t i = 0; i < reads.size(); i++)
	{
		outfile << "S\t" << reads[i].seq_id << "\t" << reads[i].sequence << std::endl;
	}

	for (size_t i = 0; i < validOverlaps.size(); i++)
	{
		outfile << "L\t" << validOverlaps[i].readname1 << "\t" << (validOverlaps[i].backward1 ? "-" : "+") << "\t" << validOverlaps[i].readname2 << "\t" << (validOverlaps[i].backward2 ? "-" : "+") << "\t" << validOverlaps[i].length1 << "M" << std::endl;
	}
}