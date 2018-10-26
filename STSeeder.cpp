#include "STSeeder.h"
#include "CommonUtils.h"

char lowercase(char c)
{
	switch(c)
	{
		case 'a':
		case 'c':
		case 't':
		case 'g':
			return c;
		case 'A': return 'a';
		case 'C': return 'c';
		case 'G': return 'g';
		case 'T': return 't';
		case '$': return '$';
	}
	return 0;
}

template <typename F>
void getMums(const std::string& sequence, const sdsl::cst_sada<>& tree, F mumCallback)
{
	size_t seqpos = 0;
	auto pos = tree.root();
	while (seqpos < sequence.size() && pos == tree.root())
	{
		pos = tree.child(tree.root(), lowercase(sequence[seqpos]));
		seqpos++;
	}
	size_t depth = 1;
	while (seqpos < sequence.size())
	{
		assert(depth <= tree.depth(pos));
		assert(depth > tree.depth(tree.parent(pos)) || pos == tree.root());
		auto c = lowercase(sequence[seqpos]);
		if (depth < tree.depth(pos))
		{
			if (tree.edge(pos, depth+1) != c)
			{
				if (tree.is_leaf(pos))
				{
					mumCallback(seqpos, tree.sn(pos), depth);
				}
				while (depth != tree.depth(pos) || tree.child(pos, c) == tree.root())
				{
					pos = tree.sl(pos);
					assert(depth > 0);
					depth--;
					while (depth <= tree.depth(tree.parent(pos)) && pos != tree.root()) pos = tree.parent(pos);
				}
				continue;
			}
			if (seqpos == sequence.size() - 1) break;
			depth++;
			seqpos++;
			continue;
		}
		auto nextpos = tree.child(pos, c);
		assert(depth == tree.depth(pos));
		if (nextpos == tree.root())
		{
			if (tree.is_leaf(pos))
			{
				mumCallback(seqpos, tree.sn(pos), depth);
			}
			while (nextpos == tree.root())
			{
				pos = tree.sl(pos);
				nextpos = tree.child(pos, c);
				assert(depth > 0);
				depth--;
				while (depth <= tree.depth(tree.parent(pos)) && pos != tree.root()) pos = tree.parent(pos);
				if (depth < tree.depth(pos)) break;
				if (pos == tree.root()) break;
			}
		}
		else
		{
			if (seqpos == sequence.size() - 1) break;
			pos = nextpos;
			seqpos++;
			depth++;
		}
	}
	if (tree.is_leaf(pos))
	{
		mumCallback(seqpos, tree.sn(pos), depth);
	}
}

STSeeder::STSeeder(const GfaGraph& graph)
{
	constructTree(graph);
}

STSeeder::STSeeder(const vg::Graph& graph)
{
	constructTree(graph);
}

void STSeeder::constructTree(const GfaGraph& graph)
{
	std::string seq;
	for (auto node : graph.nodes)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(node.first);
		nodeReverse.push_back(false);
		seq += node.second;
		seq += '$';
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(node.first);
		nodeReverse.push_back(true);
		seq += CommonUtils::ReverseComplement(node.second);
		seq += '$';
	}
	nodePositions.push_back(seq.size());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = lowercase(seq[i]);
	}
	sdsl::construct_im(tree, seq, 1);
}

void STSeeder::constructTree(const vg::Graph& graph)
{
	std::string seq;
	for (int i = 0; i < graph.node_size(); i++)
	{
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(graph.node(i).id());
		nodeReverse.push_back(false);
		seq += graph.node(i).sequence();
		seq += '$';
		nodePositions.push_back(seq.size());
		nodeIDs.push_back(graph.node(i).id());
		nodeReverse.push_back(true);
		seq += CommonUtils::ReverseComplement(graph.node(i).sequence());
		seq += '$';
	}
	nodePositions.push_back(seq.size());
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = lowercase(seq[i]);
	}
	sdsl::construct_im(tree, seq, 1);
}

size_t STSeeder::getNodeIndex(size_t indexPos) const
{
	auto next = std::upper_bound(nodePositions.begin(), nodePositions.end(), indexPos);
	assert(next != nodePositions.begin());
	size_t index = (next - nodePositions.begin()) - 1;
	assert(index < nodePositions.size()-1);
	return index;
}

std::vector<SeedHit> STSeeder::getMumSeeds(const std::string& sequence) const
{
	std::vector<SeedHit> result;
	addMumSeeds(result, sequence);
	std::sort(result.begin(), result.end(), [](const SeedHit& left, const SeedHit& right) { return left.matchLen > right.matchLen; });
	return result;
}

void STSeeder::addMumSeeds(std::vector<SeedHit>& result, const std::string& sequence) const
{
	getMums(sequence, tree, [this, &result](size_t seqPos, size_t indexPos, size_t matchLen)
	{
		auto index = getNodeIndex(indexPos);
		int nodeID = nodeIDs[index];
		size_t nodeOffset = indexPos - nodePositions[index];
		result.emplace_back(nodeID, nodeOffset, seqPos - matchLen, matchLen, nodeReverse[index]);
	});
}
