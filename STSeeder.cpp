#include <iostream>
#include "STSeeder.h"

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

STSeeder::STSeeder(const GfaGraph& graph)
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
		seq += node.second;
		seq += '$';
	}
	for (size_t i = 0; i < seq.size(); i++)
	{
		seq[i] = lowercase(seq[i]);
	}
	sdsl::construct_im(tree, seq, 1);
	std::cout << tree.size() << std::endl;
	std::cout << tree.nodes() << std::endl;
}

std::vector<SeedHit> STSeeder::getSeeds(const std::string& sequence, size_t minMatchSize) const
{
	std::string lowercaseSeq = sequence;
	for (size_t i = 0; i < sequence.size(); i++)
	{
		lowercaseSeq[i] = lowercase(sequence[i]);
	}
	std::vector<SeedHit> result;
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
		assert(lowercaseSeq.substr(seqpos - depth, depth) == sdsl::extract(tree, pos).substr(0, depth));
		auto c = lowercase(sequence[seqpos]);
		if (depth < tree.depth(pos))
		{
			if (tree.edge(pos, depth+1) != c)
			{
				if (tree.is_leaf(pos))
				{
					std::cout << "A: " << seqpos << " " << tree.sn(pos) << " " << depth << std::endl;
					std::cout << sequence.substr(seqpos - depth, depth) << std::endl;
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
				std::cout << "B: " << seqpos << " " << tree.sn(pos) << " " << depth << std::endl;
				std::cout << sequence.substr(seqpos - depth, depth) << std::endl;
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
			pos = nextpos;
			seqpos++;
			depth++;
		}
	}
	if (tree.is_leaf(pos))
	{
		std::cout << "C: " << seqpos << " " << tree.sn(pos) << " " << depth << std::endl;
		std::cout << sequence.substr(seqpos - depth, depth) << std::endl;
	}
	return result;
}
