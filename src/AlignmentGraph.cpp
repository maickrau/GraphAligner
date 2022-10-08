#include <iostream>
#include <limits>
#include <algorithm>
#include <queue>
#include "AlignmentGraph.h"
#include "CommonUtils.h"
#include "ThreadReadAssertion.h"

AlignmentGraph dummy;

AlignmentGraph AlignmentGraph::DummyGraph()
{
	return dummy;
}

AlignmentGraph::NodeEdgeIterator::NodeEdgeIterator(size_t implicitEdge, const size_t* startPointer, const size_t* endPointer) :
	implicitEdge(implicitEdge),
	startPointer(startPointer),
	endPointer(endPointer)
{
}

AlignmentGraph::EdgeIterator::EdgeIterator(size_t implicitEdge, const size_t* vecPointer) :
	implicitEdge(implicitEdge),
	vecPointer(vecPointer)
{
}

AlignmentGraph::EdgeIterator AlignmentGraph::NodeEdgeIterator::begin() const
{
	return AlignmentGraph::EdgeIterator { implicitEdge, startPointer };
}

AlignmentGraph::EdgeIterator AlignmentGraph::NodeEdgeIterator::end() const
{
	return AlignmentGraph::EdgeIterator { std::numeric_limits<size_t>::max(), endPointer };
}

size_t AlignmentGraph::NodeEdgeIterator::size() const
{
	return endPointer - startPointer + (implicitEdge != std::numeric_limits<size_t>::max() ? 1 : 0);
}

size_t AlignmentGraph::NodeEdgeIterator::operator[](size_t index) const
{
	if (implicitEdge != std::numeric_limits<size_t>::max())
	{
		if (index == 0) return implicitEdge;
		index -= 1;
	}
	return *(startPointer+index);
}

bool AlignmentGraph::EdgeIterator::operator==(const AlignmentGraph::EdgeIterator& other) const
{
	return implicitEdge == other.implicitEdge && vecPointer == other.vecPointer;
}

bool AlignmentGraph::EdgeIterator::operator!=(const AlignmentGraph::EdgeIterator& other) const
{
	return !(*this == other);
}

size_t AlignmentGraph::EdgeIterator::operator*() const
{
	if (implicitEdge != std::numeric_limits<size_t>::max()) return implicitEdge;
	return *vecPointer;
}

AlignmentGraph::EdgeIterator& AlignmentGraph::EdgeIterator::operator++()
{
	if (implicitEdge != std::numeric_limits<size_t>::max())
	{
		implicitEdge = std::numeric_limits<size_t>::max();
		return *this;
	}
	vecPointer += 1;
	return *this;
}

AlignmentGraph::EdgeIterator AlignmentGraph::EdgeIterator::operator++(int)
{
	auto result = *this;
	++*this;
	return result;
}

AlignmentGraph::AlignmentGraph() :
nodeLength(),
nodeLookup(),
nodeIDs(),
nodeSequences(),
bpSize(0),
ambiguousNodeSequences(),
firstAmbiguous(std::numeric_limits<size_t>::max()),
finalized(false)
{
}

void AlignmentGraph::ReserveNodes(size_t numNodes, size_t numSplitNodes)
{
	nodeSequences.reserve(numSplitNodes);
	ambiguousNodeSequences.reserve(numSplitNodes);
	nodeLookup.resize(numNodes);
	originalNodeSize.resize(numNodes, std::numeric_limits<size_t>::max());
	originalNodeName.resize(numNodes, "");
	nodeIDs.reserve(numSplitNodes);
	nodeLength.reserve(numSplitNodes);
	reverse.reserve(numSplitNodes);
	nodeOffset.reserve(numSplitNodes);
}

void AlignmentGraph::AddNode(int nodeId, const DNAString& sequence, const std::string& name, bool reverseNode, const std::vector<size_t>& breakpoints)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(nodeId < nodeLookup.size());
	assert(nodeId < originalNodeSize.size());
	assert(nodeId < originalNodeName.size());
	assert(nodeLookup[nodeId].size() == 0);
	assert(originalNodeSize[nodeId] == std::numeric_limits<size_t>::max());
	assert(originalNodeName[nodeId].size() == 0);
	originalNodeSize[nodeId] = sequence.size();
	originalNodeName[nodeId] = name;
	assert(breakpoints.size() >= 2);
	assert(breakpoints[0] == 0);
	assert(breakpoints.back() == sequence.size());
	for (size_t breakpoint = 1; breakpoint < breakpoints.size(); breakpoint++)
	{
		if (breakpoints[breakpoint] == breakpoints[breakpoint-1]) continue;
		assert(breakpoints[breakpoint] > breakpoints[breakpoint-1]);
		for (size_t offset = breakpoints[breakpoint-1]; offset < breakpoints[breakpoint]; offset += SPLIT_NODE_SIZE)
		{
			size_t size = SPLIT_NODE_SIZE;
			if (breakpoints[breakpoint] - offset < size) size = breakpoints[breakpoint] - offset;
			assert(size > 0);
			AddNode(nodeId, offset, sequence.substr(offset, size), reverseNode);
			if (offset > 0)
			{
				assert(nodeOffset.size() == nodeIDs.size());
				assert(nodeIDs[nodeIDs.size()-2] == nodeIDs[nodeIDs.size()-1]);
				assert(nodeOffset[nodeIDs.size()-2] + NodeLength(nodeIDs.size()-2) == nodeOffset[nodeIDs.size()-1]);
			}
		}
	}
}

void AlignmentGraph::AddNode(int nodeId, int offset, const std::string& sequence, bool reverseNode)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(sequence.size() <= SPLIT_NODE_SIZE);

	bpSize += sequence.size();
	nodeLookup[nodeId].push_back(nodeLength.size());
	nodeLength.push_back(sequence.size());
	nodeIDs.push_back(nodeId);
	reverse.push_back(reverseNode);
	nodeOffset.push_back(offset);
	NodeChunkSequence normalSeq;
	for (size_t i = 0; i < CHUNKS_IN_NODE; i++)
	{
		normalSeq[i] = 0;
	}
	AmbiguousChunkSequence ambiguousSeq;
	ambiguousSeq.A = 0;
	ambiguousSeq.C = 0;
	ambiguousSeq.G = 0;
	ambiguousSeq.T = 0;
	bool ambiguous = false;
	assert(sequence.size() <= sizeof(size_t)*8);
	for (size_t i = 0; i < sequence.size(); i++)
	{
		size_t chunk = i / BP_IN_CHUNK;
		assert(chunk < CHUNKS_IN_NODE);
		size_t offset = (i % BP_IN_CHUNK) * 2;
		switch(sequence[i])
		{
			case 'a':
			case 'A':
				ambiguousSeq.A |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)0) << offset;
				break;
			case 'c':
			case 'C':
				ambiguousSeq.C |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)1) << offset;
				break;
			case 'g':
			case 'G':
				ambiguousSeq.G |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)2) << offset;
				break;
			case 't':
			case 'T':
			case 'u':
			case 'U':
				ambiguousSeq.T |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)3) << offset;
				break;
			case 'r':
			case 'R':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'y':
			case 'Y':
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 's':
			case 'S':
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'w':
			case 'W':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'k':
			case 'K':
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'm':
			case 'M':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'b':
			case 'B':
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'd':
			case 'D':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'h':
			case 'H':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'v':
			case 'V':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'n':
			case 'N':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			default:
				assert(false);
		}
	}
	ambiguousNodes.push_back(ambiguous);
	if (ambiguous)
	{
		ambiguousNodeSequences.emplace_back(ambiguousSeq);
	}
	else
	{
		nodeSequences.emplace_back(normalSeq);
	}
	assert(nodeIDs.size() == nodeLength.size());
}

void AlignmentGraph::AddEdgeNodeId(int node_id_from, int node_id_to, size_t startOffset)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(node_id_from < nodeLookup.size());
	assert(node_id_to < nodeLookup.size());
	std::tuple<size_t, size_t, size_t> key { node_id_from, node_id_to, startOffset };
	tempConstructionOutEdges.insert(key);
}

void AlignmentGraph::buildEdges()
{
	hasImplicitOutEdge.resize(nodeLength.size(), false);
	std::vector<uint16_t> explicitInEdgeCount;
	std::vector<uint16_t> explicitOutEdgeCount;
	explicitInEdgeCount.resize(nodeLength.size(), 0);
	explicitOutEdgeCount.resize(nodeLength.size(), 0);
	size_t totalExplicitEdges = 0;
	for (size_t i = 0; i < nodeLookup.size(); i++)
	{
		for (size_t j = 1; j < nodeLookup[i].size(); j++)
		{
			size_t from = nodeLookup[i][j-1];
			size_t to = nodeLookup[i][j];
			if (to == from+1)
			{
				hasImplicitOutEdge[from] = true;
			}
			else
			{
				explicitOutEdgeCount[from] += 1;
				explicitInEdgeCount[to] += 1;
				totalExplicitEdges += 2;
				assert(explicitOutEdgeCount[from] < 1024);
				assert(explicitInEdgeCount[to] < 1024);
			}
		}
	}
	for (auto t : tempConstructionOutEdges)
	{
		size_t from = nodeLookup[std::get<0>(t)].back();
		assert(nodeOffset[from] + NodeLength(from) == originalNodeSize[std::get<0>(t)]);
		size_t to = GetUnitigNode(std::get<1>(t), std::get<2>(t));
		assert(nodeOffset[to] == std::get<2>(t));
		if (to == from+1)
		{
			hasImplicitOutEdge[from] = true;
		}
		else
		{
			explicitOutEdgeCount[from] += 1;
			explicitInEdgeCount[to] += 1;
			totalExplicitEdges += 2;
			assert(explicitOutEdgeCount[from] < 1024);
			assert(explicitInEdgeCount[to] < 1024);
		}
	}
	assert(totalExplicitEdges < 17592186044416ll); // 2^44, 44 bits for start index, 10 bits for in-edge count, 10 for out-edge count
	size_t start = 0;
	explicitEdges.resize(nodeLength.size());
	for (size_t i = 0; i < explicitEdges.size(); i++)
	{
		assert(explicitOutEdgeCount[i] < 1024);
		assert(explicitInEdgeCount[i] < 1024);
		explicitEdges[i] = (start << 20) + (explicitInEdgeCount[i] << 10) + explicitOutEdgeCount[i];
		start += explicitInEdgeCount[i] + explicitOutEdgeCount[i];
	}
	assert(start == totalExplicitEdges);
	edgeStorage.resize(totalExplicitEdges, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < nodeLookup.size(); i++)
	{
		for (size_t j = 1; j < nodeLookup[i].size(); j++)
		{
			size_t from = nodeLookup[i][j-1];
			size_t to = nodeLookup[i][j];
			if (to == from+1) continue;
			assert(explicitOutEdgeCount[from] >= 1);
			explicitOutEdgeCount[from] -= 1;
			size_t outIndex = (explicitEdges[from] >> 20) + explicitOutEdgeCount[from];
			assert(outIndex < edgeStorage.size());
			assert(edgeStorage[outIndex] == std::numeric_limits<size_t>::max());
			edgeStorage[outIndex] = to;
			assert(explicitInEdgeCount[to] >= 1);
			explicitInEdgeCount[to] -= 1;
			size_t inIndex = (explicitEdges[to] >> 20) + (explicitEdges[to] & 1023) + explicitInEdgeCount[to];
			assert(inIndex < edgeStorage.size());
			assert(edgeStorage[inIndex] == std::numeric_limits<size_t>::max());
			edgeStorage[inIndex] = from;
		}
	}
	for (auto t : tempConstructionOutEdges)
	{
		size_t from = nodeLookup[std::get<0>(t)].back();
		assert(nodeOffset[from] + NodeLength(from) == originalNodeSize[std::get<0>(t)]);
		size_t to = GetUnitigNode(std::get<1>(t), std::get<2>(t));
		assert(nodeOffset[to] == std::get<2>(t));
		if (to == from+1) continue;
		assert(explicitOutEdgeCount[from] >= 1);
		explicitOutEdgeCount[from] -= 1;
		size_t outIndex = (explicitEdges[from] >> 20) + explicitOutEdgeCount[from];
		assert(outIndex < edgeStorage.size());
		assert(edgeStorage[outIndex] == std::numeric_limits<size_t>::max());
		edgeStorage[outIndex] = to;
		assert(explicitInEdgeCount[to] >= 1);
		explicitInEdgeCount[to] -= 1;
		size_t inIndex = (explicitEdges[to] >> 20) + (explicitEdges[to] & 1023) + explicitInEdgeCount[to];
		assert(inIndex < edgeStorage.size());
		assert(edgeStorage[inIndex] == std::numeric_limits<size_t>::max());
		edgeStorage[inIndex] = from;
	}
	for (size_t i = 0; i < edgeStorage.size(); i++)
	{
		assert(edgeStorage[i] != std::numeric_limits<size_t>::max());
		assert(edgeStorage[i] < nodeLength.size());
	}
	std::set<std::tuple<size_t, size_t, size_t>> tmp;
	std::swap(tmp, tempConstructionOutEdges);
}

void AlignmentGraph::Finalize(int wordSize)
{
	assert(nodeSequences.size() + ambiguousNodeSequences.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	assert(nodeLookup.size() == originalNodeSize.size());
	assert(nodeLookup.size() == originalNodeName.size());
	for (size_t i = 0; i < nodeLookup.size(); i++)
	{
		assert(nodeLookup[i].size() > 0);
		assert(originalNodeName[i] != "");
		assert(originalNodeSize[i] != std::numeric_limits<size_t>::max());
		assert(originalNodeSize[i] != 0);
	}
	RenumberAmbiguousToEnd();
	buildEdges();
	ambiguousNodes.clear();
	findLinearizable();
	doComponentOrder();
	findChains();
	finalized = true;
	int specialNodes = 0;
	size_t edges = 0;
	assert(nodeSequences.size() + ambiguousNodeSequences.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	assert(nodeOffset.size() == nodeLength.size());
	nodeLength.shrink_to_fit();
	nodeIDs.shrink_to_fit();
	reverse.shrink_to_fit();
	nodeSequences.shrink_to_fit();
	ambiguousNodeSequences.shrink_to_fit();
#ifndef NDEBUG
	for (size_t i = 0; i < nodeLookup.size(); i++)
	{
		for (size_t j = 1; j < nodeLookup[i].size(); j++)
		{
			assert(nodeOffset[nodeLookup[i][j-1]] < nodeOffset[nodeLookup[i][j]]);
		}
	}
#endif
}

std::pair<bool, size_t> AlignmentGraph::findBubble(const size_t start, const std::vector<bool>& ignorableTip)
{
	std::vector<size_t> S;
	S.push_back(start);
	std::unordered_set<size_t> visited;
	std::unordered_set<size_t> seen;
	seen.insert(start);
	while (S.size() > 0)
	{
		const size_t v = S.back();
		S.pop_back();
		assert(seen.count(v) == 1);
		seen.erase(v);
		assert(visited.count(v) == 0);
		visited.insert(v);
		if (OutNeighbors(v).size() == 0) return std::make_pair(false, 0);
		for (const size_t u : OutNeighbors(v))
		{
			if (ignorableTip[u]) continue;
			if (u == v) continue;
			if (u == start) return std::make_pair(false, 0);
			assert(visited.count(u) == 0);
			seen.insert(u);
			bool hasNonvisitedParent = false;
			for (const size_t w : InNeighbors(u))
			{
				if (w == u) continue;
				if (!ignorableTip[w] && visited.count(w) == 0)
				{
					hasNonvisitedParent = true;
					break;
				}
			}
			if (!hasNonvisitedParent) S.push_back(u);
		}
		if (S.size() == 1 && seen.size() == 1 && seen.count(S[0]) == 1)
		{
			const size_t t = S.back();
			for (const size_t u : OutNeighbors(t))
			{
				if (u == start) return std::make_pair(false, 0);
			}
			return std::make_pair(true, t);
		}
	}
	return std::make_pair(false, 0);
}

size_t find(std::vector<size_t>& parent, size_t item)
{
	if (parent[item] == item) return item;
	std::vector<size_t> stack;
	stack.push_back(item);
	while (parent[stack.back()] != stack.back()) stack.push_back(parent[stack.back()]);
	for (size_t i : stack) parent[i] = stack.back();
	return stack.back();
}

void merge(std::vector<size_t>& parent, std::vector<size_t>& rank, size_t left, size_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	if (rank[left] < rank[right])
	{
		std::swap(left, right);
	}
	parent[right] = left;
	if (rank[left] == rank[right]) rank[left] += 1;
}

void AlignmentGraph::chainBubble(const size_t start, const std::vector<bool>& ignorableTip, std::vector<size_t>& rank)
{
	bool hasBubble;
	size_t bubbleEnd;
	std::tie(hasBubble, bubbleEnd) = findBubble(start, ignorableTip);
	if (!hasBubble) return;
	std::unordered_set<size_t> visited;
	std::vector<size_t> stack;
	stack.push_back(start);
	merge(chainNumber, rank, BigraphNodeID(start), BigraphNodeID(bubbleEnd));
	while (stack.size() > 0)
	{
		const size_t top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		if (ignorableTip[top]) continue;
		visited.insert(top);
		merge(chainNumber, rank, BigraphNodeID(start), BigraphNodeID(top));
		for (const auto neighbor : OutNeighbors(top))
		{
			if (visited.count(neighbor) == 1) continue;
			if (neighbor == bubbleEnd) continue;
			stack.push_back(neighbor);
		}
	}
}

void AlignmentGraph::fixChainApproxPos(const size_t start)
{
	assert(chainApproxPos[start] == std::numeric_limits<size_t>::max());
	assert(std::numeric_limits<size_t>::max() / SPLIT_NODE_SIZE > nodeLength.size());
	assert(std::numeric_limits<size_t>::max() > nodeLength.size() * SPLIT_NODE_SIZE * 2);
	std::vector<std::pair<size_t, size_t>> stack;
	size_t chain = chainNumber[start];
	stack.emplace_back(start, (nodeLength.size() + 5) * SPLIT_NODE_SIZE);
	while (stack.size() > 0)
	{
		size_t v;
		size_t dist;
		std::tie(v, dist) = stack.back();
		stack.pop_back();
		if (chainApproxPos[v] != std::numeric_limits<size_t>::max()) continue;
		chainApproxPos[v] = dist;
		for (const size_t diNode : nodeLookup[v])
		{
			for (const size_t u : OutNeighbors(diNode))
			{
				if (chainNumber[BigraphNodeID(u)] != chain) continue;
				if (chainApproxPos[BigraphNodeID(u)] != std::numeric_limits<size_t>::max()) continue;
				assert(nodeOffset[u] < dist);
				assert(dist - nodeOffset[u] < std::numeric_limits<size_t>::max() - nodeOffset[diNode] - NodeLength(diNode));
				assert(std::numeric_limits<size_t>::max() - NodeLength(u) > dist);
				stack.emplace_back(u, dist + nodeOffset[diNode] + NodeLength(diNode) - nodeOffset[u]);
			}
			for (const size_t u : InNeighbors(diNode))
			{
				if (chainNumber[BigraphNodeID(u)] != chain) continue;
				if (chainApproxPos[BigraphNodeID(u)] != std::numeric_limits<size_t>::max()) continue;
				assert(dist + nodeOffset[diNode] > nodeOffset[u] + NodeLength(u));
				stack.emplace_back(u, dist + nodeOffset[diNode] - (nodeOffset[u] + NodeLength(u)));
			}
		}
	}
}

phmap::flat_hash_map<size_t, std::unordered_set<size_t>> AlignmentGraph::chainTips(std::vector<size_t>& rank, std::vector<bool>& ignorableTip)
{
	assert(componentNumber.size() == NodeSize());
	std::vector<size_t> order;
	order.reserve(NodeSize());
	for (size_t i = 0; i < NodeSize(); i++)
	{
		order.push_back(i);
	}
	std::sort(order.begin(), order.end(), [this](size_t left, size_t right) { return componentNumber[left] < componentNumber[right]; });
	std::vector<bool> fwTipComponent;
	fwTipComponent.resize(componentNumber[order.back()]+1, true);
	for (size_t ind = order.size()-1; ind < order.size(); ind--)
	{
		size_t i = order[ind];
		if (!fwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : OutNeighbors(i))
		{
			assert(componentNumber[neighbor] >= componentNumber[i]);
			if (componentNumber[neighbor] == componentNumber[i])
			{
				fwTipComponent[componentNumber[i]] = false;
				break;
			}
			if (!fwTipComponent[componentNumber[neighbor]])
			{
				fwTipComponent[componentNumber[i]] = false;
				break;
			}
		}
	}
	for (size_t ind = order.size()-1; ind < order.size(); ind--)
	{
		size_t i = order[ind];
		if (!fwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : OutNeighbors(i))
		{
			assert(fwTipComponent[componentNumber[neighbor]]);
			merge(chainNumber, rank, BigraphNodeID(i), BigraphNodeID(neighbor));
		}
	}
	std::vector<bool> bwTipComponent;
	bwTipComponent.resize(componentNumber[order.back()]+1, true);
	for (size_t ind = 0; ind < order.size(); ind++)
	{
		size_t i = order[ind];
		if (!bwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : InNeighbors(i))
		{
			assert(componentNumber[neighbor] <= componentNumber[i]);
			if (componentNumber[neighbor] == componentNumber[i])
			{
				bwTipComponent[componentNumber[i]] = false;
				break;
			}
			if (!bwTipComponent[componentNumber[neighbor]])
			{
				bwTipComponent[componentNumber[i]] = false;
				break;
			}
		}
	}
	for (size_t ind = 0; ind < order.size(); ind++)
	{
		size_t i = order[ind];
		if (!bwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : InNeighbors(i))
		{
			assert(bwTipComponent[componentNumber[neighbor]]);
			merge(chainNumber, rank, BigraphNodeID(i), BigraphNodeID(neighbor));
		}
	}
	phmap::flat_hash_map<size_t, std::unordered_set<size_t>> result;
	for (size_t i = 0; i < NodeSize(); i++)
	{
		if (bwTipComponent[componentNumber[i]] || fwTipComponent[componentNumber[i]])
		{
			ignorableTip[i] = true;
		}
		if (bwTipComponent[componentNumber[i]])
		{
			for (auto neighbor : OutNeighbors(i))
			{
				if (BigraphNodeID(neighbor) == BigraphNodeID(i)) continue;
				result[BigraphNodeID(i)].insert(BigraphNodeID(neighbor));
			}
		}
		if (fwTipComponent[componentNumber[i]])
		{
			for (auto neighbor : InNeighbors(i))
			{
				if (BigraphNodeID(neighbor) == BigraphNodeID(i)) continue;
				result[BigraphNodeID(i)].insert(BigraphNodeID(neighbor));
			}
		}
	}
	return result;
}

void AlignmentGraph::chainCycles(std::vector<size_t>& rank, std::vector<bool>& ignorableTip)
{
	for (size_t i = 0; i < nodeLength.size(); i++)
	{
		size_t uniqueFwNeighbor = std::numeric_limits<size_t>::max();
		for (auto u : OutNeighbors(i))
		{
			assert(u < nodeLength.size());
			if (ignorableTip[u]) continue;
			if (u == i) continue;
			if (uniqueFwNeighbor == std::numeric_limits<size_t>::max())
			{
				uniqueFwNeighbor = u;
			}
			else
			{
				assert(u != uniqueFwNeighbor);
				uniqueFwNeighbor = std::numeric_limits<size_t>::max()-1;
			}
		}
		size_t uniqueBwNeighbor = std::numeric_limits<size_t>::max();
		for (auto u : InNeighbors(i))
		{
			assert(u < nodeLength.size());
			if (ignorableTip[u]) continue;
			if (u == i) continue;
			if (uniqueBwNeighbor == std::numeric_limits<size_t>::max())
			{
				uniqueBwNeighbor = u;
			}
			else if (u != uniqueBwNeighbor)
			{
				uniqueBwNeighbor = std::numeric_limits<size_t>::max()-1;
			}
		}
		if (uniqueFwNeighbor != uniqueBwNeighbor) continue;
		if (uniqueFwNeighbor == std::numeric_limits<size_t>::max()) continue;
		if (uniqueFwNeighbor == std::numeric_limits<size_t>::max()-1) continue;
		if (uniqueBwNeighbor == std::numeric_limits<size_t>::max()) continue;
		if (uniqueBwNeighbor == std::numeric_limits<size_t>::max()-1) continue;
		ignorableTip[i] = true;
		assert(uniqueBwNeighbor == uniqueFwNeighbor);
		merge(chainNumber, rank, BigraphNodeID(i), BigraphNodeID(uniqueFwNeighbor));
	}
}

void AlignmentGraph::findChains()
{
	chainNumber.resize(nodeLookup.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chainNumber.size(); i++)
	{
		chainNumber[i] = i;
	}
	std::vector<bool> ignorableTip;
	ignorableTip.resize(nodeLength.size(), false);
	std::vector<size_t> rank;
	rank.resize(chainNumber.size(), 0);
	auto tipChainers = chainTips(rank, ignorableTip);
	chainCycles(rank, ignorableTip);
	for (size_t i = 0; i < nodeLookup.size(); i++)
	{
		chainBubble(nodeLookup[i].back(), ignorableTip, rank);
	}
	for (auto& pair : tipChainers)
	{
		size_t uniqueNeighbor = std::numeric_limits<size_t>::max();
		for (auto n : pair.second)
		{
			if (uniqueNeighbor == std::numeric_limits<size_t>::max())
			{
				uniqueNeighbor = n;
			}
			if (uniqueNeighbor != n)
			{
				uniqueNeighbor = std::numeric_limits<size_t>::max()-1;
				break;
			}
		}
		assert(uniqueNeighbor != std::numeric_limits<size_t>::max());
		if (uniqueNeighbor == std::numeric_limits<size_t>::max()-1) continue;
		merge(chainNumber, rank, pair.first, *pair.second.begin());
	}
	{
		std::vector<size_t> tmp;
		std::vector<bool> tmp2;
		std::swap(rank, tmp);
		std::swap(ignorableTip, tmp2);
	}
	for (size_t i = 0; i < chainNumber.size(); i++)
	{
		find(chainNumber, i);
	}
	chainApproxPos.resize(nodeLookup.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chainApproxPos.size(); i++)
	{
		if (chainApproxPos[i] == std::numeric_limits<size_t>::max()) fixChainApproxPos(i);
	}
}

void AlignmentGraph::findLinearizable()
{
	linearizable.resize(nodeLength.size(), false);
	std::vector<bool> checked;
	checked.resize(nodeLength.size(), false);
	std::vector<size_t> stack;
	std::vector<bool> onStack;
	onStack.resize(nodeLength.size(), false);
	for (size_t node = 0; node < nodeLength.size(); node++)
	{
		if (checked[node]) continue;
		if (InNeighbors(node).size() != 1)
		{
			checked[node] = true;
			continue;
		}
		checked[node] = true;
		assert(InNeighbors(node).size() == 1);
		assert(stack.size() == 0);
		stack.push_back(node);
		onStack[node] = true;
		while (true)
		{
			assert(stack.size() <= nodeLength.size());
			if (InNeighbors(stack.back()).size() != 1)
			{
				for (size_t i = 0; i < stack.size()-1; i++)
				{
					assert(InNeighbors(stack[i]).size() == 1);
					checked[stack[i]] = true;
					linearizable[stack[i]] = true;
					onStack[stack[i]] = false;
				}
				linearizable[stack.back()] = false;
				checked[stack.back()] = true;
				onStack[stack.back()] = false;
				stack.clear();
				break;
			}
			assert(InNeighbors(stack.back()).size() == 1);
			if (checked[stack.back()])
			{
				for (size_t i = 0; i < stack.size()-1; i++)
				{
					assert(InNeighbors(stack[i]).size() == 1);
					checked[stack[i]] = true;
					linearizable[stack[i]] = true;
					onStack[stack[i]] = false;
				}
				linearizable[stack.back()] = false;
				checked[stack.back()] = true;
				onStack[stack.back()] = false;
				stack.clear();
				break;
			}
			assert(InNeighbors(stack.back()).size() == 1);
			auto neighbor = InNeighbors(stack.back())[0];
			if (neighbor == node)
			{
				for (size_t i = 0; i < stack.size(); i++)
				{
					checked[stack[i]] = true;
					linearizable[stack[i]] = false;
					onStack[stack[i]] = false;
				}
				stack.clear();
				break;
			}
			if (onStack[neighbor])
			{
				assert(neighbor != node);
				size_t i = stack.size();
				for (; i > 0; i--)
				{
					if (stack[i] == neighbor) break;
					checked[stack[i]] = true;
					linearizable[stack[i]] = false;
					onStack[stack[i]] = false;
				}
				for (size_t j = 0; j < i; j++)
				{
					checked[stack[j]] = true;
					linearizable[stack[j]] = true;
					onStack[stack[j]] = false;
				}
				stack.clear();
				break;
			}
			stack.push_back(neighbor);
			onStack[neighbor] = true;
		}
	}
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
size_t AlignmentGraph::NodeLength(size_t index) const
{
	return nodeLength[index];
}

char AlignmentGraph::NodeSequences(size_t node, size_t pos) const
{
	assert(pos < NodeLength(node));
	if (node < firstAmbiguous)
	{
		assert(node < nodeSequences.size());
		size_t chunk = pos / BP_IN_CHUNK;
		size_t offset = (pos % BP_IN_CHUNK) * 2;
		return "ACGT"[(nodeSequences[node][chunk] >> offset) & 3];
	}
	else
	{
		assert(node >= firstAmbiguous);
		assert(node - firstAmbiguous < ambiguousNodeSequences.size());
		assert(pos < sizeof(size_t) * 8);
		bool A = (ambiguousNodeSequences[node - firstAmbiguous].A >> pos) & 1;
		bool C = (ambiguousNodeSequences[node - firstAmbiguous].C >> pos) & 1;
		bool G = (ambiguousNodeSequences[node - firstAmbiguous].G >> pos) & 1;
		bool T = (ambiguousNodeSequences[node - firstAmbiguous].T >> pos) & 1;
		assert(A + C + G + T >= 1);
		assert(A + C + G + T <= 4);
		if ( A && !C && !G && !T) return 'A';
		if (!A &&  C && !G && !T) return 'C';
		if (!A && !C &&  G && !T) return 'G';
		if (!A && !C && !G &&  T) return 'T';
		if ( A && !C &&  G && !T) return 'R';
		if (!A &&  C && !G &&  T) return 'Y';
		if (!A &&  C &&  G && !T) return 'S';
		if ( A && !C && !G &&  T) return 'W';
		if (!A && !C &&  G &&  T) return 'K';
		if ( A &&  C && !G && !T) return 'M';
		if (!A &&  C &&  G &&  T) return 'B';
		if ( A && !C &&  G &&  T) return 'D';
		if ( A &&  C && !G &&  T) return 'H';
		if ( A &&  C &&  G && !T) return 'V';
		if ( A &&  C &&  G &&  T) return 'N';
		assert(false);
		return 'N';
	}
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
AlignmentGraph::NodeChunkSequence AlignmentGraph::NodeChunks(size_t index) const
{
	assert(index < nodeSequences.size());
	return nodeSequences[index];
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
AlignmentGraph::AmbiguousChunkSequence AlignmentGraph::AmbiguousNodeChunks(size_t index) const
{
	assert(index >= firstAmbiguous);
	assert(index - firstAmbiguous < ambiguousNodeSequences.size());
	return ambiguousNodeSequences[index - firstAmbiguous];
}

size_t AlignmentGraph::NodeSize() const
{
	return nodeLength.size();
}

class NodeWithDistance
{
public:
	NodeWithDistance(size_t node, bool start, size_t distance) : node(node), start(start), distance(distance) {};
	bool operator>(const NodeWithDistance& other) const
	{
		return distance > other.distance;
	}
	size_t node;
	bool start;
	size_t distance;
};

size_t AlignmentGraph::GetUnitigNode(int nodeId, size_t offset) const
{
	assert(nodeLookup[nodeId].size() > 0);
	//guess the index
	size_t index = nodeLookup[nodeId].size() * ((double)offset / (double)originalNodeSize.at(nodeId));
	if (index >= nodeLookup[nodeId].size()) index = nodeLookup[nodeId].size()-1;
	//go to the exact index
	while (index < nodeLookup[nodeId].size()-1 && (nodeOffset[nodeLookup[nodeId][index]] + NodeLength(nodeLookup[nodeId][index]) <= offset)) index++;
	while (index > 0 && (nodeOffset[nodeLookup[nodeId][index]] > offset)) index--;
	assert(index != nodeLookup[nodeId].size());
	size_t result = nodeLookup[nodeId][index];
	assert(BigraphNodeID(result) == nodeId);
	assert(nodeOffset[result] <= offset);
	assert(nodeOffset[result] + NodeLength(result) > offset);
	return result;
}

std::pair<int, size_t> AlignmentGraph::GetReversePosition(int nodeId, size_t offset) const
{
	assert(nodeId < nodeLookup.size());
	size_t originalSize = originalNodeSize[nodeId];
	assert(offset < originalSize);
	size_t newOffset = originalSize - offset - 1;
	assert(newOffset < originalSize);
	int reverseNodeId;
	if (nodeId % 2 == 0)
	{
		reverseNodeId = (nodeId / 2) * 2 + 1;
	}
	else
	{
		reverseNodeId = (nodeId / 2) * 2;
	}
	return std::make_pair(reverseNodeId, newOffset);
}

AlignmentGraph::MatrixPosition::MatrixPosition() :
	node(0),
	nodeOffset(0),
	seqPos(0)
{
}

AlignmentGraph::MatrixPosition::MatrixPosition(size_t node, size_t nodeOffset, size_t seqPos) :
	node(node),
	nodeOffset(nodeOffset),
	seqPos(seqPos)
{
}

bool AlignmentGraph::MatrixPosition::operator==(const AlignmentGraph::MatrixPosition& other) const
{
	return node == other.node && nodeOffset == other.nodeOffset && seqPos == other.seqPos;
}

bool AlignmentGraph::MatrixPosition::operator!=(const AlignmentGraph::MatrixPosition& other) const
{
	return !(*this == other);
}

std::string AlignmentGraph::BigraphNodeName(int nodeId) const
{
	assert(nodeId < originalNodeName.size());
	assert(originalNodeName[nodeId] != "");
	return originalNodeName[nodeId];
}

std::vector<size_t> renumber(const std::vector<size_t>& vec, const std::vector<size_t>& renumbering)
{
	std::vector<size_t> result;
	result.reserve(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		assert(vec[i] < renumbering.size());
		result.push_back(renumbering[vec[i]]);
	}
	return result;
}

template <typename T>
std::vector<T> reorder(const std::vector<T>& vec, const std::vector<size_t>& renumbering)
{
	assert(vec.size() == renumbering.size());
	std::vector<T> result;
	result.resize(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		result[renumbering[i]] = vec[i];
	}
	return result;
}

size_t AlignmentGraph::BigraphNodeSize(int nodeId) const
{
	assert(nodeId < originalNodeSize.size());
	assert(originalNodeSize[nodeId] != std::numeric_limits<size_t>::max());
	return originalNodeSize[nodeId];
}

void AlignmentGraph::RenumberAmbiguousToEnd()
{
	assert(nodeSequences.size() + ambiguousNodeSequences.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	assert(ambiguousNodes.size() == nodeLength.size());
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	std::vector<size_t> renumbering;
	renumbering.reserve(ambiguousNodes.size());
	size_t nonAmbiguousCount = 0;
	size_t ambiguousCount = 0;
	for (size_t i = 0; i < ambiguousNodes.size(); i++)
	{
		if (!ambiguousNodes[i])
		{
			renumbering.push_back(nonAmbiguousCount);
			nonAmbiguousCount++;
		}
		else
		{
			assert(ambiguousCount < ambiguousNodes.size());
			assert(ambiguousNodes.size()-1-ambiguousCount >= nonAmbiguousCount);
			renumbering.push_back(ambiguousNodes.size()-1-ambiguousCount);
			ambiguousCount++;
		}
	}
	assert(renumbering.size() == ambiguousNodes.size());
	assert(nonAmbiguousCount + ambiguousCount == ambiguousNodes.size());
	assert(ambiguousCount == ambiguousNodeSequences.size());
	assert(nonAmbiguousCount == nodeSequences.size());
	firstAmbiguous = nonAmbiguousCount;

	if (ambiguousCount == 0) return;

	//the ambiguous nodes were added in the reverse order, reverse the sequence containers too
	std::reverse(ambiguousNodeSequences.begin(), ambiguousNodeSequences.end());

	nodeLength = reorder(nodeLength, renumbering);
	nodeOffset = reorder(nodeOffset, renumbering);
	nodeIDs = reorder(nodeIDs, renumbering);
	reverse = reorder(reverse, renumbering);
	for (size_t i = 0; i < nodeLookup.size(); i++)
	{
		nodeLookup[i] = renumber(nodeLookup[i], renumbering);
	}

#ifndef NDEBUG
	for (size_t i = 0; i < nodeLookup.size(); i++)
	{
		size_t foundSize = 0;
		std::set<size_t> offsets;
		size_t lastOffset = 0;
		for (auto node : nodeLookup[i])
		{
			assert(offsets.count(nodeOffset[node]) == 0);
			assert(offsets.size() == 0 || nodeOffset[node] > lastOffset);
			lastOffset = nodeOffset[node];
			offsets.insert(nodeOffset[node]);
			assert(BigraphNodeID(node) == i);
			foundSize += NodeLength(node);
		}
		assert(foundSize == originalNodeSize[i]);
	}
#endif
}

void AlignmentGraph::doComponentOrder()
{
	std::vector<std::tuple<size_t, int, size_t>> callStack;
	size_t i = 0;
	std::vector<size_t> index;
	std::vector<size_t> lowlink;
	std::vector<bool> onStack;
	std::vector<size_t> stack;
	index.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	lowlink.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	onStack.resize(nodeLength.size(), false);
	size_t checknode = 0;
	size_t nextComponent = 0;
	componentNumber.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	while (true)
	{
		if (callStack.size() == 0)
		{
			while (checknode < nodeLength.size() && index[checknode] != std::numeric_limits<size_t>::max())
			{
				checknode++;
			}
			if (checknode == nodeLength.size()) break;
			callStack.emplace_back(checknode, 0, 0);
			checknode++;
		}
		auto top = callStack.back();
		const size_t v = std::get<0>(top);
		int state = std::get<1>(top);
		size_t w;
		size_t neighborI = std::get<2>(top);
		callStack.pop_back();
		switch(state)
		{
			case 0:
				assert(index[v] == std::numeric_limits<size_t>::max());
				assert(lowlink[v] == std::numeric_limits<size_t>::max());
				assert(!onStack[v]);
				index[v] = i;
				lowlink[v] = i;
				i += 1;
				stack.push_back(v);
				onStack[v] = true;
			startloop:
				[[fallthrough]];
			case 1:
				if (neighborI >= OutNeighbors(v).size()) goto endloop;
				assert(neighborI < OutNeighbors(v).size());
				w = OutNeighbors(v)[neighborI];
				if (index[w] == std::numeric_limits<size_t>::max())
				{
					assert(lowlink[w] == std::numeric_limits<size_t>::max());
					assert(!onStack[w]);
					callStack.emplace_back(v, 2, neighborI);
					callStack.emplace_back(w, 0, 0);
					continue;
				}
				else if (onStack[w])
				{
					lowlink[v] = std::min(lowlink[v], index[w]);
					neighborI += 1;
					goto startloop;
				}
				else
				{
					neighborI += 1;
					goto startloop;
				}
			case 2:
				assert(neighborI < OutNeighbors(v).size());
				w = OutNeighbors(v)[neighborI];
				assert(index[w] != std::numeric_limits<size_t>::max());
				assert(lowlink[w] != std::numeric_limits<size_t>::max());
				lowlink[v] = std::min(lowlink[v], lowlink[w]);
				neighborI++;
				goto startloop;
			endloop:
			case 3:
				if (lowlink[v] == index[v])
				{
					do
					{
						w = stack.back();
						stack.pop_back();
						onStack[w] = false;
						componentNumber[w] = nextComponent;
					} while (w != v);
					nextComponent++;
				}
		}
	}
	assert(stack.size() == 0);
	for (size_t i = 0; i < componentNumber.size(); i++)
	{
		assert(componentNumber[i] != std::numeric_limits<size_t>::max());
		assert(componentNumber[i] <= nextComponent-1);
		componentNumber[i] = nextComponent-1-componentNumber[i];
	}
#ifdef EXTRACORRECTNESSASSERTIONS
	for (size_t i = 0; i < nodeLength.size(); i++)
	{
		for (auto neighbor : OutNeighbors(i))
		{
			assert(componentNumber[neighbor] >= componentNumber[i]);
		}
	}
#endif
}

size_t AlignmentGraph::ComponentSize() const
{
	return componentNumber.size();
}

size_t AlignmentGraph::SizeInBP() const
{
	return bpSize;
}

size_t AlignmentGraph::BigraphNodeID(size_t directedNodeId) const
{
	return nodeIDs[directedNodeId];
}

AlignmentGraph::NodeEdgeIterator AlignmentGraph::OutNeighbors(size_t nodeId) const
{
	size_t explicitInfo = explicitEdges[nodeId];
	size_t implicitEdge = std::numeric_limits<size_t>::max();
	if (hasImplicitOutEdge[nodeId]) implicitEdge = nodeId+1;
	size_t startIndex = explicitInfo >> 20;
	size_t count = explicitInfo & 1023; 
	assert(startIndex + count <= edgeStorage.size());
	AlignmentGraph::NodeEdgeIterator result { implicitEdge, edgeStorage.data()+startIndex, edgeStorage.data()+startIndex+count };
	return result;
}

AlignmentGraph::NodeEdgeIterator AlignmentGraph::InNeighbors(size_t nodeId) const
{
	size_t explicitInfo = explicitEdges[nodeId];
	size_t implicitEdge = std::numeric_limits<size_t>::max();
	if (nodeId > 0 && hasImplicitOutEdge[nodeId-1]) implicitEdge = nodeId-1;
	size_t startIndex = (explicitInfo >> 20) + (explicitInfo & 1023);
	size_t count = (explicitInfo >> 10) & 1023;
	assert(startIndex + count <= edgeStorage.size());
	assert(implicitEdge == std::numeric_limits<size_t>::max() || implicitEdge < nodeLength.size());
	AlignmentGraph::NodeEdgeIterator result { implicitEdge, edgeStorage.data()+startIndex, edgeStorage.data()+startIndex+count };
	return result;
}

size_t AlignmentGraph::BigraphNodeCount() const
{
	return originalNodeName.size();
}
