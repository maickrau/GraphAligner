#include <iostream>
#include <limits>
#include <algorithm>
#include <queue>
#include "AlignmentGraph.h"
#include "CommonUtils.h"
#include "ThreadReadAssertion.h"

bool isNumber(const std::string name)
{
	assert(name.size() >= 1);
	if (name.size() == 0 && name == "0") return true;
	if (name[0] < '1') return false;
	if (name[0] > '9') return false;
	for (size_t i = 1; i < name.size(); i++)
	{
		if (name[i] < '0') return false;
		if (name[i] > '9') return false;
	}
	return true;
}

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
	bpSize(0),
	firstAmbiguous(std::numeric_limits<size_t>::max()),
	finalized(false),
	originalNodeName(),
	bigraphIntermediateList(),
	originalNodeSize(),
	chainNumber(),
	chainApproxPos(),
	reverse(),
	componentNumber(),
	lastDinodeLength(),
	firstDinodeOffset(),
	intermediateBigraphNodeIDs(),
	intermediateDinodesStart(),
	intermediateInEdges(),
	intermediateOutEdges(),
	firstOfIntermediates(),
	nodeSequences(),
	ambiguousNodeSequences(),
	allNodeNamesAreNumbers(true)
{
}

void AlignmentGraph::ReserveNodes(size_t numBigraphNodes, size_t numDigraphNodes)
{
	originalNodeName.reserve(numBigraphNodes);
	bigraphIntermediateList.reserve(numBigraphNodes);
	originalNodeSize.reserve(numBigraphNodes);
	chainNumber.reserve(numBigraphNodes);
	chainApproxPos.reserve(numBigraphNodes);
	reverse.reserve(numBigraphNodes);
	nodeSequences.reserve(numDigraphNodes);
}

void AlignmentGraph::AddNode(size_t bigraphNodeId, const DNAString& sequence, const std::string& name, bool reverseNode, const std::vector<size_t>& breakpoints)
{
	if (!isNumber(name)) allNodeNamesAreNumbers = false;
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(bigraphNodeId == BigraphNodeCount());
	originalNodeName.push_back(name);
	bigraphIntermediateList.emplace_back();
	originalNodeSize.push_back(sequence.size());
	reverse.push_back(reverseNode);
	assert(breakpoints.size() >= 2);
	assert(breakpoints[0] == 0);
	assert(breakpoints.back() == sequence.size());
	bpSize += sequence.size();
	size_t lastAdded = std::numeric_limits<size_t>::max();
	for (size_t breakpoint = 1; breakpoint < breakpoints.size(); breakpoint++)
	{
		if (breakpoints[breakpoint] == breakpoints[breakpoint-1]) continue;
		size_t thisNode = addIntermediateNodes(bigraphNodeId, sequence, breakpoints[breakpoint-1], breakpoints[breakpoint]);
		if (lastAdded != std::numeric_limits<size_t>::max())
		{
			intermediateOutEdges[lastAdded].push_back(thisNode);
			intermediateInEdges[thisNode].push_back(lastAdded);
		}
		lastAdded = thisNode;
	}
}

size_t AlignmentGraph::intermediateNodeCount() const
{
	return intermediateDinodesStart.size();
}

size_t AlignmentGraph::addIntermediateNodes(size_t bigraphNodeId, const DNAString& sequence, size_t start, size_t end)
{
	size_t thisNode = intermediateNodeCount();
	assert(end > start);
	bool currentlyAmbiguous = false;
	for (size_t offset = start; offset < end; offset += SPLIT_NODE_SIZE)
	{
		size_t size = SPLIT_NODE_SIZE;
		if (offset+size >= end) size = end - offset;
		assert(size > 0);
		bool ambiguous = false;
		std::string seq = sequence.substr(offset, size);
		for (size_t i = 0; i < size; i++)
		{
			switch(seq[i])
			{
				case 'a':
				case 'A':
				case 'c':
				case 'C':
				case 'g':
				case 'G':
				case 't':
				case 'T':
				case 'u':
				case 'U':
					break;
				default:
					ambiguous = true;
					break;
			}
		}
		if (offset > start && ambiguous != currentlyAmbiguous)
		{
			sequence.rewindIterators(size);
			size_t newNode = addIntermediateNodes(bigraphNodeId, sequence, offset, end);
			assert(std::find(intermediateOutEdges[thisNode].begin(), intermediateOutEdges[thisNode].end(), newNode) == intermediateOutEdges[thisNode].end());
			assert(std::find(intermediateInEdges[newNode].begin(), intermediateInEdges[newNode].end(), thisNode) == intermediateInEdges[newNode].end());
			intermediateOutEdges[thisNode].push_back(newNode);
			intermediateInEdges[newNode].push_back(thisNode);
			return thisNode;
		}
		if (offset == start)
		{
			currentlyAmbiguous = ambiguous;
			assert(intermediateNodeCount() == thisNode);
			lastDinodeLength.push_back(size);
			firstDinodeOffset.push_back(start);
			intermediateBigraphNodeIDs.push_back(bigraphNodeId);
			if (ambiguous)
			{
				intermediateDinodesStart.push_back(ambiguousNodeSequences.size() + std::numeric_limits<size_t>::max()/2);
			}
			else
			{
				intermediateDinodesStart.push_back(nodeSequences.size());
			}
			intermediateInEdges.emplace_back();
			intermediateOutEdges.emplace_back();
			bigraphIntermediateList[bigraphNodeId].push_back(thisNode);
			assert(intermediateNodeCount() == thisNode+1);
		}
		if (ambiguous)
		{
			AddAmbiguousDinode(seq);
		}
		else
		{
			AddNormalDinode(seq);
		}
		lastDinodeLength[thisNode] = seq.size();
	}
	return thisNode;
}

void AlignmentGraph::AddAmbiguousDinode(const std::string& sequence)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(sequence.size() <= SPLIT_NODE_SIZE);
	AmbiguousChunkSequence ambiguousSeq;
	ambiguousSeq.A = 0;
	ambiguousSeq.C = 0;
	ambiguousSeq.G = 0;
	ambiguousSeq.T = 0;
	assert(sequence.size() <= sizeof(size_t)*8);
	for (size_t i = 0; i < sequence.size(); i++)
	{
		size_t chunk = i / BP_IN_CHUNK;
		assert(chunk < CHUNKS_IN_NODE);
		switch(sequence[i])
		{
			case 'a':
			case 'A':
				ambiguousSeq.A |= ((size_t)1) << (i);
				break;
			case 'c':
			case 'C':
				ambiguousSeq.C |= ((size_t)1) << (i);
				break;
			case 'g':
			case 'G':
				ambiguousSeq.G |= ((size_t)1) << (i);
				break;
			case 't':
			case 'T':
			case 'u':
			case 'U':
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			case 'r':
			case 'R':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				break;
			case 'y':
			case 'Y':
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			case 's':
			case 'S':
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				break;
			case 'w':
			case 'W':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			case 'k':
			case 'K':
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			case 'm':
			case 'M':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				break;
			case 'b':
			case 'B':
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			case 'd':
			case 'D':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			case 'h':
			case 'H':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			case 'v':
			case 'V':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				break;
			case 'n':
			case 'N':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				break;
			default:
				assert(false);
		}
	}
	ambiguousNodeSequences.emplace_back(ambiguousSeq);
}

void AlignmentGraph::AddNormalDinode(const std::string& sequence)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(sequence.size() <= SPLIT_NODE_SIZE);
	NodeChunkSequence normalSeq;
	for (size_t i = 0; i < CHUNKS_IN_NODE; i++)
	{
		normalSeq[i] = 0;
	}
	for (size_t i = 0; i < sequence.size(); i++)
	{
		size_t chunk = i / BP_IN_CHUNK;
		assert(chunk < CHUNKS_IN_NODE);
		size_t offset = (i % BP_IN_CHUNK) * 2;
		switch(sequence[i])
		{
			case 'a':
			case 'A':
				normalSeq[chunk] |= ((size_t)0) << offset;
				break;
			case 'c':
			case 'C':
				normalSeq[chunk] |= ((size_t)1) << offset;
				break;
			case 'g':
			case 'G':
				normalSeq[chunk] |= ((size_t)2) << offset;
				break;
			case 't':
			case 'T':
			case 'u':
			case 'U':
				normalSeq[chunk] |= ((size_t)3) << offset;
				break;
			default:
				assert(false);
		}
	}
	nodeSequences.emplace_back(normalSeq);
}

void AlignmentGraph::AddEdgeNodeId(size_t bigraphNodeIdFrom, size_t bigraphNodeIdTo, size_t startOffset)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(bigraphNodeIdFrom < bigraphIntermediateList.size());
	assert(bigraphNodeIdTo < bigraphIntermediateList.size());
	assert(bigraphIntermediateList[bigraphNodeIdFrom].size() > 0);
	assert(bigraphIntermediateList[bigraphNodeIdTo].size() > 0);
	size_t intermediateFrom = bigraphIntermediateList[bigraphNodeIdFrom].back();
	size_t intermediateTo = std::numeric_limits<size_t>::max();
	for (auto node : bigraphIntermediateList[bigraphNodeIdTo])
	{
		if (firstDinodeOffset[node] == startOffset)
		{
			intermediateTo = node;
			break;
		}
	}
	assert(intermediateTo != std::numeric_limits<size_t>::max());
	if (std::find(intermediateOutEdges[intermediateFrom].begin(), intermediateOutEdges[intermediateFrom].end(), intermediateTo) == intermediateOutEdges[intermediateFrom].end()) intermediateOutEdges[intermediateFrom].push_back(intermediateTo);
	if (std::find(intermediateInEdges[intermediateTo].begin(), intermediateInEdges[intermediateTo].end(), intermediateFrom) == intermediateInEdges[intermediateTo].end()) intermediateInEdges[intermediateTo].push_back(intermediateFrom);
}

void AlignmentGraph::Finalize(int wordSize)
{
	assert(originalNodeName.size() == BigraphNodeCount());
	assert(bigraphIntermediateList.size() == BigraphNodeCount());
	assert(originalNodeSize.size() == BigraphNodeCount());
	assert(reverse.size() == BigraphNodeCount());
	assert(chainNumber.size() == 0);
	assert(chainApproxPos.size() == 0);
	assert(intermediateOutEdges.size() == intermediateInEdges.size());
	for (size_t i = 0; i < intermediateOutEdges.size(); i++)
	{
		for (auto neighbor : intermediateOutEdges[i])
		{
			assert(std::find(intermediateInEdges[neighbor].begin(), intermediateInEdges[neighbor].end(), i) != intermediateInEdges[neighbor].end());
		}
	}
	for (size_t i = 0; i < bigraphIntermediateList.size(); i++)
	{
		assert(bigraphIntermediateList[i].size() > 0);
		assert(originalNodeName[i] != "");
		assert(originalNodeSize[i] != std::numeric_limits<size_t>::max());
		assert(originalNodeSize[i] != 0);
	}
	RenumberAmbiguousToEnd();
	makeDinodeIntermediateMapping();
	doComponentOrder();
	findChains();
	sparsenComponentNumbers();
	replaceIntermediateEdgesWithDinodes();
	finalized = true;

	assert(chainNumber.size() == BigraphNodeCount());
	assert(chainApproxPos.size() == BigraphNodeCount());
	assert(nodeSequences.size() + ambiguousNodeSequences.size() == NodeSize());
	assert(firstOfIntermediates.size() == NodeSize()+1);

	assert(originalNodeName.size() == BigraphNodeCount());
	assert(bigraphIntermediateList.size() == BigraphNodeCount());
	assert(originalNodeSize.size() == BigraphNodeCount());
	assert(chainNumber.size() == BigraphNodeCount());
	assert(chainApproxPos.size() == BigraphNodeCount());
	assert(reverse.size() == BigraphNodeCount());
	assert(componentNumber.size() == lastDinodeLength.size());
	assert(componentNumber.size() == firstDinodeOffset.size());
	assert(componentNumber.size() == intermediateBigraphNodeIDs.size());
	assert(componentNumber.size() == intermediateDinodesStart.size());
	assert(componentNumber.size() == intermediateInEdges.size());
	assert(componentNumber.size() == intermediateOutEdges.size());
	for (size_t i = 0; i < BigraphNodeCount(); i++)
	{
		size_t intermediateLenSum = 0;
		for (auto interNode : bigraphIntermediateList[i])
		{
			assert(intermediateBigraphNodeIDs[interNode] == i);
			intermediateLenSum += intermediateNodeLength(interNode);
		}
		for (size_t j = 1; j < bigraphIntermediateList[i].size(); j++)
		{
			assert(firstDinodeOffset[bigraphIntermediateList[i][j]] > firstDinodeOffset[bigraphIntermediateList[i][j-1]]);
		}
		assert(intermediateLenSum == BigraphNodeSize(i));
	}
	for (size_t i = 0; i < intermediateBigraphNodeIDs.size(); i++)
	{
		assert(i == 0 || intermediateDinodesStart[i] > intermediateDinodesStart[i-1]);
		size_t bigraphNodeId = intermediateBigraphNodeIDs[i];
		assert(std::find(bigraphIntermediateList[bigraphNodeId].begin(), bigraphIntermediateList[bigraphNodeId].end(), i) != bigraphIntermediateList[bigraphNodeId].end());
		for (auto neighborDinode : intermediateOutEdges[i])
		{
			size_t neighbor = digraphToIntermediate(neighborDinode);
			assert(componentNumber[neighbor] >= componentNumber[i]);
			size_t last = intermediateLastDinode(i);
			assert(std::find(intermediateInEdges[neighbor].begin(), intermediateInEdges[neighbor].end(), last) != intermediateInEdges[neighbor].end());
		}
		for (auto neighborDinode : intermediateInEdges[i])
		{
			size_t neighbor = digraphToIntermediate(neighborDinode);
			assert(componentNumber[neighbor] <= componentNumber[i]);
			size_t first = intermediateDinodesStart[i];
			assert(std::find(intermediateOutEdges[neighbor].begin(), intermediateOutEdges[neighbor].end(), first) != intermediateOutEdges[neighbor].end());
		}
		assert(lastDinodeLength[i] == NodeLength(intermediateLastDinode(i)));
	}
	for (size_t i = 0; i < NodeSize(); i++)
	{
		for (auto neighbor : OutNeighbors(i))
		{
			assert(neighbor < NodeSize());
			auto other = InNeighbors(neighbor);
			assert(std::find(other.begin(), other.end(), i) != other.end());
		}
		for (auto neighbor : InNeighbors(i))
		{
			assert(neighbor < NodeSize());
			auto other = OutNeighbors(neighbor);
			assert(std::find(other.begin(), other.end(), i) != other.end());
		}
	}
}

void AlignmentGraph::makeDinodeIntermediateMapping()
{
	firstOfIntermediates.resize(nodeSequences.size() + ambiguousNodeSequences.size() + 1);
	for (size_t i = 0; i < intermediateDinodesStart.size(); i++)
	{
		assert(!firstOfIntermediates.get(intermediateDinodesStart[i]));
		firstOfIntermediates.set(intermediateDinodesStart[i], true);
	}
	firstOfIntermediates.set(nodeSequences.size() + ambiguousNodeSequences.size(), true);
	firstOfIntermediates.buildRanks();
}

void AlignmentGraph::sparsenComponentNumbers()
{
	partOfStronglyConnectedComponent.resize(componentNumber.size(), false);
	for (size_t i = 0; i < intermediateOutEdges.size(); i++)
	{
		for (auto edge : intermediateOutEdges[i])
		{
			assert(componentNumber[edge] >= componentNumber[i]);
			if (componentNumber[edge] == componentNumber[i])
			{
				partOfStronglyConnectedComponent[i] = true;
				partOfStronglyConnectedComponent[edge] = true;
			}
		}
		for (auto edge : intermediateInEdges[i])
		{
			assert(componentNumber[edge] <= componentNumber[i]);
			if (componentNumber[edge] == componentNumber[i])
			{
				partOfStronglyConnectedComponent[i] = true;
				partOfStronglyConnectedComponent[edge] = true;
			}
		}
	}
	std::vector<size_t> order;
	order.reserve(componentNumber.size());
	for (size_t i = 0; i < componentNumber.size(); i++)
	{
		order.push_back(i);
	}
	std::sort(order.begin(), order.end(), [this](size_t left, size_t right) { return componentNumber[left] < componentNumber[right]; });
	size_t extraAdd = 0;
	size_t nextAdd = 0;
	size_t lastNumber = std::numeric_limits<size_t>::max();
	size_t lastRawNumber = std::numeric_limits<size_t>::max();
	for (size_t i : order)
	{
		assert(componentNumber[i] >= lastRawNumber || lastRawNumber == std::numeric_limits<size_t>::max());
		if (componentNumber[i] != lastRawNumber)
		{
			extraAdd += nextAdd;
			nextAdd = 0;
		}
		else
		{
			assert(partOfStronglyConnectedComponent[i]);
		}
		lastRawNumber = componentNumber[i];
		componentNumber[i] += extraAdd;
		assert(componentNumber[i] >= lastNumber || lastNumber == std::numeric_limits<size_t>::max());
		nextAdd += intermediateDinodesCount(i);
		lastNumber = componentNumber[i];
	}
}

void AlignmentGraph::replaceIntermediateEdgesWithDinodes()
{
	for (size_t i = 0; i < intermediateOutEdges.size(); i++)
	{
		for (auto neighbor : intermediateOutEdges[i])
		{
			assert(std::find(intermediateInEdges[neighbor].begin(), intermediateInEdges[neighbor].end(), i) != intermediateInEdges[neighbor].end());
		}
	}
	for (size_t i = 0; i < intermediateInEdges.size(); i++)
	{
		for (auto neighbor : intermediateInEdges[i])
		{
			assert(std::find(intermediateOutEdges[neighbor].begin(), intermediateOutEdges[neighbor].end(), i) != intermediateOutEdges[neighbor].end());
		}
	}
	for (size_t i = 0; i < intermediateOutEdges.size(); i++)
	{
		for (size_t j = 0; j < intermediateOutEdges[i].size(); j++)
		{
			intermediateOutEdges[i][j] = intermediateDinodesStart[intermediateOutEdges[i][j]];
		}
	}
	for (size_t i = 0; i < intermediateInEdges.size(); i++)
	{
		for (size_t j = 0; j < intermediateInEdges[i].size(); j++)
		{
			intermediateInEdges[i][j] = intermediateLastDinode(intermediateInEdges[i][j]);
		}
	}
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
		if (intermediateOutEdges[v].size() == 0) return std::make_pair(false, 0);
		for (const size_t u : intermediateOutEdges[v])
		{
			if (ignorableTip[u]) continue;
			if (u == v) continue;
			if (u == start) return std::make_pair(false, 0);
			assert(visited.count(u) == 0);
			seen.insert(u);
			bool hasNonvisitedParent = false;
			for (const size_t w : intermediateInEdges[u])
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
			for (const size_t u : intermediateOutEdges[t])
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
	merge(chainNumber, rank, intermediateBigraphNodeIDs[start], intermediateBigraphNodeIDs[bubbleEnd]);
	while (stack.size() > 0)
	{
		const size_t top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		if (ignorableTip[top]) continue;
		visited.insert(top);
		merge(chainNumber, rank, intermediateBigraphNodeIDs[start], intermediateBigraphNodeIDs[top]);
		for (const auto neighbor : intermediateOutEdges[top])
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
	assert(std::numeric_limits<size_t>::max() / SPLIT_NODE_SIZE > NodeSize());
	assert(std::numeric_limits<size_t>::max() > NodeSize() * SPLIT_NODE_SIZE * 2);
	std::vector<std::pair<size_t, size_t>> stack;
	size_t chain = chainNumber[start];
	stack.emplace_back(start, std::numeric_limits<size_t>::max()/2);
	while (stack.size() > 0)
	{
		size_t v;
		size_t dist;
		std::tie(v, dist) = stack.back();
		stack.pop_back();
		assert(v < chainApproxPos.size());
		if (chainApproxPos[v] != std::numeric_limits<size_t>::max()) continue;
		chainApproxPos[v] = dist;
		for (const size_t interNode : bigraphIntermediateList[v])
		{
			for (const size_t u : intermediateOutEdges[interNode])
			{
				if (chainNumber[intermediateBigraphNodeIDs[u]] != chain) continue;
				if (chainApproxPos[intermediateBigraphNodeIDs[u]] != std::numeric_limits<size_t>::max()) continue;
				assert(firstDinodeOffset[u] < dist);
				assert(dist - firstDinodeOffset[u] < std::numeric_limits<size_t>::max() - firstDinodeOffset[interNode] - intermediateNodeLength(interNode));
				assert(std::numeric_limits<size_t>::max() - intermediateNodeLength(u) > dist);
				stack.emplace_back(intermediateBigraphNodeIDs[u], dist + firstDinodeOffset[interNode] + intermediateNodeLength(interNode) - firstDinodeOffset[u]);
			}
			for (const size_t u : intermediateInEdges[interNode])
			{
				if (chainNumber[intermediateBigraphNodeIDs[u]] != chain) continue;
				if (chainApproxPos[intermediateBigraphNodeIDs[u]] != std::numeric_limits<size_t>::max()) continue;
				assert(dist + firstDinodeOffset[interNode] > firstDinodeOffset[u] + intermediateNodeLength(u));
				stack.emplace_back(intermediateBigraphNodeIDs[u], dist + firstDinodeOffset[interNode] - (firstDinodeOffset[u] + intermediateNodeLength(u)));
			}
		}
	}
}

size_t AlignmentGraph::intermediateNodeLength(size_t intermediate) const
{
	size_t bigraphNodeId = intermediateBigraphNodeIDs[intermediate];
	for (size_t i = 0; i < bigraphIntermediateList[bigraphNodeId].size(); i++)
	{
		if (bigraphIntermediateList[bigraphNodeId][i] != intermediate) continue;
		if (i+1 < bigraphIntermediateList[bigraphNodeId].size())
		{
			return firstDinodeOffset[bigraphIntermediateList[bigraphNodeId][i+1]] - firstDinodeOffset[bigraphIntermediateList[bigraphNodeId][i]];
		}
		return BigraphNodeSize(bigraphNodeId) - firstDinodeOffset[bigraphIntermediateList[bigraphNodeId][i]];
	}
	assert(false);
	return 0;
}

phmap::flat_hash_map<size_t, std::unordered_set<size_t>> AlignmentGraph::chainTips(std::vector<size_t>& rank, std::vector<bool>& ignorableTip)
{
	assert(componentNumber.size() == intermediateBigraphNodeIDs.size());
	std::vector<size_t> order;
	order.reserve(intermediateBigraphNodeIDs.size());
	for (size_t i = 0; i < intermediateBigraphNodeIDs.size(); i++)
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
		for (auto neighbor : intermediateOutEdges[i])
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
		for (auto neighbor : intermediateOutEdges[i])
		{
			assert(fwTipComponent[componentNumber[neighbor]]);
			merge(chainNumber, rank, intermediateBigraphNodeIDs[i], intermediateBigraphNodeIDs[neighbor]);
		}
	}
	std::vector<bool> bwTipComponent;
	bwTipComponent.resize(componentNumber[order.back()]+1, true);
	for (size_t ind = 0; ind < order.size(); ind++)
	{
		size_t i = order[ind];
		if (!bwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : intermediateInEdges[i])
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
		for (auto neighbor : intermediateInEdges[i])
		{
			assert(bwTipComponent[componentNumber[neighbor]]);
			merge(chainNumber, rank, intermediateBigraphNodeIDs[i], intermediateBigraphNodeIDs[neighbor]);
		}
	}
	phmap::flat_hash_map<size_t, std::unordered_set<size_t>> result;
	for (size_t i = 0; i < intermediateBigraphNodeIDs.size(); i++)
	{
		if (bwTipComponent[componentNumber[i]] || fwTipComponent[componentNumber[i]])
		{
			ignorableTip[i] = true;
		}
		if (bwTipComponent[componentNumber[i]])
		{
			for (auto neighbor : intermediateOutEdges[i])
			{
				if (intermediateBigraphNodeIDs[neighbor] == intermediateBigraphNodeIDs[i]) continue;
				result[intermediateBigraphNodeIDs[i]].insert(intermediateBigraphNodeIDs[neighbor]);
			}
		}
		if (fwTipComponent[componentNumber[i]])
		{
			for (auto neighbor : intermediateInEdges[i])
			{
				if (intermediateBigraphNodeIDs[neighbor] == intermediateBigraphNodeIDs[i]) continue;
				result[intermediateBigraphNodeIDs[i]].insert(intermediateBigraphNodeIDs[neighbor]);
			}
		}
	}
	return result;
}

void AlignmentGraph::chainCycles(std::vector<size_t>& rank, std::vector<bool>& ignorableTip)
{
	for (size_t i = 0; i < intermediateBigraphNodeIDs.size(); i++)
	{
		size_t uniqueFwNeighbor = std::numeric_limits<size_t>::max();
		for (auto u : intermediateOutEdges[i])
		{
			assert(u < intermediateBigraphNodeIDs.size());
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
		for (auto u : intermediateInEdges[i])
		{
			assert(u < intermediateBigraphNodeIDs.size());
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
		merge(chainNumber, rank, intermediateBigraphNodeIDs[i], intermediateBigraphNodeIDs[uniqueFwNeighbor]);
	}
}

void AlignmentGraph::findChains()
{
	chainNumber.resize(BigraphNodeCount(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chainNumber.size(); i++)
	{
		chainNumber[i] = i;
	}
	std::vector<bool> ignorableTip;
	ignorableTip.resize(intermediateInEdges.size(), false);
	std::vector<size_t> rank;
	rank.resize(chainNumber.size(), 0);
	auto tipChainers = chainTips(rank, ignorableTip);
	chainCycles(rank, ignorableTip);
	for (size_t i = 0; i < bigraphIntermediateList.size(); i++)
	{
		chainBubble(bigraphIntermediateList[i].back(), ignorableTip, rank);
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
	for (size_t i = 0; i < chainNumber.size(); i++)
	{
		find(chainNumber, i);
	}
	chainApproxPos.resize(BigraphNodeCount(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chainApproxPos.size(); i++)
	{
		if (chainApproxPos[i] == std::numeric_limits<size_t>::max()) fixChainApproxPos(i);
	}
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
size_t AlignmentGraph::NodeLength(size_t digraphNodeId) const
{
	if (firstOfIntermediates.get(digraphNodeId+1))
	{
		return lastDinodeLength[digraphToIntermediate(digraphNodeId)];
	}
	return SPLIT_NODE_SIZE;
}

char AlignmentGraph::NodeSequences(size_t digraphNodeId, size_t pos) const
{
	assert(pos < NodeLength(digraphNodeId));
	if (digraphNodeId < firstAmbiguous)
	{
		assert(digraphNodeId < nodeSequences.size());
		size_t chunk = pos / BP_IN_CHUNK;
		size_t offset = (pos % BP_IN_CHUNK) * 2;
		return "ACGT"[(nodeSequences[digraphNodeId][chunk] >> offset) & 3];
	}
	else
	{
		assert(digraphNodeId >= firstAmbiguous);
		assert(digraphNodeId - firstAmbiguous < ambiguousNodeSequences.size());
		assert(pos < sizeof(size_t) * 8);
		bool A = (ambiguousNodeSequences[digraphNodeId - firstAmbiguous].A >> pos) & 1;
		bool C = (ambiguousNodeSequences[digraphNodeId - firstAmbiguous].C >> pos) & 1;
		bool G = (ambiguousNodeSequences[digraphNodeId - firstAmbiguous].G >> pos) & 1;
		bool T = (ambiguousNodeSequences[digraphNodeId - firstAmbiguous].T >> pos) & 1;
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
	return nodeSequences.size() + ambiguousNodeSequences.size();
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

size_t AlignmentGraph::GetDigraphNode(size_t bigraphNodeId, size_t offset) const
{
	assert(offset < originalNodeSize[bigraphNodeId]);
	assert(bigraphIntermediateList[bigraphNodeId].size() > 0);
	size_t intermediateNode = bigraphIntermediateList[bigraphNodeId].back();
	for (size_t i = 1; i < bigraphIntermediateList[bigraphNodeId].size(); i++)
	{
		if (firstDinodeOffset[bigraphIntermediateList[bigraphNodeId][i]] > offset)
		{
			intermediateNode = bigraphIntermediateList[bigraphNodeId][i-1];
			break;
		}
	}
	assert(offset >= firstDinodeOffset[intermediateNode]);
	size_t extraOffset = offset - firstDinodeOffset[intermediateNode];
	return intermediateDinodesStart[intermediateNode] + extraOffset / SPLIT_NODE_SIZE;
}

std::pair<size_t, size_t> AlignmentGraph::GetReverseDigraphPosition(size_t digraphNodeId, size_t offset) const
{
	size_t bigraphNodeId = BigraphNodeID(digraphNodeId);
	size_t bigraphOffset = offset + NodeOffset(digraphNodeId);
	auto reverseBigraphPos = GetReversePosition(bigraphNodeId, bigraphOffset);
	size_t revdigraphNodeId = GetDigraphNode(reverseBigraphPos.first, reverseBigraphPos.second);
	size_t revoffset = reverseBigraphPos.second - NodeOffset(revdigraphNodeId);
	assert(revoffset < NodeLength(revdigraphNodeId));
	return std::make_pair(revdigraphNodeId, revoffset);
}

std::pair<size_t, size_t> AlignmentGraph::GetReversePosition(size_t bigraphNodeId, size_t offset) const
{
	assert(bigraphNodeId < BigraphNodeCount());
	size_t originalSize = originalNodeSize[bigraphNodeId];
	assert(offset < originalSize);
	size_t newOffset = originalSize - offset - 1;
	assert(newOffset < originalSize);
	size_t reverseNodeId;
	if (bigraphNodeId % 2 == 0)
	{
		reverseNodeId = (bigraphNodeId / 2) * 2 + 1;
	}
	else
	{
		reverseNodeId = (bigraphNodeId / 2) * 2;
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

std::string AlignmentGraph::BigraphNodeName(size_t bigraphNodeId) const
{
	assert(bigraphNodeId < originalNodeName.size());
	assert(originalNodeName[bigraphNodeId] != "");
	return originalNodeName[bigraphNodeId];
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

size_t AlignmentGraph::BigraphNodeSize(size_t bigraphNodeId) const
{
	assert(bigraphNodeId < originalNodeSize.size());
	assert(originalNodeSize[bigraphNodeId] != std::numeric_limits<size_t>::max());
	return originalNodeSize[bigraphNodeId];
}

void AlignmentGraph::RenumberAmbiguousToEnd()
{
	assert(!finalized);
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	std::vector<size_t> renumbering;
	renumbering.resize(intermediateDinodesStart.size(), std::numeric_limits<size_t>::max());
	size_t normalCount = 0;
	size_t ambiguousCount = 0;
	for (size_t i = 0; i < intermediateDinodesStart.size(); i++)
	{
		if (intermediateDinodesStart[i] >= std::numeric_limits<size_t>::max()/2)
		{
			ambiguousCount += 1;
		}
		else
		{
			normalCount += 1;
		}
	}
	assert(ambiguousCount + normalCount == renumbering.size());
	size_t normalDone = 0;
	size_t ambiguousDone = 0;
	for (size_t i = 0; i < intermediateDinodesStart.size(); i++)
	{
		if (intermediateDinodesStart[i] >= std::numeric_limits<size_t>::max()/2)
		{
			renumbering[i] = normalCount + ambiguousDone;
			ambiguousDone += 1;
		}
		else
		{
			renumbering[i] = normalDone;
			normalDone += 1;
		}
		assert(renumbering[i] < renumbering.size());
	}
	for (size_t i = 0; i < renumbering.size(); i++)
	{
		assert(renumbering[i] != std::numeric_limits<size_t>::max());
	}
	assert(normalDone == normalCount);
	assert(ambiguousDone == ambiguousCount);
	assert(normalDone + ambiguousDone == intermediateDinodesStart.size());
	firstAmbiguous = nodeSequences.size();
	if (ambiguousCount == 0) return;

	for (size_t i = 0; i < bigraphIntermediateList.size(); i++)
	{
		bigraphIntermediateList[i] = renumber(bigraphIntermediateList[i], renumbering);
	}
	lastDinodeLength = reorder(lastDinodeLength, renumbering);
	firstDinodeOffset = reorder(firstDinodeOffset, renumbering);
	intermediateBigraphNodeIDs = reorder(intermediateBigraphNodeIDs, renumbering);
	intermediateDinodesStart = reorder(intermediateDinodesStart, renumbering);
	intermediateInEdges = reorder(intermediateInEdges, renumbering);
	intermediateOutEdges = reorder(intermediateOutEdges, renumbering);
	for (size_t i = 0; i < intermediateInEdges.size(); i++)
	{
		intermediateInEdges[i] = renumber(intermediateInEdges[i], renumbering);
		intermediateOutEdges[i] = renumber(intermediateOutEdges[i], renumbering);
	}
	for (size_t i = 0; i < intermediateDinodesStart.size(); i++)
	{
		assert(i == 0 || intermediateDinodesStart[i] > intermediateDinodesStart[i-1]);
		if (intermediateDinodesStart[i] >= std::numeric_limits<size_t>::max()/2)
		{
			intermediateDinodesStart[i] -= std::numeric_limits<size_t>::max()/2;
			intermediateDinodesStart[i] += nodeSequences.size();
		}
		assert(i == 0 || intermediateDinodesStart[i] > intermediateDinodesStart[i-1]);
	}
}

void AlignmentGraph::doComponentOrder()
{
	for (size_t i = 0; i < intermediateOutEdges.size(); i++)
	{
		for (auto neighbor : intermediateOutEdges[i])
		{
			assert(std::find(intermediateInEdges[neighbor].begin(), intermediateInEdges[neighbor].end(), i) != intermediateInEdges[neighbor].end());
		}
	}
	std::vector<std::tuple<size_t, int, size_t>> callStack;
	size_t i = 0;
	std::vector<size_t> index;
	std::vector<size_t> lowlink;
	std::vector<bool> onStack;
	std::vector<size_t> stack;
	index.resize(intermediateBigraphNodeIDs.size(), std::numeric_limits<size_t>::max());
	lowlink.resize(intermediateBigraphNodeIDs.size(), std::numeric_limits<size_t>::max());
	onStack.resize(intermediateBigraphNodeIDs.size(), false);
	size_t checknode = 0;
	size_t nextComponent = 0;
	componentNumber.resize(intermediateBigraphNodeIDs.size(), std::numeric_limits<size_t>::max());
	while (true)
	{
		if (callStack.size() == 0)
		{
			while (checknode < intermediateBigraphNodeIDs.size() && index[checknode] != std::numeric_limits<size_t>::max())
			{
				checknode++;
			}
			if (checknode == intermediateBigraphNodeIDs.size()) break;
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
				if (neighborI >= intermediateOutEdges[v].size()) goto endloop;
				assert(neighborI < intermediateOutEdges[v].size());
				w = intermediateOutEdges[v][neighborI];
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
				assert(neighborI < intermediateOutEdges[v].size());
				w = intermediateOutEdges[v][neighborI];
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
	for (size_t i = 0; i < intermediateOutEdges.size(); i++)
	{
		for (auto neighbor : intermediateOutEdges[i])
		{
			assert(componentNumber[neighbor] >= componentNumber[i]);
		}
		for (auto neighbor : intermediateInEdges[i])
		{
			assert(componentNumber[neighbor] <= componentNumber[i]);
		}
	}
}

size_t AlignmentGraph::ComponentSize() const
{
	return NodeSize();
}

size_t AlignmentGraph::SizeInBP() const
{
	return bpSize;
}

size_t AlignmentGraph::BigraphNodeID(size_t digraphNodeId) const
{
	assert(Finalized());
	size_t intermediate = digraphToIntermediate(digraphNodeId);
	return intermediateBigraphNodeIDs[intermediate];
}

AlignmentGraph::NodeEdgeIterator AlignmentGraph::OutNeighbors(size_t digraphNodeId) const
{
	assert(Finalized());
	assert(digraphNodeId < NodeSize());
	if (!firstOfIntermediates.get(digraphNodeId+1))
	{
		return AlignmentGraph::NodeEdgeIterator { digraphNodeId+1, nullptr, nullptr };
	}
	size_t intermediate = digraphToIntermediate(digraphNodeId);
	return AlignmentGraph::NodeEdgeIterator { std::numeric_limits<size_t>::max(), intermediateOutEdges[intermediate].data(), intermediateOutEdges[intermediate].data() + intermediateOutEdges[intermediate].size() };
}

AlignmentGraph::NodeEdgeIterator AlignmentGraph::InNeighbors(size_t digraphNodeId) const
{
	assert(digraphNodeId < NodeSize());
	if (!firstOfIntermediates.get(digraphNodeId))
	{
		return AlignmentGraph::NodeEdgeIterator { digraphNodeId - 1, nullptr, nullptr };
	}
	size_t intermediate = digraphToIntermediate(digraphNodeId);
	return AlignmentGraph::NodeEdgeIterator { std::numeric_limits<size_t>::max(), intermediateInEdges[intermediate].data(), intermediateInEdges[intermediate].data() + intermediateInEdges[intermediate].size() };
}

size_t AlignmentGraph::BigraphNodeCount() const
{
	return originalNodeName.size();
}

size_t AlignmentGraph::NodeOffset(size_t digraphNodeId) const
{
	assert(Finalized());
	size_t intermediate = digraphToIntermediate(digraphNodeId);
	// assert(intermediateDinodesStart[intermediate] <= digraphNodeId);
	// assert(intermediateDinodesCount(intermediate) + intermediateDinodesStart[intermediate] > digraphNodeId);
	size_t result = firstDinodeOffset[intermediate] + (digraphNodeId - intermediateDinodesStart[intermediate]) * SPLIT_NODE_SIZE;
	// assert(result < firstDinodeOffset[intermediate] + intermediateNodeLength(intermediate));
	return result;
}

size_t AlignmentGraph::ChainApproxPos(size_t bigraphNodeId) const
{
	return chainApproxPos[bigraphNodeId];
}

size_t AlignmentGraph::ChainNumber(size_t bigraphNodeId) const
{
	return chainNumber[bigraphNodeId];
}

size_t AlignmentGraph::ComponentNumber(size_t digraphNodeId) const
{
	assert(Finalized());
	size_t intermediate = digraphToIntermediate(digraphNodeId);
	assert(digraphNodeId >= intermediateDinodesStart[intermediate]);
	size_t extra = 0;
	if (!partOfStronglyConnectedComponent[intermediate])
	{
		extra = digraphNodeId - intermediateDinodesStart[intermediate];
	}
	return componentNumber[intermediate] + extra;
}

bool AlignmentGraph::Reverse(size_t digraphNodeId) const
{
	assert(Finalized());
	size_t intermediate = digraphToIntermediate(digraphNodeId);
	size_t bigraphNodeId = intermediateBigraphNodeIDs[intermediate];
	return reverse[bigraphNodeId];
}

bool AlignmentGraph::Linearizable(size_t digraphNodeId) const
{
	return !firstOfIntermediates.get(digraphNodeId);
}

bool AlignmentGraph::Finalized() const
{
	return finalized;
}

size_t AlignmentGraph::FirstAmbiguous() const
{
	return firstAmbiguous;
}

size_t AlignmentGraph::digraphToIntermediate(size_t digraphNodeId) const
{
	assert(digraphNodeId < NodeSize());
	assert(firstOfIntermediates.size() == NodeSize()+1);
	size_t result = firstOfIntermediates.rankOne(digraphNodeId) + (firstOfIntermediates.get(digraphNodeId) ? 1 : 0) - 1;
	assert(result < intermediateDinodesStart.size());
	// assert(digraphNodeId >= intermediateDinodesStart[result]);
	// assert(digraphNodeId <= intermediateLastDinode(result));
	return result;
}

size_t AlignmentGraph::intermediateLastDinode(size_t intermediate) const
{
	size_t result = intermediateDinodesStart[intermediate];
	result += 1;
	while (!firstOfIntermediates.get(result))
	{
		result += 1;
	}
	result -= 1;
	assert(digraphToIntermediate(result) == intermediate);
	assert(result == NodeSize()-1 || digraphToIntermediate(result+1) != intermediate);
	return result;
}

size_t AlignmentGraph::intermediateDinodesCount(size_t intermediate) const
{
	return intermediateLastDinode(intermediate) - intermediateDinodesStart[intermediate] + 1;
}

std::string AlignmentGraph::BigraphNodeSeq(size_t bigraphNodeId) const
{
	std::string result;
	result.reserve(BigraphNodeSize(bigraphNodeId));
	for (auto interNode : bigraphIntermediateList[bigraphNodeId])
	{
		size_t end = intermediateLastDinode(interNode)+1;
		for (size_t digraphNodeId = intermediateDinodesStart[interNode]; digraphNodeId < end; digraphNodeId++)
		{
			size_t len = NodeLength(digraphNodeId);
			for (size_t i = 0; i < len; i++)
			{
				result.push_back(NodeSequences(digraphNodeId, i));
			}
		}
	}
	assert(result.size() == BigraphNodeSize(bigraphNodeId));
	return result;
}

bool AlignmentGraph::AllNodeNamesAreNumbers() const
{
	return allNodeNamesAreNumbers;
}
