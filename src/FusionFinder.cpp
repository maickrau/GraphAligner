#include <mutex>
#include <thread>
#include <iostream>
#include <regex>
#include <fstream>
#include <string>
#include "BigraphToDigraph.h"
#include "vg.pb.h"
#include "AlignmentGraph.h"
#include "Aligner.h"
#include "GfaGraph.h"
#include "CommonUtils.h"
#include "fastqloader.h"
#include "GraphAlignerWrapper.h"
#include "ThreadReadAssertion.h"
#include "MummerSeeder.h"

struct FusionAlignment
{
	FusionAlignment() {}
	FusionAlignment(std::shared_ptr<vg::Alignment> alignment, std::string leftGene, std::string rightGene, int scoreDifference, std::string corrected) :
		alignment(alignment),
		leftGene(leftGene),
		rightGene(rightGene),
		scoreDifference(scoreDifference),
		corrected(corrected)
	{}
	std::shared_ptr<vg::Alignment> alignment;
	std::string leftGene;
	std::string rightGene;
	int scoreDifference;
	std::string corrected;
};

std::string geneFromTranscript(std::string transcript)
{
	std::regex generegex("[_ ]gene:(ENSG\\d{11}\\.\\d{1,2})[_ ]");
	std::smatch match;
	std::regex_search(transcript, match, generegex);
	assert(match.ready());
	assert(!match.empty());
	assert(match.size() >= 2);
	assert(match[1].matched);
	return std::string { match[1].first, match[1].second };
}

std::vector<std::pair<std::string, std::string>> loadPutativeFusions(std::string filename, int minPutativeSupport)
{
	std::ifstream file { filename };
	std::vector<std::pair<std::string, std::string>> result;
	while (file.good())
	{
		std::string left = "", right = "";
		int support = 0;
		file >> left >> right >> support;
		if (!file.good()) break;
		if (left == right) continue;
		if (support >= minPutativeSupport) result.emplace_back(left, right);
	}
	return result;
}

std::unordered_map<std::string, std::unordered_set<int>> getGeneBelongers(const std::vector<vg::Alignment>& alns, const GfaGraph& graph)
{
	std::unordered_map<std::string, std::unordered_set<int>> result;
	for (auto aln : alns)
	{
		std::string gene = geneFromTranscript(aln.name());
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			assert(graph.nodes.count(aln.path().mapping(i).position().node_id()) == 1);
			result[gene].insert(aln.path().mapping(i).position().node_id());
		}
	}
	return result;
}

GfaGraph getNonfusionGraph(const std::string& gene, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers)
{
	assert(geneBelongers.count(gene) == 1);
	return graph.GetSubgraph(geneBelongers.at(gene));
}

GfaGraph getFusionGraph(const std::string& leftGene, const std::string& rightGene, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers)
{
	assert(geneBelongers.count(leftGene) == 1);
	assert(geneBelongers.count(rightGene) == 1);
	GfaGraph result;
	result.edgeOverlap = 0;
	result.nodes[0] = 'N';
	result.nodes[1] = 'N';
	result.nodes[2] = 'N';
	result.nodes[3] = 'N';
	result.originalNodeName[0] = "DUMMY_MIDDLE";
	result.originalNodeName[1] = "DUMMY_MIDDLE";
	result.originalNodeName[2] = "DUMMY_MIDDLE";
	result.originalNodeName[3] = "DUMMY_MIDDLE";
	int nextNodeId = 4;
	for (int leftOrientationI = 0; leftOrientationI <= 1; leftOrientationI++)
	{
		for (int rightOrientationI = 0; rightOrientationI <= 1; rightOrientationI++)
		{
			int subgraphNumber = leftOrientationI * 2 + rightOrientationI;
			bool leftOrientation = leftOrientationI == 1 ? true : false;
			bool rightOrientation = rightOrientationI == 1 ? true : false;
			std::unordered_map<int, int> nodeStart;
			std::unordered_map<int, int> nodeEnd;
			assert(geneBelongers.count(leftGene) == 1);
			assert(geneBelongers.at(leftGene).size() >= 1);
			for (auto node : geneBelongers.at(leftGene))
			{
				nodeStart[node] = nextNodeId;
				auto seq = graph.nodes.at(node);
				for (size_t i = 0; i < seq.size(); i++)
				{
					if (i > 0) result.edges[NodePos { nextNodeId-1, true }].emplace_back(nextNodeId, true);
					result.nodes[nextNodeId] = seq[i];
					result.originalNodeName[nextNodeId] = graph.originalNodeName.at(node);
					result.edges[NodePos { nextNodeId, leftOrientation }].emplace_back(subgraphNumber, true);
					nextNodeId++;
				}
				nodeEnd[node] = nextNodeId-1;
			}
			for (auto node : geneBelongers.at(leftGene))
			{
				assert(nodeEnd.count(node) == 1);
				int nodeid = nodeEnd.at(node);
				if (graph.edges.count(NodePos { node, true }) == 1)
				{
					for (auto edge : graph.edges.at(NodePos { node, true }))
					{
						assert(graph.nodes.count(edge.id) == 1);
						if (geneBelongers.at(leftGene).count(edge.id) == 0) continue;
						assert(edge.end);
						assert(nodeStart.count(edge.id) == 1);
						int targetNodeId = nodeStart.at(edge.id);
						result.edges[NodePos { nodeid, true }].emplace_back(targetNodeId, true);
					}
				}
			}
			nodeStart.clear();
			nodeEnd.clear();
			for (auto node : geneBelongers.at(rightGene))
			{
				nodeStart[node] = nextNodeId;
				auto seq = graph.nodes.at(node);
				for (size_t i = 0; i < seq.size(); i++)
				{
					if (i > 0) result.edges[NodePos { nextNodeId-1, true }].emplace_back(nextNodeId, true);
					result.nodes[nextNodeId] = seq[i];
					result.originalNodeName[nextNodeId] = graph.originalNodeName.at(node);
					result.edges[NodePos { subgraphNumber, true }].emplace_back(nextNodeId, rightOrientation);
					nextNodeId++;
				}
				nodeEnd[node] = nextNodeId-1;
			}
			for (auto node : geneBelongers.at(rightGene))
			{
				assert(nodeEnd.count(node) == 1);
				int nodeid = nodeEnd.at(node);
				if (graph.edges.count(NodePos { node, true }) == 1)
				{
					for (auto edge : graph.edges.at(NodePos { node, true }))
					{
						assert(graph.nodes.count(edge.id) == 1);
						if (geneBelongers.at(rightGene).count(edge.id) == 0) continue;
						assert(edge.end);
						assert(nodeStart.count(edge.id) == 1);
						int targetNodeId = nodeStart.at(edge.id);
						result.edges[NodePos { nodeid, true }].emplace_back(targetNodeId, true);
					}
				}
			}
		}
	}
	return result;
}

std::string getCorrected(const vg::Alignment& aln, const GfaGraph& graph)
{
	std::string result;
	for (int i = 0; i < aln.path().mapping_size(); i++)
	{
		for (int j = 0; j < aln.path().mapping(i).edit_size(); j++)
		{
			std::string n = graph.nodes.at(aln.path().mapping(i).position().node_id());
			if (aln.path().mapping(i).position().is_reverse()) n = CommonUtils::ReverseComplement(n);
			result += n.substr(aln.path().mapping(i).position().offset(), aln.path().mapping(i).edit(j).from_length());
		}
	}
	return result;
}

void addBestAlnsOnePair(std::unordered_map<size_t, FusionAlignment>& bestAlns, const std::string& leftGene, const std::string& rightGene, const GfaGraph& fusiongraph, const std::vector<FastQ>& reads, const std::vector<size_t>& readIndices, double maxScoreFraction, int minFusionLen, bool optimal)
{
	auto alignmentGraph = DirectedGraph::BuildFromGFA(fusiongraph);
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, 100, true };
	for (auto readIndex : readIndices)
	{
		const auto& read = reads[readIndex];
		try
		{
			AlignmentResult alignments;
			if (optimal)
			{
				alignments = AlignOneWayDijkstra(alignmentGraph, read.seq_id, read.sequence, true, reusableState, true, false);
			}
			else
			{
				alignments = AlignOneWay(alignmentGraph, read.seq_id, read.sequence, 100, 100, true, reusableState, true, true, false, false, 0, 0, 0);
			}
			if (bestAlns.count(readIndex) == 1 && alignments.alignments[0].alignmentScore > bestAlns.at(readIndex).alignment->score()) continue;
			AddAlignment(read.seq_id, read.sequence, alignments.alignments[0]);
			replaceDigraphNodeIdsWithOriginalNodeIds(*alignments.alignments[0].alignment, alignmentGraph);
			if (alignments.alignments[0].alignment->score() > read.sequence.size() * maxScoreFraction) continue;
			int leftAlnSize = 0;
			int rightAlnSize = 0;
			bool crossedDummy = false;
			for (int i = 0; i < alignments.alignments[0].alignment->path().mapping_size(); i++)
			{
				if (alignments.alignments[0].alignment->path().mapping(i).position().name().substr(0, 12) == "DUMMY_MIDDLE")
				{
					crossedDummy = true;
					continue;
				}
				if (!crossedDummy)
				{
					leftAlnSize += alignments.alignments[0].alignment->path().mapping(i).edit(0).to_length();
				}
				else
				{
					rightAlnSize += alignments.alignments[0].alignment->path().mapping(i).edit(0).to_length();
				}
			}
			if (leftAlnSize < minFusionLen || rightAlnSize < minFusionLen) continue;
			if (bestAlns.count(readIndex) == 0 || alignments.alignments[0].alignment->score() < bestAlns.at(readIndex).alignment->score())
			{
				bestAlns[readIndex] = FusionAlignment { alignments.alignments[0].alignment, leftGene, rightGene, 0, getCorrected(*alignments.alignments[0].alignment, fusiongraph) };
			}
		}
		catch (ThreadReadAssertion::AssertionFailure& e)
		{
			reusableState.clear();
		}
	}
}

std::vector<FusionAlignment> getBestAlignments(const std::vector<std::pair<std::string, std::string>>& putativeFusions, const std::unordered_map<std::string, std::vector<size_t>>& hasSeeds, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers, const std::vector<FastQ>& allReads, double maxScoreFraction, int minFusionLen, int fusionPenalty, size_t numThreads, std::unordered_map<std::string, std::unordered_set<size_t>>& readsInNonfusionGraph)
{
	std::cerr << "get fusions" << std::endl;
	std::vector<std::thread> threads;
	size_t nextPair = 0;
	std::mutex nextPairMutex;
	std::vector<std::unordered_map<size_t, FusionAlignment>> bestFusionAlnsPerThread;
	std::mutex readsInNonfusionGraphMutex;
	bestFusionAlnsPerThread.resize(numThreads);
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&putativeFusions, &readsInNonfusionGraph, &readsInNonfusionGraphMutex, &bestFusionAlnsPerThread, &allReads, thread, maxScoreFraction, minFusionLen, &nextPair, &graph, &geneBelongers, &hasSeeds, &nextPairMutex]()
		{
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					if (nextPair == putativeFusions.size()) break;
					i = nextPair;
					std::cerr << "fusion " << nextPair << "/" << putativeFusions.size() << " " << putativeFusions[i].first << "-" << putativeFusions[i].second << std::endl;
					nextPair += 1;
				}
				assert(putativeFusions[i].first != putativeFusions[i].second);
				std::unordered_set<size_t> readsHere;
				// if (hasSeeds.count(putativeFusions[i].first) == 0) continue;
				// if (hasSeeds.count(putativeFusions[i].second) == 0) continue;
				// readsHere.resize(hasSeeds.at(putativeFusions[i].first).size() + hasSeeds.at(putativeFusions[i].second).size(), 0);
				// auto final = std::set_intersection(hasSeeds.at(putativeFusions[i].first).begin(), hasSeeds.at(putativeFusions[i].first).end(), hasSeeds.at(putativeFusions[i].second).begin(), hasSeeds.at(putativeFusions[i].second).end(), readsHere.begin());
				// readsHere.resize(final - readsHere.begin());
				if (hasSeeds.count(putativeFusions[i].first) == 1) readsHere.insert(hasSeeds.at(putativeFusions[i].first).begin(), hasSeeds.at(putativeFusions[i].first).end());
				if (hasSeeds.count(putativeFusions[i].second) == 1) readsHere.insert(hasSeeds.at(putativeFusions[i].second).begin(), hasSeeds.at(putativeFusions[i].second).end());
				{
					std::lock_guard<std::mutex> lock { readsInNonfusionGraphMutex };
					readsInNonfusionGraph[putativeFusions[i].first].insert(readsHere.begin(), readsHere.end());
					readsInNonfusionGraph[putativeFusions[i].second].insert(readsHere.begin(), readsHere.end());
				}
				std::vector<size_t> readIndices { readsHere.begin(), readsHere.end() };
				if (readIndices.size() == 0) continue;
				auto fusiongraph = getFusionGraph(putativeFusions[i].first, putativeFusions[i].second, graph, geneBelongers);
				addBestAlnsOnePair(bestFusionAlnsPerThread[thread], putativeFusions[i].first, putativeFusions[i].second, fusiongraph, allReads, readIndices, maxScoreFraction, minFusionLen, false);
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					std::cerr << "finished fusion " << i << "/" << putativeFusions.size() << " " << putativeFusions[i].first << "-" << putativeFusions[i].second << std::endl;
				}
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	threads.clear();
	std::cerr << "merge fusion results" << std::endl;
	std::unordered_map<size_t, FusionAlignment> bestFusionAlns;
	for (size_t i = 0; i < numThreads; i++)
	{
		for (auto pair : bestFusionAlnsPerThread[i])
		{
			if (bestFusionAlns.count(pair.first) == 0 || pair.second.alignment->score() < bestFusionAlns.at(pair.first).alignment->score())
			{
				bestFusionAlns[pair.first] = pair.second;
			}
		}
	}
	std::vector<std::unordered_map<size_t, FusionAlignment>> bestNonfusionAlnsPerThread;
	bestNonfusionAlnsPerThread.resize(numThreads);
	std::cerr << "get forbidden genes" << std::endl;
	std::unordered_set<std::string> forbiddenGenes;
	size_t totalCount = 0;
	size_t geneCount = 0;
	for (const auto& pair : readsInNonfusionGraph)
	{
		if (pair.second.size() > 0) geneCount += 1;
		totalCount += pair.second.size();
	}
	assert(geneCount > 0);
	size_t cutoff = (double)totalCount * 200.0 / (double)geneCount;
	for (const auto& pair : readsInNonfusionGraph)
	{
		if (pair.second.size() > cutoff) forbiddenGenes.insert(pair.first);
	}
	std::cerr << forbiddenGenes.size() << " forbidden genes:";
	for (auto gene : forbiddenGenes)
	{
		std::cerr << " " << gene;
	}
	std::cerr << std::endl;
	std::cerr << "get nonfusions" << std::endl;
	nextPair = 0;
	std::vector<std::string> nonFusionGenes;
	for (auto pair : readsInNonfusionGraph)
	{
		nonFusionGenes.push_back(pair.first);
	}
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&nonFusionGenes, &bestFusionAlns, fusionPenalty, &forbiddenGenes, &bestNonfusionAlnsPerThread, &allReads, thread, maxScoreFraction, minFusionLen, &nextPair, &graph, &geneBelongers, &readsInNonfusionGraph, &nextPairMutex]()
		{
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					if (nextPair == nonFusionGenes.size()) break;
					i = nextPair;
					std::cerr << "nonfusion " << nextPair << "/" << nonFusionGenes.size() << " " << nonFusionGenes[i] << " " << readsInNonfusionGraph.at(nonFusionGenes[i]).size() << " reads" << std::endl;
					nextPair += 1;
				}
				if (forbiddenGenes.count(nonFusionGenes[i]) == 1)
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					std::cerr << "forbidden nonfusion " << i << "/" << nonFusionGenes.size() << " " << nonFusionGenes[i] << std::endl;
					continue;
				}
				auto nonfusiongraph = getNonfusionGraph(nonFusionGenes[i], graph, geneBelongers);
				assert(readsInNonfusionGraph.count(nonFusionGenes[i]) == 1);
				std::vector<size_t> readIndices { readsInNonfusionGraph.at(nonFusionGenes[i]).begin(), readsInNonfusionGraph.at(nonFusionGenes[i]).end() };
				for (size_t j = readIndices.size()-1; j < readIndices.size(); j++)
				{
					if (bestFusionAlns.count(j) == 0)
					{
						std::swap(readIndices[j], readIndices.back());
						readIndices.pop_back();
						continue;
					}
					if (bestNonfusionAlnsPerThread[thread].count(j) == 1 && bestNonfusionAlnsPerThread[thread].at(j).alignment->score() <= bestFusionAlns.at(j).alignment->score() + fusionPenalty)
					{
						std::swap(readIndices[j], readIndices.back());
						readIndices.pop_back();
						continue;
					}
				}
				addBestAlnsOnePair(bestNonfusionAlnsPerThread[thread], nonFusionGenes[i], nonFusionGenes[i], nonfusiongraph, allReads, readIndices, 1, 0, false);
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					std::cerr << "finished nonfusion " << i << "/" << nonFusionGenes.size() << " " << nonFusionGenes[i] << std::endl;
				}
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	std::unordered_map<size_t, FusionAlignment> bestNonfusionAlns;
	std::cerr << "merge nonfusion results" << std::endl;
	for (size_t i = 0; i < numThreads; i++)
	{
		for (auto pair : bestNonfusionAlnsPerThread[i])
		{
			if (bestNonfusionAlns.count(pair.first) == 0 || pair.second.alignment->score() < bestNonfusionAlns.at(pair.first).alignment->score())
			{
				bestNonfusionAlns[pair.first] = pair.second;
			}
		}
	}
	std::cerr << "filter fusion by nonfusion" << std::endl;
	std::vector<FusionAlignment> result;
	for (auto aln : bestFusionAlns)
	{
		if (bestNonfusionAlns.count(aln.first) == 1)
		{
			if (bestNonfusionAlns.at(aln.first).alignment->score() <= aln.second.alignment->score() + fusionPenalty)
			{
				continue;
			}
			else
			{
				result.emplace_back(aln.second.alignment, aln.second.leftGene, aln.second.rightGene, aln.second.alignment->score() - bestNonfusionAlns.at(aln.first).alignment->score(), aln.second.corrected);
			}
		}
		else
		{
			result.emplace_back(aln.second.alignment, aln.second.leftGene, aln.second.rightGene, aln.second.alignment->sequence().size() - aln.second.alignment->score(), aln.second.corrected);
		}
	}
	return result;
}

void writeFusions(const std::vector<FusionAlignment>& result, std::ofstream& file)
{
	for (auto aln : result)
	{
		auto fusionaln = aln.alignment;
		int fusionIndex = -1;
		size_t leftLen = 0;
		size_t rightLen = 0;
		for (int i = 0; i < fusionaln->path().mapping_size(); i++)
		{
			if (fusionaln->path().mapping(i).position().name().substr(0, 12) == "DUMMY_MIDDLE")
			{
				fusionIndex = i;
				continue;
			}
			if (fusionIndex == -1)
			{
				leftLen += fusionaln->path().mapping(i).edit(0).to_length();
			}
			else
			{
				rightLen += fusionaln->path().mapping(i).edit(0).to_length();
			}
		}
		assert(fusionIndex != -1);
		assert(fusionIndex > 0);
		assert(fusionIndex < fusionaln->path().mapping_size() - 1);
		std::string leftName = fusionaln->path().mapping(fusionIndex-1).position().name();
		std::string rightName = fusionaln->path().mapping(fusionIndex+1).position().name();
		bool leftReverse = fusionaln->path().mapping(fusionIndex-1).position().is_reverse();
		bool rightReverse = fusionaln->path().mapping(fusionIndex+1).position().is_reverse();
		for (int i = fusionIndex-1; i >= 0; i--)
		{
			if (fusionaln->path().mapping(i).position().name() != leftName)
			{
				leftName = fusionaln->path().mapping(i).position().name();
				leftReverse = fusionaln->path().mapping(i).position().is_reverse();
				break;
			}
		}
		for (int i = fusionIndex+1; i < fusionaln->path().mapping_size(); i++)
		{
			if (fusionaln->path().mapping(i).position().name() != rightName)
			{
				rightName = fusionaln->path().mapping(i).position().name();
				rightReverse = fusionaln->path().mapping(i).position().is_reverse();
				break;
			}
		}
		if (fusionaln->path().mapping(fusionIndex).position().is_reverse())
		{
			std::swap(leftName, rightName);
			std::swap(leftReverse, rightReverse);
			leftReverse = !leftReverse;
			rightReverse = !rightReverse;
			std::swap(aln.leftGene, aln.rightGene);
		}
		file << fusionaln->name() << "\t" << ((double)fusionaln->score() / (double)fusionaln->sequence().size()) << "\t" << aln.scoreDifference << "\t" << aln.leftGene << "\t" << aln.rightGene << "\t" << leftLen << "\t" << leftName << "\t" << (leftReverse ? "-" : "+") << "\t" << rightName << "\t" << (rightReverse ? "-" : "+") << "\t" << rightLen << std::endl;
	}
}

std::unordered_map<std::string, std::vector<size_t>> getIntSeeds(const std::unordered_map<std::string, std::vector<std::string>>& hasSeeds, const std::vector<FastQ>& reads)
{
	std::unordered_map<std::string, std::vector<size_t>> result;
	for (size_t i = 0; i < reads.size(); i++)
	{
		if (hasSeeds.count(reads[i].seq_id) == 0) continue;
		for (const auto& gene : hasSeeds.at(reads[i].seq_id))
		{
			result[gene].push_back(i);
		}
	}
	return result;
}

std::unordered_map<std::string, std::vector<std::string>> loadPartialToTranscripts(std::string filename)
{
	std::ifstream file { filename };
	std::regex splitter("([^\\t]+)_pair\\d+_\\d+\\t([^\\t]+)\\t1");
	std::unordered_map<std::string, std::vector<std::string>> result;
	while (file.good())
	{
		std::string line;
		std::getline(file, line);
		if (!file.good()) break;
		std::smatch match;
		std::regex_search(line, match, splitter);
		if (match.empty()) continue;
		assert(match.size() == 3);
		std::string read { match[1].first, match[1].second };
		std::string transcript { match[2].first, match[2].second };
		result[read].push_back(geneFromTranscript(transcript));
	}
	return result;
}

void writeCorrected(const std::vector<FusionAlignment>& result, const GfaGraph& graph, std::ofstream& file)
{
	for (auto aln : result)
	{
		file << ">" << aln.alignment->name() << std::endl;
		file << aln.corrected << std::endl;
	}
}

// fake declarations because MinimizerSeeder.cpp is also linked
size_t charToInt(char c);
std::vector<bool> getValidChars();
extern std::vector<bool> validChar;

template <typename CallbackF>
void iterateKmers(const std::string& str, size_t kmerLength, CallbackF callback)
{
	if (str.size() < kmerLength) return;
	const size_t mask = ~(0xFFFFFFFFFFFFFFFF << (kmerLength * 2));
	assert(mask == pow(4, kmerLength)-1);
	size_t offset = 0;
start:
	while (offset < str.size() && !validChar[str[offset]]) offset++;
	if (offset + kmerLength > str.size()) return;
	size_t kmer = 0;
	for (size_t i = 0; i < kmerLength; i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset += i;
			goto start;
		}
		kmer <<= 2;
		kmer |= charToInt(str[offset+i]);
	}
	callback(offset + kmerLength-1, kmer);
	for (size_t i = kmerLength; offset+i < str.size(); i++)
	{
		if (!validChar[str[offset+i]])
		{
			offset += i;
			goto start;
		}
		kmer <<= 2;
		kmer &= mask;
		kmer |= charToInt(str[offset + i]);
		callback(offset + i, kmer);
	}
}

std::vector<std::vector<size_t>> getExtraGeneMatchKmerIndex(const std::vector<vg::Alignment>& transcripts)
{
	std::vector<std::vector<size_t>> kmerIndex;
	std::vector<std::pair<std::string, size_t>> transcriptptrs;
	for (size_t i = 0; i < transcripts.size(); i++)
	{
		transcriptptrs.emplace_back(geneFromTranscript(transcripts[i].name()), i);
	}
	std::sort(transcriptptrs.begin(), transcriptptrs.end(), [](const std::pair<std::string, size_t>& left, const std::pair<std::string, size_t>& right) { return left.first < right.first; });
	kmerIndex.resize(pow(4, 11));
	std::unordered_set<size_t> currentKmers;
	std::string currentGene;
	size_t currentNameIndex = -1;
	for (auto& pair : transcriptptrs)
	{
		if (geneFromTranscript(transcripts[pair.second].name()) != currentGene)
		{
			if (currentNameIndex != -1)
			{
				for (auto kmer : currentKmers)
				{
					kmerIndex[kmer].push_back(currentNameIndex);
				}
			}
			currentKmers.clear();
			currentGene = geneFromTranscript(transcripts[pair.second].name());
			currentNameIndex = pair.second;
		}
		iterateKmers(transcripts[pair.second].sequence(), 11, [&currentKmers](size_t offset, size_t kmer)
		{
			currentKmers.insert(kmer);
		});
		iterateKmers(CommonUtils::ReverseComplement(transcripts[pair.second].sequence()), 11, [&currentKmers](size_t offset, size_t kmer)
		{
			currentKmers.insert(kmer);
		});
	}
	for (auto kmer : currentKmers)
	{
		kmerIndex[kmer].push_back(currentNameIndex);
	}
	return kmerIndex;
}

std::unordered_map<std::string, std::unordered_set<size_t>> getExtraGeneMatches(const std::vector<vg::Alignment>& transcripts, const std::vector<std::vector<size_t>>& kmerIndex, const std::vector<FastQ>& reads, const size_t numThreads)
{
	std::unordered_map<std::string, std::unordered_set<size_t>> result;
	std::mutex resultMutex;
	std::mutex readMutex;
	size_t nextRead = 0;
	std::vector<std::thread> threads;
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&nextRead, &transcripts, &reads, &result, &resultMutex, &readMutex, &kmerIndex]()
		{
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> guard { readMutex };
					i = nextRead;
					nextRead += 1;
				}
				if (i >= reads.size()) break;
				std::unordered_map<size_t, size_t> lastMatch;
				std::unordered_map<size_t, size_t> matchSize;
				iterateKmers(reads[i].sequence, 11, [&lastMatch, &matchSize, &kmerIndex](size_t offset, size_t kmer)
				{
					for (auto index : kmerIndex[kmer])
					{
						size_t addition = std::min(offset - lastMatch[index], (size_t)11);
						lastMatch[index] = offset;
						matchSize[index] += addition;
					}
				});
				for (auto pair : matchSize)
				{
					std::vector<std::string> insertions;
					if (pair.second >= 1000 || pair.second >= reads[i].sequence.size() * .25)
					{
						insertions.push_back(geneFromTranscript(transcripts[pair.first].name()));
					}
					{
						std::lock_guard<std::mutex> guard { resultMutex };
						for (const auto& gene : insertions)
						{
							result[gene].insert(i);
						}
					}
				}
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	return result;
}

template <typename F>
void runInReadChunks(std::string filename, size_t chunkSize, F callback)
{
	std::vector<FastQ> vec;
	std::cerr << "stream reads from " << filename << std::endl;
	FastQ::streamFastqFromFile(filename, false, [&vec, &callback, chunkSize](FastQ read)
	{
		vec.emplace_back(std::move(read));
		if (vec.size() == chunkSize)
		{
			callback(vec);
			vec.clear();
		}
		assert(vec.size() < chunkSize);
	});
	std::cerr << "done streaming reads, size " << vec.size() << std::endl;
	if (vec.size() > 0)
	{
		callback(vec);
		vec.clear();
	}
	std::cerr << "done reading" << std::endl;
}

int main(int argc, char** argv)
{
	std::cerr << "Fusion finder " << VERSION << std::endl;

	std::string graphFile { argv[1] };
	std::string putativeFusionsFile { argv[2] };
	std::string partialMatrixFile { argv[3] };
	std::string transcriptAlignmentFile { argv[4] };
	std::string readFile { argv[5] };
	int minPutativeSupport = std::stoi(argv[6]);
	double maxScoreFraction = std::stod(argv[7]);
	int minFusionLen = std::stoi(argv[8]);
	int fusionPenalty = std::stoi(argv[9]);
	int numThreads = std::stoi(argv[10]);
	std::string resultFusionFile { argv[11] };
	std::string correctedReadsFile { argv[12] };
	int chunkSize = std::stoi(argv[13]);

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(graphFile);
	std::cerr << "load putative fusions" << std::endl;
	auto putativeFusions = loadPutativeFusions(putativeFusionsFile, minPutativeSupport);
	std::cerr << "load transcript alignments" << std::endl;
	auto transcripts = CommonUtils::LoadVGAlignments(transcriptAlignmentFile);
	std::cerr << "get gene belongers" << std::endl;
	auto geneBelongers = getGeneBelongers(transcripts, graph);
	std::cerr << "load partial assignments" << std::endl;
	auto hasSeeds = loadPartialToTranscripts(partialMatrixFile);
	std::cerr << "build extra gene-match index" << std::endl;
	auto kmerIndex = getExtraGeneMatchKmerIndex(transcripts);
	std::cerr << "run in chunks" << std::endl;
	size_t readsProcessed = 0;
	std::ofstream correctedOut { correctedReadsFile };
	std::ofstream resultFusionsOut { resultFusionFile };
	runInReadChunks(readFile, chunkSize, [&readsProcessed, &correctedOut, &resultFusionsOut, &kmerIndex, &transcripts, numThreads, &putativeFusions, &hasSeeds, &graph, &geneBelongers, maxScoreFraction, minFusionLen, fusionPenalty](const std::vector<FastQ>& reads)
	{
		std::cerr << "reads " << readsProcessed << " - " << readsProcessed + reads.size() << std::endl;
		readsProcessed += reads.size();
		std::cerr << "get chunk extra gene-matches" << std::endl;
		auto extraGeneMatches = getExtraGeneMatches(transcripts, kmerIndex, reads, numThreads);
		std::cerr << "get chunk partial assignments" << std::endl;
		auto intSeeds = getIntSeeds(hasSeeds, reads);
		std::cerr << "get chunk alns" << std::endl;
		auto bestAlns = getBestAlignments(putativeFusions, intSeeds, graph, geneBelongers, reads, maxScoreFraction, minFusionLen, fusionPenalty, numThreads, extraGeneMatches);
		std::cerr << "write chunk fusions" << std::endl;
		writeFusions(bestAlns, resultFusionsOut);
		std::cerr << "write chunk corrected reads" << std::endl;
		writeCorrected(bestAlns, graph, correctedOut);
	});
}