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

GfaGraph getNonfusionGraph(std::string gene, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers)
{
	assert(geneBelongers.count(gene) == 1);
	return graph.GetSubgraph(geneBelongers.at(gene));
}

GfaGraph getFusionGraph(std::string leftGene, std::string rightGene, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers)
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
	for (int leftOrientationI = 0; leftOrientationI <= 1; leftOrientationI++)
	{
		for (int rightOrientationI = 0; rightOrientationI <= 1; rightOrientationI++)
		{
			int subgraphNumber = leftOrientationI * 2 + rightOrientationI;
			bool leftOrientation = leftOrientationI == 1 ? true : false;
			bool rightOrientation = rightOrientationI == 1 ? true : false;
			assert(geneBelongers.count(leftGene) == 1);
			assert(geneBelongers.at(leftGene).size() >= 1);
			for (auto node : geneBelongers.at(leftGene))
			{
				int nodeid = 4 + subgraphNumber * graph.nodes.size() + node;
				assert(graph.nodes.count(node) == 1);
				assert(graph.originalNodeName.count(node) == 1);
				result.nodes[nodeid] = graph.nodes.at(node);
				result.originalNodeName[nodeid] = graph.originalNodeName.at(node);
				result.edges[NodePos { nodeid, leftOrientation }].emplace_back(subgraphNumber, true);
				if (graph.edges.count(NodePos { node, true }) == 1)
				{
					for (auto edge : graph.edges.at(NodePos { node, true }))
					{
						assert(graph.nodes.count(edge.id) == 1);
						if (geneBelongers.at(leftGene).count(edge.id) == 0) continue;
						assert(edge.end);
						int targetNodeId = 4 + subgraphNumber * graph.nodes.size() + edge.id;
						result.edges[NodePos { nodeid, true }].emplace_back(targetNodeId, true);
					}
				}
			}
			assert(geneBelongers.count(rightGene) == 1);
			assert(geneBelongers.at(rightGene).size() >= 1);
			for (auto node : geneBelongers.at(rightGene))
			{
				int nodeid = 4 + (4 + subgraphNumber) * graph.nodes.size() + node;
				assert(graph.nodes.count(node) == 1);
				assert(graph.originalNodeName.count(node) == 1);
				result.nodes[nodeid] = graph.nodes.at(node);
				result.originalNodeName[nodeid] = graph.originalNodeName.at(node);
				result.edges[NodePos { subgraphNumber, true }].emplace_back(nodeid, rightOrientation);
				if (graph.edges.count(NodePos { node, true }) == 1)
				{
					for (auto edge : graph.edges.at(NodePos { node, true }))
					{
						assert(graph.nodes.count(edge.id) == 1);
						if (geneBelongers.at(rightGene).count(edge.id) == 0) continue;
						assert(edge.end);
						int targetNodeId = 4 + (4 + subgraphNumber) * graph.nodes.size() + edge.id;
						result.edges[NodePos { nodeid, true }].emplace_back(targetNodeId, true);
					}
				}
			}
		}
	}
	return result;
}

void addBestAlnsOnePair(std::unordered_map<std::string, std::tuple<AlignmentResult, std::string, std::string>>& bestAlns, std::string leftGene, std::string rightGene, const GfaGraph& fusiongraph, const std::vector<FastQ>& reads, double maxScoreFraction, int minFusionLen)
{
	auto alignmentGraph = DirectedGraph::BuildFromGFA(fusiongraph, true);
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, 1000, true };
	for (auto read : reads)
	{
		try
		{
			auto alignments = AlignOneWay(alignmentGraph, read.seq_id, read.sequence, 1000, 1000, true, reusableState, true, true);
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
			if (bestAlns.count(read.seq_id) == 0 || alignments.alignments[0].alignment->score() < std::get<0>(bestAlns.at(read.seq_id)).alignments[0].alignment->score())
			{
				bestAlns[read.seq_id] = std::make_tuple(alignments, leftGene, rightGene);
			}
		}
		catch (ThreadReadAssertion::AssertionFailure& e)
		{
			reusableState.clear();
		}
	}
}

std::vector<std::tuple<std::shared_ptr<vg::Alignment>, std::string, std::string, int>> getBestAlignments(const std::vector<std::pair<std::string, std::string>>& putativeFusions, const std::unordered_map<std::string, std::vector<size_t>>& hasSeeds, const GfaGraph& graph, const std::unordered_map<std::string, std::unordered_set<int>>& geneBelongers, const std::vector<FastQ>& allReads, double maxScoreFraction, int minFusionLen, int fusionPenalty, size_t numThreads)
{
	std::cerr << "get fusions" << std::endl;
	std::vector<std::thread> threads;
	size_t nextPair = 0;
	std::mutex nextPairMutex;
	std::vector<std::unordered_map<std::string, std::tuple<AlignmentResult, std::string, std::string>>> bestFusionAlnsPerThread;
	std::unordered_map<std::string, std::unordered_set<size_t>> readsInNonfusionGraph;
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
					std::cerr << "fusion " << nextPair << "/" << putativeFusions.size() << std::endl;
					nextPair += 1;
				}
				auto fusiongraph = getFusionGraph(putativeFusions[i].first, putativeFusions[i].second, graph, geneBelongers);
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
				std::vector<FastQ> reads;
				for (auto index : readsHere)
				{
					reads.push_back(allReads[index]);
				}
				addBestAlnsOnePair(bestFusionAlnsPerThread[thread], putativeFusions[i].first, putativeFusions[i].second, fusiongraph, reads, maxScoreFraction, minFusionLen);
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	threads.clear();
	std::cerr << "merge fusion results" << std::endl;
	std::unordered_map<std::string, std::tuple<AlignmentResult, std::string, std::string>> bestFusionAlns;
	for (size_t i = 0; i < numThreads; i++)
	{
		for (auto pair : bestFusionAlnsPerThread[i])
		{
			if (bestFusionAlns.count(pair.first) == 0 || std::get<0>(pair.second).alignments[0].alignment->score() < std::get<0>(bestFusionAlns.at(pair.first)).alignments[0].alignment->score())
			{
				bestFusionAlns[pair.first] = pair.second;
			}
		}
	}
	std::vector<std::unordered_map<std::string, std::tuple<AlignmentResult, std::string, std::string>>> bestNonfusionAlnsPerThread;
	bestNonfusionAlnsPerThread.resize(numThreads);
	std::cerr << "get nonfusions" << std::endl;
	nextPair = 0;
	std::vector<std::string> fusionGenes;
	for (auto pair : readsInNonfusionGraph)
	{
		fusionGenes.push_back(pair.first);
	}
	for (size_t thread = 0; thread < numThreads; thread++)
	{
		threads.emplace_back([&fusionGenes, &bestNonfusionAlnsPerThread, &allReads, thread, maxScoreFraction, minFusionLen, &nextPair, &graph, &geneBelongers, &readsInNonfusionGraph, &nextPairMutex]()
		{
			while (true)
			{
				size_t i = 0;
				{
					std::lock_guard<std::mutex> lock { nextPairMutex };
					if (nextPair == fusionGenes.size()) break;
					i = nextPair;
					std::cerr << "nonfusion " << nextPair << "/" << fusionGenes.size() << std::endl;
					nextPair += 1;
				}
				auto nonfusiongraph = getNonfusionGraph(fusionGenes[i], graph, geneBelongers);
				assert(readsInNonfusionGraph.count(fusionGenes[i]) == 1);
				std::vector<FastQ> reads;
				for (auto index : readsInNonfusionGraph.at(fusionGenes[i]))
				{
					reads.push_back(allReads[index]);
				}
				addBestAlnsOnePair(bestNonfusionAlnsPerThread[thread], fusionGenes[i], fusionGenes[i], nonfusiongraph, reads, 1, 0);
			}
		});
	}
	for (size_t i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}
	std::unordered_map<std::string, std::tuple<AlignmentResult, std::string, std::string>> bestNonfusionAlns;
	std::cerr << "merge nonfusion results" << std::endl;
	for (size_t i = 0; i < numThreads; i++)
	{
		for (auto pair : bestNonfusionAlnsPerThread[i])
		{
			if (bestNonfusionAlns.count(pair.first) == 0 || std::get<0>(pair.second).alignments[0].alignment->score() < std::get<0>(bestNonfusionAlns.at(pair.first)).alignments[0].alignment->score())
			{
				bestNonfusionAlns[pair.first] = pair.second;
			}
		}
	}
	std::cerr << "filter fusion by nonfusion" << std::endl;
	std::vector<std::tuple<std::shared_ptr<vg::Alignment>, std::string, std::string, int>> result;
	for (auto aln : bestFusionAlns)
	{
		if (bestNonfusionAlns.count(aln.first) == 1)
		{
			if (std::get<0>(bestNonfusionAlns.at(aln.first)).alignments[0].alignment->score() <= std::get<0>(aln.second).alignments[0].alignment->score() + fusionPenalty)
			{
				continue;
			}
			else
			{
				result.emplace_back(std::get<0>(aln.second).alignments[0].alignment, std::get<1>(aln.second), std::get<2>(aln.second), std::get<0>(aln.second).alignments[0].alignment->score() - std::get<0>(bestNonfusionAlns.at(aln.first)).alignments[0].alignment->score());
			}
		}
		else
		{
			result.emplace_back(std::get<0>(aln.second).alignments[0].alignment, std::get<1>(aln.second), std::get<2>(aln.second), std::get<0>(aln.second).alignments[0].alignment->sequence().size() - std::get<0>(aln.second).alignments[0].alignment->score());
		}
	}
	return result;
}

void writeFusions(const std::vector<std::tuple<std::shared_ptr<vg::Alignment>, std::string, std::string, int>>& result, std::string filename)
{
	std::ofstream file { filename };
	for (auto pair : result)
	{
		auto fusionaln = std::get<0>(pair);
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
		if (fusionaln->path().mapping(fusionIndex).position().is_reverse())
		{
			std::swap(leftName, rightName);
			std::swap(leftReverse, rightReverse);
			leftReverse = !leftReverse;
			rightReverse = !rightReverse;
			std::swap(std::get<1>(pair), std::get<2>(pair));
		}
		file << fusionaln->name() << "\t" << ((double)fusionaln->score() / (double)fusionaln->sequence().size()) << "\t" << std::get<3>(pair) << "\t" << std::get<1>(pair) << "\t" << std::get<2>(pair) << "\t" << leftLen << "\t" << leftName << "\t" << (leftReverse ? "-" : "+") << "\t" << rightName << "\t" << (rightReverse ? "-" : "+") << "\t" << rightLen << std::endl;
	}
}

std::unordered_map<std::string, std::vector<size_t>> loadPartialToTranscripts(std::string filename, const std::vector<FastQ>& reads)
{
	std::unordered_map<std::string, size_t> readIndex;
	for (size_t i = 0; i < reads.size(); i++)
	{
		readIndex[reads[i].seq_id] = i;
	}
	std::ifstream file { filename };
	std::regex splitter("([^\\t]+)_pair\\d+_\\d+\\t([^\\t]+)\\t1");
	std::unordered_map<std::string, std::vector<size_t>> result;
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
		assert(readIndex.count(read) == 1);
		result[geneFromTranscript(transcript)].push_back(readIndex[read]);
	}
	return result;
}

void printCorrected(std::ofstream& str, const std::unordered_map<std::string, std::string>& mapping, const GfaGraph& graph, const vg::Alignment& v)
{
	str << ">" << v.name() << std::endl;
	for (int i = 0; i < v.path().mapping_size(); i++)
	{
		std::string sequence;
		if (mapping.count(v.path().mapping(i).position().name()) == 0)
		{
			assert(v.path().mapping(i).position().name() == "DUMMY_MIDDLE");
			sequence = "N";
		}
		else
		{
			sequence = mapping.at(v.path().mapping(i).position().name());
		}
		int len = 0;
		for (int j = 0; j < v.path().mapping(i).edit_size(); j++)
		{
			len += v.path().mapping(i).edit(j).from_length();
		}
		if (v.path().mapping(i).position().is_reverse())
		{
			sequence = CommonUtils::ReverseComplement(sequence);
		}
		if (v.path().mapping(i).position().offset() > 0)
		{
			sequence = sequence.substr(v.path().mapping(i).position().offset());
		}
		sequence = sequence.substr(0, len);
		str << sequence;
	}
	str << std::endl;
}

void writeCorrected(const std::vector<std::tuple<std::shared_ptr<vg::Alignment>, std::string, std::string, int>>& result, const GfaGraph& graph, std::string filename)
{
	std::ofstream file { filename };
	std::unordered_map<std::string, std::string> mapping;
	for (auto node : graph.nodes)
	{
		assert(graph.originalNodeName.count(node.first) == 1);
		mapping[graph.originalNodeName.at(node.first)] = node.second;
	}
	for (auto aln : result)
	{
		printCorrected(file, mapping, graph, *std::get<0>(aln));
	}
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

	std::cerr << "load graph" << std::endl;
	auto graph = GfaGraph::LoadFromFile(graphFile);
	std::cerr << "load putative fusions" << std::endl;
	auto putativeFusions = loadPutativeFusions(putativeFusionsFile, minPutativeSupport);
	std::cerr << "load reads" << std::endl;
	auto reads = loadFastqFromFile(readFile);
	std::cerr << "load partial assignments" << std::endl;
	auto hasSeeds = loadPartialToTranscripts(partialMatrixFile, reads);
	std::cerr << "load transcript alignments" << std::endl;
	auto transcripts = CommonUtils::LoadVGAlignments(transcriptAlignmentFile);
	std::cerr << "get gene belongers" << std::endl;
	auto geneBelongers = getGeneBelongers(transcripts, graph);
	std::cerr << "get alns" << std::endl;
	auto bestAlns = getBestAlignments(putativeFusions, hasSeeds, graph, geneBelongers, reads, maxScoreFraction, minFusionLen, fusionPenalty, numThreads);
	std::cerr << "write fusions" << std::endl;
	writeFusions(bestAlns, resultFusionFile);
	std::cerr << "write corrected reads" << std::endl;
	writeCorrected(bestAlns, graph, correctedReadsFile);
}