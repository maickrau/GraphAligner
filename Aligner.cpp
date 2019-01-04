#include <mutex>
#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include "Aligner.h"
#include "CommonUtils.h"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerCommon.h"
#include "GraphAlignerBitvectorFull.h"
#include "GraphAlignerCellbycellFull.h"

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

void wabiExperiments(AlignerParams params)
{
	std::cout << "load graph" << std::endl;
	auto alignmentGraph = DirectedGraph::StreamGFAGraphFromFile(params.graphFile);
	std::cout << "load reads" << std::endl;
	auto fastqs = loadFastqFromFile(params.fastqFile);

	std::cout << "preprocess graph" << std::endl;
	auto componentOrder = alignmentGraph.TopologicalOrderOfComponents();
	bool isAcyclic = true;
	bool isTree = true;
	for (size_t i = 0; i < componentOrder.first.size(); i++)
	{
		if (componentOrder.first[i].size() > 1)
		{
			isAcyclic = false;
			break;
		}
		for (auto neighbor : alignmentGraph.inNeighbors[componentOrder.first[i][0]])
		{
			if (neighbor == componentOrder.first[i][0])
			{
				isAcyclic = false;
				break;
			}
		}
	}
	if (isAcyclic)
	{
		for (size_t i = 0; i < alignmentGraph.NodeSize(); i++)
		{
			if (alignmentGraph.inNeighbors[i].size() > 1)
			{
				isTree = false;
				break;
			}
		}
	}
	else
	{
		isTree = false;
	}

	if (isTree)
	{
		assert(isAcyclic);
		std::cout << "The graph is linear / a tree / a forest" << std::endl;
	}
	else if (isAcyclic)
	{
		std::cout << "The graph is a DAG" << std::endl;
	}
	else
	{
		std::cout << "The graph is cyclic" << std::endl;
	}
	size_t numEdges = 0;
	for (const auto& list : alignmentGraph.inNeighbors)
	{
		numEdges += list.size();
	}
	std::cout << "Collapsed nodes: " << alignmentGraph.NodeSize() << std::endl;
	std::cout << "Collapsed edges: " << numEdges << std::endl;
	std::cout << "Nodes: " << alignmentGraph.SizeInBp() << std::endl;
	std::cout << "Edges: " << numEdges + alignmentGraph.SizeInBp() - alignmentGraph.NodeSize() << std::endl;
	std::cout << "BPs: " << alignmentGraph.SizeInBp() << std::endl;

	GraphAlignerCommon<size_t, int32_t, uint64_t>::Params alignerParams { 0, 0, alignmentGraph, 0, false, false };

	GraphAlignerBitvectorFull<size_t, int32_t, uint64_t> bv { alignerParams, componentOrder.first, componentOrder.second };
	GraphAlignerCellbycellFull<size_t, int32_t, uint64_t> cbc { alignerParams, componentOrder.first, componentOrder.second };

	std::vector<int32_t> bvScores;
	std::vector<int32_t> cbcScores;
	bvScores.reserve(fastqs.size());
	cbcScores.reserve(fastqs.size());
	std::cout << "start bitvector alignment" << std::endl;
	auto bvtimeStart = std::chrono::steady_clock::now();
	for (auto fastq : fastqs)
	{
		size_t score;
		if (params.linear)
		{
			score = bv.alignAndGetScoreLinear(fastq.sequence);
		}
		else if (isAcyclic)
		{
			score = bv.alignAndGetScoreAcyclic(fastq.sequence);
		}
		else
		{
			score = bv.alignAndGetScore(fastq.sequence);
		}
		bvScores.push_back(score);
	}
	auto bvtimeEnd = std::chrono::steady_clock::now();
	size_t bitvectorMicroseconds = std::chrono::duration_cast<std::chrono::microseconds>(bvtimeEnd - bvtimeStart).count();
	std::cout << "bitvector took " << bitvectorMicroseconds << "us" << std::endl;
	std::cout << "start cellbycell alignment" << std::endl;
	auto cbctimeStart = std::chrono::steady_clock::now();
	for (auto fastq : fastqs)
	{
		size_t score;
		if (isAcyclic)
		{
			score = cbc.alignAndGetScoreAcyclic(fastq.sequence);
		}
		else
		{
			score = cbc.alignAndGetScore(fastq.sequence);
		}
		cbcScores.push_back(score);
	}
	auto cbctimeEnd = std::chrono::steady_clock::now();
	size_t cellbycellMicroseconds = std::chrono::duration_cast<std::chrono::microseconds>(cbctimeEnd - cbctimeStart).count();
	std::cout << "cellbycell took " << cellbycellMicroseconds << "us" << std::endl;
	std::cout << "ratio: " << ((double)cellbycellMicroseconds / (double)bitvectorMicroseconds) << std::endl;
	bool scoresOK = true;
	for (size_t i = 0; i < bvScores.size(); i++)
	{
		if (bvScores[i] != cbcScores[i])
		{
			std::cout << "SCORES DON'T MATCH for read " << fastqs[i].seq_id << ": bitvector " << bvScores[i] << " cellbycell " << cbcScores[i] << std::endl;
			scoresOK = false;
		}
	}
	if (scoresOK) std::cout << "scores match" << std::endl;
}