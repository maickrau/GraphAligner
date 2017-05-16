#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "TopologicalSort.h"
#include "SubgraphFromSeed.h"
#include "GraphAligner.h"
#include "mfvs_graph.h"

class BufferedWriter : std::ostream
{
public:
	class FlushClass {};
	BufferedWriter(std::ostream& stream) : stream(stream) {};
	template <typename T>
	BufferedWriter& operator<<(T obj)
	{
		stringstream << obj;
		return *this;
	}
	BufferedWriter& operator<<(FlushClass f)
	{
		flush();
		return *this;
	}
	void flush()
	{
		stringstream << std::endl;
		stream << stringstream.str();
		stringstream.str("");
	}
	static FlushClass Flush;
private:
	std::ostream& stream;
	std::stringstream stringstream;
};

size_t GraphSizeInBp(const vg::Graph& graph)
{
	size_t result = 0;
	for (int i = 0; i < graph.node_size(); i++)
	{
		result += graph.node(i).sequence().size();
	}
	return result;
}

vg::Graph mergeGraphs(const std::vector<vg::Graph>& parts)
{
	vg::Graph newGraph;
	std::vector<const vg::Node*> allNodes;
	std::vector<const vg::Edge*> allEdges;
	for (size_t i = 0; i < parts.size(); i++)
	{
		for (int j = 0; j < parts[i].node_size(); j++)
		{
			allNodes.push_back(&parts[i].node(j));
		}
		for (int j = 0; j < parts[i].edge_size(); j++)
		{
			allEdges.push_back(&parts[i].edge(j));
		}
	}
	for (size_t i = 0; i < allNodes.size(); i++)
	{
		auto node = newGraph.add_node();
		node->set_id(allNodes[i]->id());
		node->set_sequence(allNodes[i]->sequence());
		node->set_name(allNodes[i]->name());
	}
	for (size_t i = 0; i < allEdges.size(); i++)
	{
		auto edge = newGraph.add_edge();
		edge->set_from(allEdges[i]->from());
		edge->set_to(allEdges[i]->to());
		edge->set_from_start(allEdges[i]->from_start());
		edge->set_to_end(allEdges[i]->to_end());
		edge->set_overlap(allEdges[i]->overlap());
	}
	return newGraph;
}

int numberOfVerticesOutOfOrder(const vg::Graph& vggraph)
{
	std::map<int, int> ids;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		ids[vggraph.node(i).id()] = i;
	}
	std::set<int> outOfOrderSet;
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		assert(ids.count(vggraph.edge(i).to()) > 0);
		assert(ids.count(vggraph.edge(i).from()) > 0);
		if (ids[vggraph.edge(i).to()] <= ids[vggraph.edge(i).from()]) outOfOrderSet.insert(vggraph.edge(i).to());
	}
	return outOfOrderSet.size();
}

vg::Graph OrderByFeedbackVertexset(const vg::Graph& vggraph)
{
	BufferedWriter output { std::cout };
	mfvs::Graph mfvsgraph { vggraph.node_size() };
	std::map<int, int> ids;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		mfvsgraph.addVertex(i);
		ids[vggraph.node(i).id()] = i;
		output << vggraph.node(i).id() << " ";
	}
	output << BufferedWriter::Flush;
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		mfvsgraph.addEdge(ids[vggraph.edge(i).from()], ids[vggraph.edge(i).to()]);
	}
	std::cout<<"Before the removal of MVFS, is the graph acyclic:"<< mfvsgraph.isAcyclic()<<std::endl;
	auto vertexSetvector = mfvsgraph.minimumFeedbackVertexSet();
	std::set<int> vertexset { vertexSetvector.begin(), vertexSetvector.end() };
	output << "feedback vertex set size: " << vertexSetvector.size() << BufferedWriter::Flush;
	for(int i=0; i<vertexSetvector.size(); ++i){
		std::cout << vertexSetvector[i] << ' ' ;
		mfvsgraph.deleteVertex(vertexSetvector[i]);
	}
	std::cout<<"After the removal of MVFS, is the graph acyclic:"<< mfvsgraph.isAcyclic()<<std::endl;
	vg::Graph graphWithoutVertexset;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		if (vertexset.count(i) > 0) continue;
		auto node = graphWithoutVertexset.add_node();
		node->set_id(vggraph.node(i).id());
		node->set_sequence(vggraph.node(i).sequence());
		node->set_name(vggraph.node(i).name());
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		if (vertexset.count(ids[vggraph.edge(i).from()]) > 0 || vertexset.count(ids[vggraph.edge(i).to()]) > 0) continue;
		auto edge = graphWithoutVertexset.add_edge();
		edge->set_from(vggraph.edge(i).from());
		edge->set_to(vggraph.edge(i).to());
		edge->set_from_start(vggraph.edge(i).from_start());
		edge->set_to_end(vggraph.edge(i).to_end());
		edge->set_overlap(vggraph.edge(i).overlap());
	}
	auto order = topologicalSort(graphWithoutVertexset);
	vg::Graph resultGraph;
	for (size_t i = 0; i < vertexSetvector.size(); i++)
	{
		auto node = resultGraph.add_node();
		node->set_id(vggraph.node(vertexSetvector[i]).id());
		node->set_sequence(vggraph.node(vertexSetvector[i]).sequence());
		node->set_name(vggraph.node(vertexSetvector[i]).name());
		output << node->id() << " ";
	}
	for (size_t i = 0; i < order.size(); i++)
	{
		auto node = resultGraph.add_node();
		node->set_id(graphWithoutVertexset.node(order[i]).id());
		node->set_sequence(graphWithoutVertexset.node(order[i]).sequence());
		node->set_name(graphWithoutVertexset.node(order[i]).name());
		output << node->id() << " ";
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		auto edge = resultGraph.add_edge();
		edge->set_from(vggraph.edge(i).from());
		edge->set_to(vggraph.edge(i).to());
		edge->set_from_start(vggraph.edge(i).from_start());
		edge->set_to_end(vggraph.edge(i).to_end());
		edge->set_overlap(vggraph.edge(i).overlap());
	}
	output << BufferedWriter::Flush;
#ifndef NDEBUG
	std::set<int> existingNodeIds;
	for (int i = 0; i < resultGraph.node_size(); i++)
	{
		assert(existingNodeIds.count(resultGraph.node(i).id()) == 0);
		existingNodeIds.insert(resultGraph.node(i).id());
	}
#endif
	assert(resultGraph.node_size() == vggraph.node_size());
	assert(resultGraph.edge_size() == vggraph.edge_size());
	return resultGraph;
}

void outputGraph(std::string filename, const vg::Graph& graph)
{
	std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
	std::vector<vg::Graph> writeVector {graph};
	stream::write_buffered(alignmentOut, writeVector, 0);
}

std::vector<int> getSinkNodes(const vg::Graph& graph)
{
	std::set<int> notSinkNodes;
	for (int i = 0; i < graph.edge_size(); i++)
	{
		notSinkNodes.insert(graph.edge(i).from());
	}
	std::vector<int> result;
	for (int i = 0; i < graph.node_size(); i++)
	{
		if (notSinkNodes.count(graph.node(i).id()) == 0) result.push_back(graph.node(i).id());
	}
	return result;
}

std::vector<int> getSourceNodes(const vg::Graph& graph)
{
	std::set<int> notSourceNodes;
	for (int i = 0; i < graph.edge_size(); i++)
	{
		notSourceNodes.insert(graph.edge(i).to());
	}
	std::vector<int> result;
	for (int i = 0; i < graph.node_size(); i++)
	{
		if (notSourceNodes.count(graph.node(i).id()) == 0) result.push_back(graph.node(i).id());
	}
	return result;
}

void runComponentMappings(const vg::Graph& graph, const std::vector<const FastQ*>& fastQs, const std::map<const FastQ*, std::vector<vg::Alignment>>& seedhits, std::vector<vg::Alignment>& alignments, int threadnum)
{
	BufferedWriter cerroutput {std::cerr};
	BufferedWriter coutoutput {std::cout};
	for (size_t i = 0; i < fastQs.size(); i++)
	{
		const FastQ* fastq = fastQs[i];
		cerroutput << "thread " << threadnum << " " << i << "/" << fastQs.size() << "\n";
		cerroutput << "read size " << fastq->sequence.size() << "bp" << "\n";
		cerroutput << "components: " << seedhits.at(fastq).size() << BufferedWriter::Flush;
		std::vector<std::tuple<int, int, bool, vg::Graph>> components;
		for (size_t j = 0; j < seedhits.at(fastq).size(); j++)
		{
			auto seedGraphUnordered = ExtractSubgraph(graph, seedhits.at(fastq)[j], fastq->sequence.size());
			int startpos = 0;
			int endpos = 0;
			int score = 0;
			bool forwards;
			auto seedGraph = OrderByFeedbackVertexset(seedGraphUnordered);
			cerroutput << "thread " << threadnum << " read " << i << " component " << j << "/" << seedhits.at(fastq).size() << "\n";
			cerroutput << "component size " << GraphSizeInBp(seedGraph) << "bp" << "\n";
			coutoutput << "out of order before sorting: " << numberOfVerticesOutOfOrder(seedGraphUnordered) << "\n";
			coutoutput << "out of order after sorting: " << numberOfVerticesOutOfOrder(seedGraph) << BufferedWriter::Flush;
			{
				GraphAligner<uint32_t, int32_t> forwardAlignment;
				for (int i = 0; i < seedGraph.node_size(); i++)
				{
					forwardAlignment.AddNode(seedGraph.node(i).id(), seedGraph.node(i).sequence());
				}
				for (int i = 0; i < seedGraph.edge_size(); i++)
				{
					forwardAlignment.AddEdge(seedGraph.edge(i).from(), seedGraph.edge(i).to());
				}
				forwardAlignment.Finalize();
				auto forward = forwardAlignment.GetLocalAlignmentSequencePosition(fastq->sequence);
				forwards = true;
				startpos = std::get<1>(forward);
				endpos = std::get<2>(forward);
				score = std::get<0>(forward);
			}
			{
				GraphAligner<uint32_t, int32_t> backwardAlignment;
				for (int i = seedGraph.node_size()-1; i >= 0; i--)
				{
					backwardAlignment.AddNode(seedGraph.node(i).id(), FastQ::reverseComplement(seedGraph.node(i).sequence()));
				}
				for (int i = 0; i < seedGraph.edge_size(); i++)
				{
					backwardAlignment.AddEdge(seedGraph.edge(i).to(), seedGraph.edge(i).from());
				}
				backwardAlignment.Finalize();
				auto backwards = backwardAlignment.GetLocalAlignmentSequencePosition(fastq->sequence);
				if (std::get<0>(backwards) > score)
				{
					forwards = false;
					startpos = std::get<1>(backwards);
					endpos = std::get<2>(backwards);
					score = std::get<0>(backwards);
				}
			}
			components.emplace_back(startpos, endpos, forwards, seedGraph);
		}
		std::sort(components.begin(), components.end(), [](auto& left, auto& right) { return std::get<0>(left) < std::get<0>(right); });

		GraphAligner<uint32_t, int32_t> forwardAlignment;
		std::vector<std::vector<int>> sources;
		std::vector<std::vector<int>> sinks;
		for (int i = 0; i < components.size(); i++)
		{
			auto& g = std::get<3>(components[i]);
			for (int j = 0; j < g.node_size(); j++)
			{
				forwardAlignment.AddNode(g.node(j).id(), g.node(j).sequence());
			}
			for (int j = 0; j < g.edge_size(); j++)
			{
				forwardAlignment.AddEdge(g.edge(j).from(), g.edge(j).to());
			}
			sources.emplace_back(getSourceNodes(g));
			sinks.emplace_back(getSinkNodes(g));
		}
		for (size_t i = 0; i < components.size(); i++)
		{
			for (size_t j = i+1; j < components.size(); j++)
			{
				std::vector<int> nodeIdsFromFirst;
				if (std::get<2>(components[i]))
				{
					nodeIdsFromFirst = sinks[i];
				}
				else
				{
					nodeIdsFromFirst = sources[i];
				}
				std::vector<int> nodeIdsFromSecond;
				if (std::get<2>(components[j]))
				{
					nodeIdsFromSecond = sources[j];
				}
				else
				{
					nodeIdsFromSecond = sinks[j];
				}
				for (auto firstId : nodeIdsFromFirst)
				{
					for (auto secondId : nodeIdsFromSecond)
					{
						forwardAlignment.AddEdge(firstId, secondId);
					}
				}
			}
		}
		forwardAlignment.Finalize();

		cerroutput << "thread " << threadnum << " augmented graph is " << forwardAlignment.SizeInBp() << "bp" << BufferedWriter::Flush;

		auto alignment = forwardAlignment.AlignOneWay(fastQs[i]->seq_id, fastQs[i]->sequence, false);
		auto reverseComplement = fastQs[i]->reverseComplement();
		auto backwardsalignment = forwardAlignment.AlignOneWay(reverseComplement.seq_id, reverseComplement.sequence, true);
		if (alignment.score() > backwardsalignment.score())
		{
			alignments.push_back(alignment);
		}
		else
		{
			alignments.push_back(backwardsalignment);
		}
		cerroutput << "thread " << threadnum << " successfully aligned read " << fastq->seq_id << BufferedWriter::Flush;
		std::vector<vg::Alignment> alignmentvec;
		alignmentvec.emplace_back(alignments.back());
		std::string filename;
		filename = "alignment_";
		filename += std::to_string(threadnum);
		filename += "_";
		filename += fastq->seq_id;
		filename += ".gam";
		std::replace(filename.begin(), filename.end(), '/', '_');
		cerroutput << "write to " << filename << BufferedWriter::Flush;
		std::ofstream alignmentOut { filename, std::ios::out | std::ios::binary };
		stream::write_buffered(alignmentOut, alignmentvec, 0);
		cerroutput << "write finished" << BufferedWriter::Flush;
	}
	cerroutput << "thread " << threadnum << " finished with " << alignments.size() << " alignments" << BufferedWriter::Flush;
}

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	vg::Graph graph;
	{
		std::cout << "load graph from " << argv[1] << std::endl;
		std::ifstream graphfile { argv[1], std::ios::in | std::ios::binary };
		std::vector<vg::Graph> parts;
		std::function<void(vg::Graph&)> lambda = [&parts](vg::Graph& g) {
			parts.push_back(g);
		};
		stream::for_each(graphfile, lambda);
		graph = mergeGraphs(parts);
		std::cout << "graph is " << GraphSizeInBp(graph) << " bp large" << std::endl;
	}

	auto fastqs = loadFastqFromFile(argv[2]);
	std::cout << fastqs.size() << " reads" << std::endl;

	std::map<std::string, std::vector<vg::Alignment>> seeds;
	{
		std::ifstream seedfile { argv[3], std::ios::in | std::ios::binary };
		std::function<void(vg::Alignment&)> alignmentLambda = [&seeds](vg::Alignment& a) {
			seeds[a.name()].push_back(a);
		};
		stream::for_each(seedfile, alignmentLambda);
	}

	int numThreads = std::stoi(argv[5]);
	std::vector<std::vector<const FastQ*>> readsPerThread;
	std::map<const FastQ*, std::vector<vg::Alignment>> seedHits;
	std::vector<std::vector<vg::Alignment>> resultsPerThread;
	readsPerThread.resize(numThreads);
	resultsPerThread.resize(numThreads);
	int currentThread = 0;
	for (size_t i = 0; i < fastqs.size(); i++)
	{
		if (seeds.count(fastqs[i].seq_id) == 0) continue;
		seedHits[&(fastqs[i])].insert(seedHits[&(fastqs[i])].end(), seeds[fastqs[i].seq_id].begin(), seeds[fastqs[i].seq_id].end());
		readsPerThread[currentThread].push_back(&(fastqs[i]));
		currentThread++;
		currentThread %= numThreads;
	}

	std::vector<std::thread> threads;

	for (int i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&graph, &readsPerThread, &seedHits, &resultsPerThread, i]() { runComponentMappings(graph, readsPerThread[i], seedHits, resultsPerThread[i], i); });
	}

	for (int i = 0; i < numThreads; i++)
	{
		threads[i].join();
	}

	std::vector<vg::Alignment> alignments;

	for (size_t i = 0; i < numThreads; i++)
	{
		alignments.insert(alignments.end(), resultsPerThread[i].begin(), resultsPerThread[i].end());
	}

	std::cerr << "final result has " << alignments.size() << " alignments" << std::endl;

	std::ofstream alignmentOut { argv[4], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, alignments, 0);

	return 0;
}
