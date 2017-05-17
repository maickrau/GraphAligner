#include <vector>
#include <map>
#include <cassert>
#include "vg.pb.h"
#include "TopologicalSort.h"

using namespace std;
void topological_sort_using_DFS_loop(vector<vector<size_t> >& graph, vector<size_t>& sorted);
void topological_sort_using_DFS(vector<vector<size_t> >& graph, vector<bool>& explored, size_t i, vector<size_t>& sorted, size_t& t);

std::vector<size_t> topologicalSort(const vg::Graph& vggraph)
{
	std::vector<std::vector<size_t>> graph;
	graph.resize(vggraph.node_size());

	std::map<int, int> ids;
	for (int i = 0; i < vggraph.node_size(); i++)
	{
		ids[vggraph.node(i).id()] = i;
	}
	for (int i = 0; i < vggraph.edge_size(); i++)
	{
		assert(ids[vggraph.edge(i).from()] < graph.size());
		graph[ids[vggraph.edge(i).from()]].push_back(ids[vggraph.edge(i).to()]);
	}

	vector<size_t> sorted(vggraph.node_size(), 0);
	topological_sort_using_DFS_loop(graph, sorted);
	return sorted;
}

std::vector<size_t> topologicalSort(const DirectedGraph& digraph)
{
	std::vector<std::vector<size_t>> graph;
	graph.resize(digraph.nodes.size());

	for (size_t i = 0; i < digraph.edges.size(); i++)
	{
		graph[digraph.edges[i].fromIndex].push_back(digraph.edges[i].toIndex);
	}

	vector<size_t> sorted(digraph.nodes.size(), 0);
	topological_sort_using_DFS_loop(graph, sorted);
	return sorted;
}



void topological_sort_using_DFS(vector<vector<size_t> >& graph, vector<bool>& explored, size_t i, vector<size_t>& sorted, size_t& t)
{
	explored[i] = true;

	for (size_t j = 0; j < graph[i].size(); ++j)
	{
		if (explored[graph[i][j]] == false)
		{
			topological_sort_using_DFS(graph, explored, graph[i][j], sorted, t);
		}
	}

	--t;
	sorted[t] = i;

	return;
}

void topological_sort_using_DFS_loop(vector<vector<size_t> >& graph, vector<size_t>& sorted)
{
	vector<bool> explored(graph.size(), false);
	size_t t = graph.size();

	for (size_t i = 0; i < graph.size(); ++i)
	{
		if (explored[i] == false)
		{
			topological_sort_using_DFS(graph, explored, i, sorted, t);
		}
	}

	assert(t == 0);

	return;
}
