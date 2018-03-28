#ifndef GraphAlignerVGAlignment_h
#define GraphAlignerVGAlignment_h

#include <string>
#include <vector>
#include "AlignmentGraph.h"
#include "vg.pb.h"
#include "NodeSlice.h"
#include "CommonUtils.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerCommon.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerVGAlignment
{
	typedef GraphAlignerParams<LengthType, ScoreType, Word> Params;
	typedef std::pair<LengthType, LengthType> MatrixPosition;
public:
	static AlignmentResult::AlignmentItem traceToAlignment(const Params& params, const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, size_t cellsProcessed)
	{
		vg::Alignment result;
		result.set_name(seq_id);
		result.set_score(score);
		result.set_sequence(sequence);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		if (trace.size() == 0) return emptyAlignment(0, cellsProcessed);
		size_t pos = 0;
		size_t oldNode = params.graph.IndexToNode(trace[0].first);
		while (oldNode == params.graph.dummyNodeStart)
		{
			pos++;
			if (pos == trace.size()) return emptyAlignment(0, cellsProcessed);
			assert(pos < trace.size());
			assert(trace[pos].second >= trace[pos-1].second);
			oldNode = params.graph.IndexToNode(trace[pos].first);
			assert(oldNode < params.graph.nodeIDs.size());
		}
		if (oldNode == params.graph.dummyNodeEnd) return emptyAlignment(0, cellsProcessed);
		int rank = 0;
		int oldNodeId = params.graph.nodeIDs[oldNode];
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		auto edit = vgmapping->add_edit();
		position->set_node_id(params.graph.nodeIDs[oldNode]);
		position->set_is_reverse(params.graph.reverse[oldNode]);
		position->set_offset(trace[pos].first - params.graph.NodeStart(oldNode));
		MatrixPosition btNodeStart = trace[pos];
		MatrixPosition btNodeEnd = trace[pos];
		MatrixPosition btBeforeNode = trace[pos];
		for (; pos < trace.size(); pos++)
		{
			if (params.graph.IndexToNode(trace[pos].first) == params.graph.dummyNodeEnd) break;
			if (params.graph.IndexToNode(trace[pos].first) == oldNode)
			{
				btNodeEnd = trace[pos];
				continue;
			}
			assert(trace[pos].second >= trace[pos-1].second);
			assert(params.graph.IndexToNode(btNodeEnd.first) == params.graph.IndexToNode(btNodeStart.first));
			assert(btNodeEnd.second >= btNodeStart.second);
			assert(btNodeEnd.first >= btNodeStart.first);
			auto previousNode = oldNode;
			oldNode = params.graph.IndexToNode(trace[pos].first);
			edit->set_from_length(edit->from_length() + btNodeEnd.first - btNodeStart.first + 1);
			edit->set_to_length(edit->to_length() + btNodeEnd.second - btBeforeNode.second);
			edit->set_sequence(edit->sequence() + sequence.substr(btNodeStart.second, btNodeEnd.second - btBeforeNode.second));
			btBeforeNode = btNodeEnd;
			btNodeStart = trace[pos];
			btNodeEnd = trace[pos];
			if (params.graph.nodeIDs[oldNode] != oldNodeId || params.graph.reverse[oldNode] != params.graph.reverse[previousNode] || params.graph.nodeOffset[oldNode] != params.graph.nodeOffset[previousNode] + params.graph.SPLIT_NODE_SIZE)
			{
				rank++;
				oldNodeId = params.graph.nodeIDs[oldNode];
				vgmapping = path->add_mapping();
				position = new vg::Position;
				vgmapping->set_allocated_position(position);
				vgmapping->set_rank(rank);
				position->set_offset(params.graph.nodeOffset[oldNode]);
				position->set_node_id(params.graph.nodeIDs[oldNode]);
				position->set_is_reverse(params.graph.reverse[oldNode]);
				edit = vgmapping->add_edit();
			}
		}
		edit->set_from_length(edit->from_length() + btNodeEnd.first - btNodeStart.first);
		edit->set_to_length(edit->to_length() + btNodeEnd.second - btBeforeNode.second);
		edit->set_sequence(edit->sequence() + sequence.substr(btNodeStart.second, btNodeEnd.second - btBeforeNode.second));
		AlignmentResult::AlignmentItem item { result, cellsProcessed, std::numeric_limits<size_t>::max() };
		item.alignmentStart = trace[0].second;
		item.alignmentEnd = trace.back().second;
		return item;
	}

	static AlignmentResult::AlignmentItem mergeAlignments(const Params& params, const AlignmentResult::AlignmentItem& first, const AlignmentResult::AlignmentItem& second)
	{
		assert(!first.alignmentFailed() || !second.alignmentFailed());
		if (first.alignmentFailed()) return second;
		if (second.alignmentFailed()) return first;
		if (first.alignment.path().mapping_size() == 0) return second;
		if (second.alignment.path().mapping_size() == 0) return first;
		assert(!first.alignmentFailed());
		assert(!second.alignmentFailed());
		AlignmentResult::AlignmentItem finalResult;
		finalResult.cellsProcessed = first.cellsProcessed + second.cellsProcessed;
		finalResult.elapsedMilliseconds = first.elapsedMilliseconds + second.elapsedMilliseconds;
		finalResult.alignment = first.alignment;
		finalResult.alignment.set_score(first.alignment.score() + second.alignment.score());
		int start = 0;
		auto firstEndPos = first.alignment.path().mapping(first.alignment.path().mapping_size()-1).position();
		auto secondStartPos = second.alignment.path().mapping(0).position();
		auto firstEndPosNodeId = params.graph.nodeLookup.at(firstEndPos.node_id()).back();
		auto secondStartPosNodeId = params.graph.nodeLookup.at(secondStartPos.node_id())[0];
		if (posEqual(firstEndPos, secondStartPos))
		{
			start = 1;
		}
		else if (std::find(params.graph.outNeighbors[firstEndPosNodeId].begin(), params.graph.outNeighbors[firstEndPosNodeId].end(), secondStartPosNodeId) != params.graph.outNeighbors[firstEndPosNodeId].end())
		{
			start = 0;
		}
		for (int i = start; i < second.alignment.path().mapping_size(); i++)
		{
			auto mapping = finalResult.alignment.mutable_path()->add_mapping();
			*mapping = second.alignment.path().mapping(i);
		}
		return finalResult;
	}

	static AlignmentResult::AlignmentItem emptyAlignment(size_t elapsedMilliseconds, size_t cellsProcessed)
	{
		vg::Alignment result;
		result.set_score(std::numeric_limits<decltype(result.score())>::max());
		return AlignmentResult::AlignmentItem { result, cellsProcessed, elapsedMilliseconds };
	}

	static bool posEqual(const vg::Position& pos1, const vg::Position& pos2)
	{
		return pos1.node_id() == pos2.node_id() && pos1.is_reverse() == pos2.is_reverse();
	}
};

#endif