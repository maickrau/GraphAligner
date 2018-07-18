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
	using Common = GraphAlignerCommon<LengthType, ScoreType, Word>;
	using Params = typename Common::Params;
	using MatrixPosition = typename Common::MatrixPosition;
public:

	static AlignmentResult::AlignmentItem traceToAlignment(const Params& params, const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<MatrixPosition>& trace, size_t cellsProcessed, bool reverse)
	{
		vg::Alignment result;
		result.set_name(seq_id);
		result.set_score(score);
		result.set_sequence(sequence);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		if (trace.size() == 0) return emptyAlignment(0, cellsProcessed);
		size_t pos = 0;
		size_t oldNode = trace[0].node;
		int rank = 0;
		int oldNodeId = params.graph.nodeIDs[oldNode];
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		auto edit = vgmapping->add_edit();
		position->set_node_id(params.graph.nodeIDs[oldNode]);
		if (reverse)
		{
			position->set_is_reverse(!params.graph.reverse[oldNode]);
		}
		else
		{
			position->set_is_reverse(params.graph.reverse[oldNode]);
		}
		position->set_offset(trace[pos].nodeOffset);
		MatrixPosition btNodeStart = trace[pos];
		MatrixPosition btNodeEnd = trace[pos];
		MatrixPosition btBeforeNode = trace[pos];
		for (; pos < trace.size(); pos++)
		{
			if (trace[pos].node == oldNode)
			{
				btNodeEnd = trace[pos];
				continue;
			}
			assert(trace[pos].seqPos >= trace[pos-1].seqPos);
			assert(btNodeEnd.node == btNodeStart.node);
			assert(btNodeEnd.seqPos >= btNodeStart.seqPos);
			assert(reverse || btNodeEnd.nodeOffset >= btNodeStart.nodeOffset);
			assert(!reverse || btNodeEnd.nodeOffset <= btNodeStart.nodeOffset);
			auto previousNode = oldNode;
			oldNode = trace[pos].node;
			edit->set_from_length(edit->from_length() + btNodeEnd.nodeOffset - btNodeStart.nodeOffset + 1);
			edit->set_to_length(edit->to_length() + btNodeEnd.seqPos - btBeforeNode.seqPos);
			if (reverse)
			{
				edit->set_sequence(edit->sequence() + sequence.substr(btNodeStart.seqPos, btBeforeNode.seqPos - btNodeEnd.seqPos));
			}
			else
			{
				edit->set_sequence(edit->sequence() + sequence.substr(btNodeStart.seqPos, btNodeEnd.seqPos - btBeforeNode.seqPos));
			}
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
				if (reverse)
				{
					position->set_is_reverse(!params.graph.reverse[oldNode]);
				}
				else
				{
					position->set_is_reverse(params.graph.reverse[oldNode]);
				}
				edit = vgmapping->add_edit();
			}
		}
		edit->set_from_length(edit->from_length() + btNodeEnd.nodeOffset - btNodeStart.nodeOffset);
		edit->set_to_length(edit->to_length() + btNodeEnd.seqPos - btBeforeNode.seqPos);
		if (reverse)
		{
			edit->set_sequence(edit->sequence() + sequence.substr(btNodeStart.seqPos, btBeforeNode.seqPos - btNodeEnd.seqPos));
		}
		else
		{
			edit->set_sequence(edit->sequence() + sequence.substr(btNodeStart.seqPos, btNodeEnd.seqPos - btBeforeNode.seqPos));
		}
		AlignmentResult::AlignmentItem item { result, cellsProcessed, std::numeric_limits<size_t>::max() };
		item.alignmentStart = trace[0].seqPos;
		item.alignmentEnd = trace.back().seqPos;
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