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
	struct MergedNodePos
	{
		int nodeId;
		bool reverse;
		size_t nodeOffset;
		size_t seqPos;
	};
public:

	static AlignmentResult::AlignmentItem traceToAlignment(const Params& params, const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<std::pair<MatrixPosition, bool>>& trace, size_t cellsProcessed, bool reverse)
	{
		vg::Alignment* aln = new vg::Alignment;
		std::shared_ptr<vg::Alignment> result { aln };
		result->set_name(seq_id);
		result->set_score(score);
		result->set_sequence(sequence);
		auto path = new vg::Path;
		result->set_allocated_path(path);
		if (trace.size() == 0) return emptyAlignment(0, cellsProcessed);
		size_t currentPosIndex = 0;
		MergedNodePos currentPos;
		currentPos.nodeId = trace[0].first.node;
		currentPos.reverse = (trace[0].first.node % 2) == 1;
		currentPos.nodeOffset = trace[0].first.nodeOffset;
		currentPos.seqPos = trace[0].first.seqPos;
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		auto edit = vgmapping->add_edit();
		position->set_node_id(currentPos.nodeId);
		position->set_is_reverse(currentPos.reverse);
		position->set_offset(currentPos.nodeOffset);
		MergedNodePos btNodeStart = currentPos;
		MergedNodePos btNodeEnd = currentPos;
		MergedNodePos btBeforeNode = currentPos;
		btBeforeNode.seqPos--;
		for (size_t pos = 1; pos < trace.size(); pos++)
		{
			assert(trace[pos].first.seqPos < sequence.size());
			MergedNodePos newPos;
			newPos.nodeId = trace[pos].first.node;
			newPos.reverse = (trace[pos].first.node % 2) == 1;
			newPos.nodeOffset = trace[pos].first.nodeOffset;
			newPos.seqPos = trace[pos].first.seqPos;
			if (!trace[pos-1].second)
			{
				assert(newPos.nodeId == currentPos.nodeId);
				assert(newPos.reverse == currentPos.reverse);
				btNodeEnd = newPos;
				continue;
			}
			if (newPos.nodeId == currentPos.nodeId && newPos.reverse == currentPos.reverse && newPos.nodeOffset > currentPos.nodeOffset)
			{
				btNodeEnd = newPos;
				continue;
			}

			assert(newPos.seqPos >= currentPos.seqPos);
			assert(btNodeEnd.nodeId == btNodeStart.nodeId);
			assert(btNodeEnd.seqPos >= btNodeStart.seqPos);
			assert(btNodeEnd.nodeOffset >= btNodeStart.nodeOffset);
			assert(btNodeEnd.seqPos >= btBeforeNode.seqPos || btBeforeNode.seqPos == -1);

			edit->set_to_length(edit->to_length() + btNodeEnd.seqPos - btBeforeNode.seqPos);
			if (btNodeEnd.seqPos > btBeforeNode.seqPos || btBeforeNode.seqPos == -1)
			{
				assert(btBeforeNode.seqPos < sequence.size() - 1 || btBeforeNode.seqPos == -1);
				edit->set_sequence(edit->sequence() + sequence.substr(btBeforeNode.seqPos+1, btNodeEnd.seqPos - btBeforeNode.seqPos));
			}
			assert(btNodeEnd.nodeOffset + 1 >= btNodeStart.nodeOffset);
			edit->set_from_length(edit->from_length() + btNodeEnd.nodeOffset - btNodeStart.nodeOffset + 1);

			rank++;
			currentPos = newPos;
			vgmapping = path->add_mapping();
			position = new vg::Position;
			vgmapping->set_allocated_position(position);
			vgmapping->set_rank(rank);
			position->set_offset(currentPos.nodeOffset);
			position->set_node_id(currentPos.nodeId);
			position->set_is_reverse(currentPos.reverse);
			edit = vgmapping->add_edit();

			btBeforeNode = btNodeEnd;
			btNodeStart = newPos;
			btNodeEnd = newPos;
			currentPosIndex = pos;
		}
		assert(btNodeEnd.seqPos >= btBeforeNode.seqPos || btBeforeNode.seqPos == -1);
		edit->set_to_length(edit->to_length() + btNodeEnd.seqPos - btBeforeNode.seqPos);
		if (btBeforeNode.seqPos != btNodeEnd.seqPos)
		{
			assert(btBeforeNode.seqPos < sequence.size() - 1 || btBeforeNode.seqPos == -1);
			edit->set_sequence(edit->sequence() + sequence.substr(btBeforeNode.seqPos+1, btNodeEnd.seqPos - btBeforeNode.seqPos));
		}
		assert(btNodeEnd.nodeOffset >= btNodeStart.nodeOffset);
		edit->set_from_length(edit->from_length() + btNodeEnd.nodeOffset - btNodeStart.nodeOffset + 1);
		AlignmentResult::AlignmentItem item { result, cellsProcessed, std::numeric_limits<size_t>::max() };
		item.alignmentStart = trace[0].first.seqPos;
		item.alignmentEnd = trace.back().first.seqPos;
		return item;
	}

	static AlignmentResult::AlignmentItem emptyAlignment(size_t elapsedMilliseconds, size_t cellsProcessed)
	{
		vg::Alignment* aln = new vg::Alignment;
		std::shared_ptr<vg::Alignment> result { aln };
		result->set_score(std::numeric_limits<decltype(result->score())>::max());
		return AlignmentResult::AlignmentItem { result, cellsProcessed, elapsedMilliseconds };
	}

	static bool posEqual(const vg::Position& pos1, const vg::Position& pos2)
	{
		return pos1.node_id() == pos2.node_id() && pos1.is_reverse() == pos2.is_reverse();
	}
};

#endif