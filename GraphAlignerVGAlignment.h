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
	using SeedHit = typename Common::SeedHit;
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
		vg::Alignment result;
		result.set_name(seq_id);
		result.set_score(score);
		result.set_sequence(sequence);
		auto path = new vg::Path;
		result.set_allocated_path(path);
		if (trace.size() == 0) return emptyAlignment(0, cellsProcessed);
		size_t currentPosIndex = 0;
		MergedNodePos currentPos;
		currentPos.nodeId = params.graph.nodeIDs[trace[0].first.node];
		currentPos.reverse = params.graph.reverse[trace[0].first.node];
		currentPos.nodeOffset = trace[0].first.nodeOffset + params.graph.nodeOffset[trace[0].first.node];
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
		for (size_t pos = 1; pos < trace.size(); pos++)
		{
			assert(trace[pos].first.seqPos < sequence.size());
			MergedNodePos newPos;
			newPos.nodeId = params.graph.nodeIDs[trace[pos].first.node];
			newPos.reverse = params.graph.reverse[trace[pos].first.node];
			newPos.nodeOffset = trace[pos].first.nodeOffset + params.graph.nodeOffset[trace[pos].first.node];
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

			assert(btNodeEnd.nodeOffset == params.graph.DBGOverlap || btNodeEnd.nodeOffset == params.graph.nodeOffset[params.graph.nodeLookup.at(btNodeEnd.nodeId).back()] + params.graph.NodeLength(params.graph.nodeLookup.at(btNodeEnd.nodeId).back()) - params.graph.DBGOverlap - 1 || btNodeEnd.nodeOffset == params.graph.nodeOffset[params.graph.nodeLookup.at(btNodeEnd.nodeId).back()] + params.graph.NodeLength(params.graph.nodeLookup.at(btNodeEnd.nodeId).back()) - 1);
			assert(currentPosIndex == 0 || btNodeStart.nodeOffset == 0 || btNodeStart.nodeOffset == params.graph.DBGOverlap || btNodeStart.nodeOffset == params.graph.nodeOffset[params.graph.nodeLookup.at(btNodeStart.nodeId).back()] + params.graph.NodeLength(params.graph.nodeLookup.at(btNodeStart.nodeId).back()) - params.graph.DBGOverlap - 1);
			assert(newPos.seqPos >= currentPos.seqPos);
			assert(btNodeEnd.nodeId == btNodeStart.nodeId);
			assert(btNodeEnd.seqPos >= btNodeStart.seqPos);
			assert(btNodeEnd.nodeOffset >= btNodeStart.nodeOffset);
			assert(btNodeEnd.seqPos >= btBeforeNode.seqPos);

			edit->set_to_length(edit->to_length() + btNodeEnd.seqPos - btBeforeNode.seqPos);
			if (btNodeEnd.seqPos > btBeforeNode.seqPos)
			{
				assert(btBeforeNode.seqPos < sequence.size() - 1);
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
		assert(btNodeEnd.seqPos >= btBeforeNode.seqPos);
		edit->set_to_length(edit->to_length() + btNodeEnd.seqPos - btBeforeNode.seqPos);
		if (btBeforeNode.seqPos != btNodeEnd.seqPos)
		{
			assert(btBeforeNode.seqPos < sequence.size() - 1);
			edit->set_sequence(edit->sequence() + sequence.substr(btBeforeNode.seqPos+1, btNodeEnd.seqPos - btBeforeNode.seqPos));
		}
		assert(btNodeEnd.nodeOffset >= btNodeStart.nodeOffset);
		edit->set_from_length(edit->from_length() + btNodeEnd.nodeOffset - btNodeStart.nodeOffset + 1);
		AlignmentResult::AlignmentItem item { result, cellsProcessed, std::numeric_limits<size_t>::max() };
		item.alignmentStart = trace[0].first.seqPos;
		item.alignmentEnd = trace.back().first.seqPos;
		return item;
	}

	// static AlignmentResult::AlignmentItem mergeAlignments(const Params& params, const AlignmentResult::AlignmentItem& first, const SeedHit& seedHit, const AlignmentResult::AlignmentItem& second, const std::string& sequence)
	// {
	// 	assert(!first.alignmentFailed() || !second.alignmentFailed());
	// 	int seedHitNodeId;
	// 	char seedHitChar = sequence[seedHit.seqPos];
	// 	if (seedHit.reverse)
	// 	{
	// 		seedHitNodeId = seedHit.nodeID * 2 + 1;
	// 	}
	// 	else
	// 	{
	// 		seedHitNodeId = seedHit.nodeID * 2;
	// 	}
	// 	if (first.alignmentFailed() || first.alignment.path().mapping_size() == 0)
	// 	{
	// 		//todo add middle
	// 		return second;
	// 	}
	// 	if (second.alignmentFailed() || second.alignment.path().mapping_size() == 0)
	// 	{
	// 		//todo add middle
	// 		return first;
	// 	}
	// 	assert(!first.alignmentFailed());
	// 	assert(!second.alignmentFailed());
	// 	AlignmentResult::AlignmentItem finalResult;
	// 	finalResult.cellsProcessed = first.cellsProcessed + second.cellsProcessed;
	// 	finalResult.elapsedMilliseconds = first.elapsedMilliseconds + second.elapsedMilliseconds;
	// 	finalResult.alignment = first.alignment;
	// 	finalResult.alignment.set_score(first.alignment.score() + second.alignment.score());
	// 	int start = 0;
	// 	auto firstEndPos = first.alignment.path().mapping(first.alignment.path().mapping_size()-1).position();
	// 	auto secondStartPos = second.alignment.path().mapping(0).position();
	// 	assert(secondStartPos.node_id() == seedHitNodeId);
	// 	assert(secondStartPos.is_reverse() == seedHit.reverse);
	// 	vg::Mapping* mapping;

	// 	if (seedHit.nodeOffset != params.graph.DBGOverlap)
	// 	{
	// 		assert(firstEndPos.node_id() == seedHitNodeId);
	// 		assert(firstEndPos.is_reverse() == seedHit.reverse);
	// 		mapping = finalResult.alignment.mutable_path()->mutable_mapping(finalResult.alignment.path().mapping_size()-1);
	// 	}
	// 	else if (firstEndPos.node_id() == seedHitNodeId && firstEndPos.is_reverse() == seedHit.reverse)
	// 	{
	// 		mapping = finalResult.alignment.mutable_path()->mutable_mapping(finalResult.alignment.path().mapping_size()-1);
	// 	}
	// 	else
	// 	{
	// 		mapping = finalResult.alignment.mutable_path()->add_mapping();
	// 		mapping->set_node_id(seedHitNodeId);
	// 		mapping->set_is_reverse(seedHit.reverse);
	// 		mapping->set_offset(seedHit.nodeOffset);
	// 	}

	// 	auto middlePart = mapping->add_edit();
	// 	middlePart->set_sequence(seedHitChar);
	// 	middlePart->set_from_length(1);
	// 	middlePart->set_to_length(1);

	// 	assert(seedHit.nodeOffset <= params.graph.nodeOffset[params.graph.nodeLookup[seedHitNodeId].back()] + params.graph.NodeLength(params.graph.nodeLookup[seedHitNodeId].back()) - 1)
	// 	if (seedHit.nodeOffset < params.graph.nodeOffset[params.graph.nodeLookup[seedHitNodeId].back()] + params.graph.NodeLength(params.graph.nodeLookup[seedHitNodeId].back()) - 1)
	// 	{
	// 		assert(second.alignment.path().mapping(0).node_id() == mapping->node_id());
	// 		assert(second.alignment.path().mapping(0).is_reverse() == mapping->is_reverse());
	// 	}
	// 	else if ()
	// 	{
	// 	}
	// 	else
	// 	{
	// 		mapping = finalResult.alignment.mutable_path()->add_mapping();
	// 		mapping->set_node_id(second.alignment.path().mapping(0).node_id());
	// 		mapping->set_is_reverse(second.alignment.path().mapping(0).is_reverse());
	// 		mapping->set_offset(second.alignment.path().mapping(0).offset());
	// 	}

	// 	auto secondStartEdit = mapping->add_edit();
	// 	secondStartEdit->set_sequence(second.alignment.path().mapping(0).edit(0).sequence());
	// 	secondStartEdit->set_from_length(second.alignment.path().mapping(0).edit(0).from_length());
	// 	secondStartEdit->set_to_length(second.alignment.path().mapping(0).edit(0).to_length());

	// 	for (int i = 1; i < second.alignment.path().mapping_size(); i++)
	// 	{
	// 		auto mapping = finalResult.alignment.mutable_path()->add_mapping();
	// 		*mapping = second.alignment.path().mapping(i);
	// 	}
	// 	return finalResult;
	// }

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