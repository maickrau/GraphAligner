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
	using TraceItem = typename Common::TraceItem;
	struct MergedNodePos
	{
		int nodeId;
		bool reverse;
		size_t nodeOffset;
		size_t seqPos;
	};
	enum EditType
	{
		Match,
		Mismatch,
		Insertion,
		Deletion,
		Empty
	};
public:

	static std::shared_ptr<vg::Alignment> traceToAlignment(const std::string& seq_id, const std::string& sequence, ScoreType score, const std::vector<TraceItem>& trace, size_t mappingQuality, size_t cellsProcessed, bool reverse)
	{
		if (trace.size() == 0) return nullptr;
		vg::Alignment* aln = new vg::Alignment;
		std::shared_ptr<vg::Alignment> result { aln };
		result->set_name(seq_id);
		result->set_score(score);
		result->set_sequence(sequence);
		result->set_mapping_quality(mappingQuality);
		auto path = new vg::Path;
		result->set_allocated_path(path);
		MergedNodePos currentPos;
		currentPos.nodeId = trace[0].DPposition.node;
		currentPos.reverse = (trace[0].DPposition.node % 2) == 1;
		currentPos.nodeOffset = trace[0].DPposition.nodeOffset;
		currentPos.seqPos = trace[0].DPposition.seqPos;
		int rank = 0;
		auto vgmapping = path->add_mapping();
		auto position = new vg::Position;
		vgmapping->set_allocated_position(position);
		vgmapping->set_rank(rank);
		auto edit = vgmapping->add_edit();
		EditType currentEdit = Empty;
		size_t mismatches = 0;
		size_t deletions = 0;
		size_t insertions = 0;
		size_t matches = 0;
		if (Common::characterMatch(trace[0].sequenceCharacter, trace[0].graphCharacter))
		{
			currentEdit = Match;
			edit->set_from_length(edit->from_length()+1);
			edit->set_to_length(edit->to_length()+1);
			matches += 1;
		}
		else
		{
			currentEdit = Mismatch;
			edit->set_from_length(edit->from_length()+1);
			edit->set_to_length(edit->to_length()+1);
			edit->set_sequence(std::string { sequence[0] });
			mismatches += 1;
		}
		position->set_node_id(currentPos.nodeId);
		position->set_is_reverse(currentPos.reverse);
		position->set_offset(currentPos.nodeOffset);
		for (size_t pos = 1; pos < trace.size(); pos++)
		{
			assert(trace[pos].DPposition.seqPos < sequence.size());
			MergedNodePos newPos;
			newPos.nodeId = trace[pos].DPposition.node;
			newPos.reverse = (trace[pos].DPposition.node % 2) == 1;
			newPos.nodeOffset = trace[pos].DPposition.nodeOffset;
			newPos.seqPos = trace[pos].DPposition.seqPos;
			bool insideNode = !trace[pos-1].nodeSwitch || (newPos.nodeId == currentPos.nodeId && newPos.reverse == currentPos.reverse && newPos.nodeOffset > currentPos.nodeOffset);

			assert(newPos.seqPos >= currentPos.seqPos);

			if (!insideNode)
			{
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
				currentEdit = Empty;
			}

			if (trace[pos-1].DPposition.seqPos == trace[pos].DPposition.seqPos)
			{
				if (currentEdit == Empty) currentEdit = Deletion;
				if (currentEdit != Deletion)
				{
					edit = vgmapping->add_edit();
					currentEdit = Deletion;
				}
				edit->set_from_length(edit->from_length()+1);
				deletions += 1;
			}
			else if (insideNode && trace[pos-1].DPposition.nodeOffset == trace[pos].DPposition.nodeOffset)
			{
				if (currentEdit == Empty) currentEdit = Insertion;
				if (currentEdit != Insertion)
				{
					edit = vgmapping->add_edit();
					currentEdit = Insertion;
				}
				edit->set_to_length(edit->to_length()+1);
				edit->set_sequence(edit->sequence() + trace[pos].sequenceCharacter);
				insertions += 1;
			}
			else if (Common::characterMatch(trace[pos].sequenceCharacter, trace[pos].graphCharacter))
			{
				if (currentEdit == Empty) currentEdit = Match;
				if (currentEdit != Match)
				{
					edit = vgmapping->add_edit();
					currentEdit = Match;
				}
				edit->set_from_length(edit->from_length()+1);
				edit->set_to_length(edit->to_length()+1);
				matches += 1;
			}
			else
			{
				if (currentEdit == Empty) currentEdit = Mismatch;
				if (currentEdit != Mismatch)
				{
					edit = vgmapping->add_edit();
					currentEdit = Mismatch;
				}
				edit->set_from_length(edit->from_length()+1);
				edit->set_to_length(edit->to_length()+1);
				edit->set_sequence(edit->sequence() + trace[pos].sequenceCharacter);
				mismatches += 1;
			}
			if (insideNode)
			{
				assert(trace[pos-1].nodeSwitch || newPos.nodeId == currentPos.nodeId);
				assert(trace[pos-1].nodeSwitch || newPos.reverse == currentPos.reverse);
			}
		}
		result->set_identity((double)matches / (double)(matches + mismatches + insertions + deletions));
		assert(currentEdit != Empty);
		return result;
	}

	static bool posEqual(const vg::Position& pos1, const vg::Position& pos2)
	{
		return pos1.node_id() == pos2.node_id() && pos1.is_reverse() == pos2.is_reverse();
	}
};

#endif