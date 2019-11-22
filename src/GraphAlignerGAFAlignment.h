#ifndef GraphAlignerGAFAlignment_h
#define GraphAlignerGAFAlignment_h

#include <string>
#include <vector>
#include "AlignmentGraph.h"
#include "NodeSlice.h"
#include "CommonUtils.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerCommon.h"

template <typename LengthType, typename ScoreType, typename Word>
class GraphAlignerGAFAlignment
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

	static std::string traceToAlignment(const std::string& seq_id, const std::string& sequence, const GraphAlignerCommon<size_t, int32_t, uint64_t>::OnewayTrace& tracePair, const Params& params)
	{
		auto& trace = tracePair.trace;
		if (trace.size() == 0) return nullptr;
		std::string cigar;
		std::string readName = seq_id;
		size_t readLen = sequence.size();
		size_t readStart = trace[0].DPposition.seqPos;
		size_t readEnd = trace.back().DPposition.seqPos+1;
		bool strand = true;
		std::string nodePath;
		size_t nodePathLen = 0;
		size_t nodePathStart = trace[0].DPposition.nodeOffset;
		size_t nodePathEnd = 0;
		size_t matches = 0;
		size_t blockLength = trace.size();
		int mappingQuality = 255;

		MergedNodePos currentPos;
		currentPos.nodeId = trace[0].DPposition.node;
		currentPos.reverse = (trace[0].DPposition.node % 2) == 1;
		currentPos.nodeOffset = trace[0].DPposition.nodeOffset;
		currentPos.seqPos = trace[0].DPposition.seqPos;
		EditType currentEdit = Empty;
		size_t mismatches = 0;
		size_t deletions = 0;
		size_t insertions = 0;
		size_t editLength = 0;
		if (Common::characterMatch(trace[0].sequenceCharacter, trace[0].graphCharacter))
		{
			currentEdit = Match;
			editLength = 1;
			matches += 1;
		}
		else
		{
			currentEdit = Mismatch;
			editLength = 1;
			mismatches += 1;
		}
		nodePath += posToString(currentPos, params);
		nodePathLen += params.graph.originalNodeSize.at(currentPos.nodeId);
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
				size_t skippedBefore = params.graph.originalNodeSize.at(currentPos.nodeId) - 1 - trace[pos-1].DPposition.nodeOffset;
				currentPos = newPos;
				nodePath += posToString(currentPos, params);
				assert(trace[pos].DPposition.nodeOffset < params.graph.originalNodeSize.at(currentPos.nodeId));
				size_t skippedAfter = trace[pos].DPposition.nodeOffset;
				nodePathLen += params.graph.originalNodeSize.at(currentPos.nodeId) - (skippedBefore + skippedAfter);
			}

			if (trace[pos-1].DPposition.seqPos == trace[pos].DPposition.seqPos)
			{
				if (currentEdit == Empty) currentEdit = Deletion;
				if (currentEdit != Deletion)
				{
					cigar += cigarItem(editLength, currentEdit);
					currentEdit = Deletion;
					editLength = 0;
				}
				editLength += 1;
				deletions += 1;
			}
			else if (insideNode && trace[pos-1].DPposition.nodeOffset == trace[pos].DPposition.nodeOffset)
			{
				if (currentEdit == Empty) currentEdit = Insertion;
				if (currentEdit != Insertion)
				{
					cigar += cigarItem(editLength, currentEdit);
					currentEdit = Insertion;
					editLength = 0;
				}
				editLength += 1;
				insertions += 1;
			}
			else if (Common::characterMatch(trace[pos].sequenceCharacter, trace[pos].graphCharacter))
			{
				if (currentEdit == Empty) currentEdit = Match;
				if (currentEdit != Match)
				{
					cigar += cigarItem(editLength, currentEdit);
					currentEdit = Match;
					editLength = 0;
				}
				editLength += 1;
				matches += 1;
			}
			else
			{
				if (currentEdit == Empty) currentEdit = Mismatch;
				if (currentEdit != Mismatch)
				{
					cigar += cigarItem(editLength, currentEdit);
					currentEdit = Mismatch;
					editLength = 0;
				}
				editLength += 1;
				mismatches += 1;
			}
			if (insideNode)
			{
				assert(trace[pos-1].nodeSwitch || newPos.nodeId == currentPos.nodeId);
				assert(trace[pos-1].nodeSwitch || newPos.reverse == currentPos.reverse);
			}
		}

		assert(matches + mismatches + deletions + insertions == trace.size());
		cigar += cigarItem(editLength, currentEdit);

		nodePathEnd = nodePathLen - (params.graph.originalNodeSize.at(trace.back().DPposition.node) - 1 - trace.back().DPposition.nodeOffset);

		std::stringstream sstr;
		sstr << readName << "\t" << readLen << "\t" << readStart << "\t" << readEnd << "\t" << (strand ? "+" : "-") << "\t" << nodePath << "\t" << nodePathLen << "\t" << nodePathStart << "\t" << nodePathEnd << "\t" << matches << "\t" << blockLength << "\t" << mappingQuality << "\t" << "cg:Z:" << cigar;
		return sstr.str();
	}

private:

	static std::string posToString(MergedNodePos pos, const Params& params)
	{
		std::string result;
		if (pos.reverse)
		{
			result += "<";
		}
		else
		{
			result += ">";
		}
		std::string nodeName = params.graph.originalNodeName.at(pos.nodeId);
		if (nodeName == "") nodeName = std::to_string(pos.nodeId/2);
		result += nodeName;
		return result;
	}

	static std::string cigarItem(size_t editLength, EditType type)
	{
		if (editLength == 0) return "";
		std::string result = std::to_string(editLength);
		switch(type)
		{
			case Match:
				result += "M";
				break;
			case Mismatch:
				result += "M";
				break;
			case Insertion:
				result += "I";
				break;
			case Deletion:
				result += "D";
				break;
			case Empty:
			default:
				return "";
		}
		return result;
	}
};

#endif