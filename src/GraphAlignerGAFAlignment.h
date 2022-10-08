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
		MatchOrMismatch,
		Insertion,
		Deletion,
		Empty
	};
public:

	static std::string traceToAlignment(const std::string& seq_id, const std::string& sequence, const GraphAlignerCommon<size_t, int64_t, uint64_t>::OnewayTrace& tracePair, double alignmentXScore, int mappingQuality, const Params& params, bool cigarMatchMismatchMerge, const bool includecigar)
	{
		auto& trace = tracePair.trace;
		if (trace.size() == 0) return nullptr;
		std::stringstream cigar;
		std::string readName = seq_id;
		size_t readLen = sequence.size();
		size_t readStart = trace[0].DPposition.seqPos;
		size_t readEnd = trace.back().DPposition.seqPos+1;
		bool strand = true;
		std::stringstream nodePath;
		size_t nodePathLen = 0;
		size_t nodePathStart = trace[0].DPposition.nodeOffset;
		size_t nodePathEnd = 0;
		size_t matches = 0;
		size_t blockLength = trace.size();

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
		if (cigarMatchMismatchMerge)
		{
			currentEdit = MatchOrMismatch;
			editLength = 1;
			if (Common::characterMatch(trace[0].sequenceCharacter, trace[0].graphCharacter))
			{
				matches += 1;
			}
			else
			{
				mismatches += 1;
			}
		}
		else if (Common::characterMatch(trace[0].sequenceCharacter, trace[0].graphCharacter))
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
		addPosToString(nodePath, currentPos, params);
		nodePathLen += params.graph.BigraphNodeSize(currentPos.nodeId);
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
				size_t skippedBefore = params.graph.BigraphNodeSize(currentPos.nodeId) - 1 - trace[pos-1].DPposition.nodeOffset;
				currentPos = newPos;
				addPosToString(nodePath, currentPos, params);
				assert(trace[pos].DPposition.nodeOffset < params.graph.BigraphNodeSize(currentPos.nodeId));
				size_t skippedAfter = trace[pos].DPposition.nodeOffset;
				nodePathLen += params.graph.BigraphNodeSize(currentPos.nodeId) - (skippedBefore + skippedAfter);
			}

			if (trace[pos-1].DPposition.seqPos == trace[pos].DPposition.seqPos)
			{
				if (currentEdit == Empty) currentEdit = Deletion;
				if (currentEdit != Deletion)
				{
					if (includecigar) addCigarItem(cigar, editLength, currentEdit);
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
					if (includecigar) addCigarItem(cigar, editLength, currentEdit);
					currentEdit = Insertion;
					editLength = 0;
				}
				editLength += 1;
				insertions += 1;
			}
			else if (cigarMatchMismatchMerge)
			{
				if (currentEdit == Empty) currentEdit = MatchOrMismatch;
				if (currentEdit != MatchOrMismatch)
				{
					if (includecigar) addCigarItem(cigar, editLength, currentEdit);
					currentEdit = MatchOrMismatch;
					editLength = 0;
				}
				editLength += 1;
				if (Common::characterMatch(trace[pos].sequenceCharacter, trace[pos].graphCharacter))
				{
					matches += 1;
				}
				else
				{
					mismatches += 1;
				}
			}
			else if (Common::characterMatch(trace[pos].sequenceCharacter, trace[pos].graphCharacter))
			{
				if (currentEdit == Empty) currentEdit = Match;
				if (currentEdit != Match)
				{
					if (includecigar) addCigarItem(cigar, editLength, currentEdit);
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
					if (includecigar) addCigarItem(cigar, editLength, currentEdit);
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
		if (includecigar) addCigarItem(cigar, editLength, currentEdit);

		nodePathEnd = nodePathLen - (params.graph.BigraphNodeSize(trace.back().DPposition.node) - 1 - trace.back().DPposition.nodeOffset);

		std::stringstream sstr;
		sstr << readName << "\t" << readLen << "\t" << readStart << "\t" << readEnd << "\t" << (strand ? "+" : "-") << "\t" << nodePath.str() << "\t" << nodePathLen << "\t" << nodePathStart << "\t" << nodePathEnd << "\t" << matches << "\t" << blockLength << "\t" << mappingQuality;
		sstr << "\t" << "NM:i:" << (mismatches + deletions + insertions);
		if (alignmentXScore != -1) sstr << "\t" << "AS:f:" << alignmentXScore;
		sstr << "\t" << "dv:f:" << 1.0-((double)matches / (double)(matches + mismatches + deletions + insertions));
		sstr << "\t" << "id:f:" << ((double)matches / (double)(matches + mismatches + deletions + insertions));
		if (includecigar) sstr << "\t" << "cg:Z:" << cigar.str();
		return sstr.str();
	}

private:

	static void addPosToString(std::stringstream& str, MergedNodePos pos, const Params& params)
	{
		if (pos.reverse)
		{
			str << "<";
		}
		else
		{
			str << ">";
		}
		std::string nodeName = params.graph.BigraphNodeName(pos.nodeId);
		str << nodeName;
	}

	static void addCigarItem(std::stringstream& str, size_t editLength, EditType type)
	{
		if (editLength == 0) return;
		str << editLength;
		switch(type)
		{
			case MatchOrMismatch:
				str << "M";
				break;
			case Match:
				str << "=";
				break;
			case Mismatch:
				str << "X";
				break;
			case Insertion:
				str << "I";
				break;
			case Deletion:
				str << "D";
				break;
			case Empty:
			default:
				return;
		}
	}
};

#endif