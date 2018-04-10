#include <string>
#include <fstream>
#include "CommonUtils.h"
#include "stream.hpp"

int main(int argc, char** argv)
{
	std::string inputMappingFile { argv[1] };
	std::string inputAlignmentFile { argv[2] };
	std::string outputAlignmentFile { argv[3] };

	std::vector<std::string> mapping;
	mapping.push_back("");
	{
		std::ifstream file {inputMappingFile};
		while (file.good())
		{
			std::string name;
			getline(file, name);
			if (file.good()) mapping.push_back(name);
		}
	}

	auto alns = CommonUtils::LoadVGAlignments(inputAlignmentFile);
	for (auto& aln : alns)
	{
		for (int i = 0; i < aln.path().mapping_size(); i++)
		{
			auto position = aln.mutable_path()->mutable_mapping(i)->mutable_position();
			assert(position->node_id() > 0);
			assert(position->node_id() < mapping.size());
			position->set_name(mapping[position->node_id()]);
		}
	}

	std::ofstream outfile { outputAlignmentFile,  std::ios::out | std::ios::binary };
	stream::write_buffered(outfile, alns, 0);
}
