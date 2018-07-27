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



	std::vector<vg::Alignment> writeAlns;
	std::ifstream alnfile { inputAlignmentFile, std::ios::in | std::ios::binary };
	std::ofstream outfile { outputAlignmentFile,  std::ios::out | std::ios::binary };
	std::function<void(vg::Alignment&)> lambda = [&writeAlns, &outfile, &mapping](vg::Alignment& g) {
		for (int i = 0; i < g.path().mapping_size(); i++)
		{
			auto position = g.mutable_path()->mutable_mapping(i)->mutable_position();
			assert(position->node_id() > 0);
			assert(position->node_id() < mapping.size());
			position->set_name(mapping[position->node_id()]);
		}
		writeAlns.push_back(g);
		if (writeAlns.size() > 10000)
		{
			stream::write_buffered(outfile, writeAlns, 0);
			writeAlns.clear();
		}
	};
	stream::for_each(alnfile, lambda);
}
