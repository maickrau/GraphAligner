#include <fstream>
#include "fastqloader.h"
#include "vg.pb.h"
#include "stream.hpp"

int main(int argc, char** argv)
{
	auto fastqs = loadFastqFromFile(argv[1]);
	std::ifstream gremfile {argv[2], std::ios::in};
	std::vector<vg::Alignment> output;
	do
	{
		std::string line;
		std::getline(gremfile, line);
		if (line == "") break;
		int commapos = line.find(',', 0);
		int nodeid = std::stoi(line.substr(0, commapos));
		int readid = std::stoi(line.substr(commapos+1));
		std::string readname = fastqs[readid].seq_id;
		vg::Alignment alignment;
		vg::Path* path = new vg::Path();
		alignment.set_name(readname);
		alignment.set_allocated_path(path);
		auto mapping = path->add_mapping();
		vg::Position* pos = new vg::Position();
		mapping->set_allocated_position(pos);
		pos->set_node_id(nodeid);
		output.push_back(alignment);
	} while (gremfile.good());

	std::ofstream alignmentOut { argv[3], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, output, 0);
}
