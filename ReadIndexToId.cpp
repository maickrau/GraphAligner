#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include "fastqloader.h"
#include "vg.pb.h"
#include "stream.hpp"

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

int main(int argc, char** argv)
{
	auto fastqs = loadFastqFromFile(argv[1]);
	std::ifstream gremfile {argv[2], std::ios::in};
	std::vector<vg::Alignment> output;
	std::set<std::pair<int, int>> existing;
	do
	{
		std::string line;
		std::getline(gremfile, line);
		if (line == "") break;
		std::istringstream iss {line};
		std::vector<std::string> parts = split(line, ',');
		// std::cerr << parts[0] << std::endl;
		// std::cerr << parts[1] << std::endl;
		// std::cerr << parts[2] << std::endl;
		int nodeid = std::stoi(parts[0]);
		int readid = std::stoi(parts[1]);
		int readpos = std::stoi(parts[2]);
		if (existing.count(std::make_pair(nodeid, readid)) > 0) continue;
		existing.insert(std::make_pair(nodeid, readid));
		std::string readname = fastqs[readid].seq_id;
		vg::Alignment alignment;
		vg::Path* path = new vg::Path();
		alignment.set_name(readname);
		alignment.set_query_position(readpos);
		alignment.set_allocated_path(path);
		auto mapping = path->add_mapping();
		vg::Position* pos = new vg::Position();
		pos->set_node_id(nodeid);
		mapping->set_allocated_position(pos);
		output.push_back(alignment);
	} while (gremfile.good());

	std::ofstream alignmentOut { argv[3], std::ios::out | std::ios::binary };
	stream::write_buffered(alignmentOut, output, 0);
}
