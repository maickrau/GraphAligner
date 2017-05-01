#include <iostream>
#include "vg.pb.h"

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;
	vg::Graph graph;
	if (!graph.ParseFromIstream(&std::cin)) {
		std::cerr << "Graph load failed." << std::endl;
		return 1;
	}
	if (!graph.SerializeToOstream(&std::cout)) {
		std::cerr << "Graph write failed." << std::endl;
		return 1;
	}
	return 0;
}
