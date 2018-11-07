#include <boost/program_options.hpp>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <limits>
#include "Aligner.h"
#include "stream.hpp"
#include "ThreadReadAssertion.h"
#include "ByteStuff.h"

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	std::cout << "Branch " << GITBRANCH << " commit " << GITCOMMIT << " " << GITDATE << std::endl;
	std::cerr << "Branch " << GITBRANCH << " commit " << GITCOMMIT << " " << GITDATE << std::endl;

#ifndef NOBUILTINPOPCOUNT
	if (__builtin_cpu_supports("popcnt") == 0)
	{
		std::cerr << "CPU does not support builtin popcount operation" << std::endl;
		std::cerr << "recompile with -DNOBUILTINPOPCOUNT" << std::endl;
		std::abort();
	}
#endif

	struct sigaction act;
	act.sa_handler = ThreadReadAssertion::signal;
	sigemptyset(&act.sa_mask);
	act.sa_flags = 0;
	sigaction(SIGSEGV, &act, 0);

	boost::program_options::options_description mandatory("Mandatory parameters");
	mandatory.add_options()
		("graph,g", boost::program_options::value<std::string>(), "input graph (.gfa / .vg)")
		("reads,f", boost::program_options::value<std::string>(), "input reads (fasta or fastq, uncompressed or gzipped)")
		("bandwidth,b", boost::program_options::value<int>(), "default bandwidth (int)")
		("threads,t", boost::program_options::value<int>(), "number of threads (int)")
		("alignments-out,a", boost::program_options::value<std::string>(), "output alignment file (.gam)")
		;
	boost::program_options::options_description seeding("Seeding");
	seeding.add_options()
		("seeds-file,s", boost::program_options::value<std::string>(), "external seeds (.gam)")
		("seeds-mum", boost::program_options::value<int>(), "n longest maximal unique matches of fully contained in a node (int)")
		("seeds-mum-cache-prefix", boost::program_options::value<std::string>(), "store the mum seeding index to the disk for reuse, or reuse it if it exists (filename prefix)")
		("seeds-first-full-rows", boost::program_options::value<int>(), "no seeding, instead calculate the first arg rows fully. VERY SLOW except on tiny graphs (int)")
	;
	boost::program_options::options_description optional("Optional parameters");
	optional.add_options()
		("help,h", "help message")
		("ramp-bandwidth,B", boost::program_options::value<int>(), "ramp bandwidth (int)")
		("tangle-effort,C", boost::program_options::value<int>(), "tangle effort limit, higher results in slower but more accurate alignments, default is unlimited (int)")
		("max-alignments,A", boost::program_options::value<int>(), "return up to arg non-overlapping alignments instead of all alignments (int)")
		("quiet,q", "don't print progress messages")
		("sloppy-optimizations,u", "use speed-up heuristics which might result in missing alignments")
		("low-memory,l", "use less memory with a higher runtime")
	;
	boost::program_options::options_description hidden("don't use these unless you know what you're doing");
	hidden.add_options()
		("subgraph-extraction-heuristic,S", "")
	;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(mandatory).add(seeding).add(optional).add(hidden);

	boost::program_options::variables_map vm;
	try
	{
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdline_options), vm);
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << "run with option -h for help" << std::endl;
		std::exit(1);
	}
	boost::program_options::notify(vm);

	if (vm.count("help"))
	{
		std::cerr << mandatory << std::endl << seeding << std::endl << optional << std::endl;
		std::exit(0);
	}

	AlignerParams params;
	params.graphFile = "";
	params.fastqFile = "";
	params.seedFile = "";
	params.outputAlignmentFile = "";
	params.numThreads = 0;
	params.initialBandwidth = 0;
	params.rampBandwidth = 0;
	params.dynamicRowStart = 0;
	params.maxCellsPerSlice = std::numeric_limits<decltype(params.maxCellsPerSlice)>::max();
	params.quietMode = false;
	params.sloppyOptimizations = false;
	params.lowMemory = false;
	params.maxAlns = 0;
	params.useSubgraph = false;
	params.mumCount = 0;
	params.seederCachePrefix = "";

	if (vm.count("graph")) params.graphFile = vm["graph"].as<std::string>();
	if (vm.count("reads")) params.fastqFile = vm["reads"].as<std::string>();
	if (vm.count("alignments-out")) params.outputAlignmentFile = vm["alignments-out"].as<std::string>();
	if (vm.count("threads")) params.numThreads = vm["threads"].as<int>();
	if (vm.count("bandwidth")) params.initialBandwidth = vm["bandwidth"].as<int>();

	if (vm.count("seeds-file")) params.seedFile = vm["seeds-file"].as<std::string>();
	if (vm.count("seeds-mum")) params.mumCount = vm["seeds-mum"].as<int>();
	if (vm.count("seeds-mum-cache-prefix")) params.seederCachePrefix = vm["seeds-mum-cache-prefix"].as<std::string>();
	if (vm.count("seeds-first-full-rows")) params.dynamicRowStart = vm["seeds-first-full-rows"].as<int>();

	if (vm.count("ramp-bandwidth")) params.rampBandwidth = vm["ramp-bandwidth"].as<int>();
	if (vm.count("tangle-effort")) params.maxCellsPerSlice = vm["tangle-effort"].as<int>();
	if (vm.count("max-alignments")) params.maxAlns = vm["max-alignments"].as<int>();
	if (vm.count("quiet")) params.quietMode = true;
	if (vm.count("sloppy-optimizations")) params.sloppyOptimizations = true;
	if (vm.count("low-memory")) params.lowMemory = true;

	if (vm.count("subgraph-extraction-heuristic")) params.useSubgraph = true;

	bool paramError = false;

	if (params.graphFile == "")
	{
		std::cerr << "graph file must be given" << std::endl;
		paramError = true;
	}
	if (params.fastqFile == "")
	{
		std::cerr << "read file must be given" << std::endl;
		paramError = true;
	}
	if (params.outputAlignmentFile == "")
	{
		std::cerr << "alignments-out must be given" << std::endl;
		paramError = true;
	}
	if (params.dynamicRowStart % 64 != 0)
	{
		std::cerr << "first-full-rows has to be a multiple of 64" << std::endl;
		paramError = true;
	}
	if (params.numThreads < 1)
	{
		std::cerr << "number of threads must be >= 1" << std::endl;
		paramError = true;
	}
	if (params.initialBandwidth < 1)
	{
		std::cerr << "default bandwidth must be >= 1" << std::endl;
		paramError = true;
	}
	if (params.rampBandwidth != 0 && params.rampBandwidth <= params.initialBandwidth)
	{
		std::cerr << "ramp bandwidth must be higher than default bandwidth" << std::endl;
		paramError = true;
	}
	int pickedSeedingMethods = ((params.dynamicRowStart != 0) ? 1 : 0) + ((params.seedFile != "") ? 1 : 0) + ((params.mumCount != 0) ? 1 : 0);
	if (pickedSeedingMethods == 0)
	{
		std::cerr << "pick a seeding method" << std::endl;
		paramError = true;
	}
	if (pickedSeedingMethods > 1)
	{
		std::cerr << "pick only one seeding method" << std::endl;
		paramError = true;
	}

	if (paramError)
	{
		std::cerr << "run with option -h for help" << std::endl;
		std::exit(1);
	}

	alignReads(params);

	return 0;
}
