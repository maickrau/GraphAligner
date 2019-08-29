#include <omp.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <csignal>
#include "Aligner.h"
#include "stream.hpp"
#include "ThreadReadAssertion.h"

int main(int argc, char** argv)
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	std::cout << "GraphAligner " << VERSION << std::endl;
	std::cerr << "GraphAligner " << VERSION << std::endl;

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
		("reads,f", boost::program_options::value<std::vector<std::string>>()->multitoken(), "input reads (fasta or fastq, uncompressed or gzipped)")
		("alignments-out,a", boost::program_options::value<std::string>(), "output alignment file (.gam/.json)")
		("corrected-out", boost::program_options::value<std::string>(), "output corrected reads file (.fa/.fa.gz)")
		("corrected-clipped-out", boost::program_options::value<std::string>(), "output corrected clipped reads file (.fa/.fa.gz)")
	;
	boost::program_options::options_description general("General parameters");
	general.add_options()
		("help,h", "help message")
		("version", "print version")
		("threads,t", boost::program_options::value<size_t>(), "number of threads (int) (default 1)")
		("verbose", "print progress messages")
		("all-alignments", "return all alignments instead of the best non-overlapping alignments")
		("try-all-seeds", "extend all seeds instead of a reasonable looking subset")
		("global-alignment", "force the read to be aligned end-to-end even if the alignment score is poor")
	;
	boost::program_options::options_description seeding("Seeding");
	seeding.add_options()
		("seeds-minimizer-count", boost::program_options::value<size_t>(), "arg least common minimizers per chunk fully contained in a node (int) (-1 for all)")
		("seeds-minimizer-length", boost::program_options::value<size_t>(), "k-mer length for minimizer seeding (int)")
		("seeds-minimizer-windowsize", boost::program_options::value<size_t>(), "window size for minimizer seeding (int)")
		("seeds-minimizer-chunksize", boost::program_options::value<size_t>(), "chunk size for minimizer seeding (int)")
		("seeds-mum-count", boost::program_options::value<size_t>(), "arg longest maximal unique matches fully contained in a node (int) (-1 for all)")
		("seeds-mem-count", boost::program_options::value<size_t>(), "arg longest maximal exact matches fully contained in a node (int) (-1 for all)")
		("seeds-mxm-length", boost::program_options::value<size_t>(), "minimum length for maximal unique / exact matches (int)")
		("seeds-mxm-cache-prefix", boost::program_options::value<std::string>(), "store the mum/mem seeding index to the disk for reuse, or reuse it if it exists (filename prefix)")
		("seeds-file,s", boost::program_options::value<std::vector<std::string>>()->multitoken(), "external seeds (.gam)")
		("seeds-first-full-rows", boost::program_options::value<int>(), "no seeding, instead calculate the first arg rows fully. VERY SLOW except on tiny graphs (int)")
	;
	boost::program_options::options_description alignment("Extension");
	alignment.add_options()
		("bandwidth,b", boost::program_options::value<size_t>(), "alignment bandwidth (int)")
		("ramp-bandwidth,B", boost::program_options::value<size_t>(), "ramp bandwidth (int)")
		("tangle-effort,C", boost::program_options::value<size_t>(), "tangle effort limit, higher results in slower but more accurate alignments (int) (-1 for unlimited)")
		("high-memory", "use slightly less CPU but a lot more memory")
	;
	boost::program_options::options_description hidden("hidden");
	hidden.add_options()
		("precise-clipping", "clip the alignment ends more precisely. Recommended for Illumina reads")
	;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(mandatory).add(general).add(seeding).add(alignment).add(hidden);

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
		std::cerr << mandatory << std::endl << general << std::endl << seeding;
		std::cerr << "defaults are --seeds-minimizer-count 5 --seeds-minimizer-length 19 --seeds-minimizer-windowsize 30 --seeds-minimizer-chunksize 100" << std::endl << std::endl;
		std::cerr << alignment;
		std::cerr << "defaults are -b 5 -B 10 -C 10000" << std::endl << std::endl;
		std::exit(0);
	}
	if (vm.count("version"))
	{
		std::cout << "Version " << VERSION << std::endl;
		std::exit(0);
	}

	AlignerParams params;
	params.graphFile = "";
	params.outputAlignmentFile = "";
	params.outputCorrectedFile = "";
	params.outputCorrectedClippedFile = "";
	params.numThreads = 1;
	params.initialBandwidth = 0;
	params.rampBandwidth = 0;
	params.dynamicRowStart = 0;
	params.maxCellsPerSlice = std::numeric_limits<decltype(params.maxCellsPerSlice)>::max();
	params.verboseMode = false;
	params.tryAllSeeds = false;
	params.highMemory = false;
	params.mxmLength = 20;
	params.mumCount = 0;
	params.memCount = 0;
	params.seederCachePrefix = "";
	params.outputAllAlns = false;
	params.forceGlobal = false;
	params.outputJSON = false;
	params.compressCorrected = false;
	params.compressClipped = false;
	params.preciseClipping = false;
	params.minimizerCount = 0;
	params.minimizerLength = 19;
	params.minimizerWindowSize = 30;
	params.minimizerChunkSize = 100;

	if (vm.count("graph")) params.graphFile = vm["graph"].as<std::string>();
	if (vm.count("reads")) params.fastqFiles = vm["reads"].as<std::vector<std::string>>();
	if (vm.count("alignments-out")) params.outputAlignmentFile = vm["alignments-out"].as<std::string>();
	if (vm.count("corrected-out")) params.outputCorrectedFile = vm["corrected-out"].as<std::string>();
	if (vm.count("corrected-clipped-out")) params.outputCorrectedClippedFile = vm["corrected-clipped-out"].as<std::string>();
	if (vm.count("threads")) params.numThreads = vm["threads"].as<size_t>();
	if (vm.count("bandwidth")) params.initialBandwidth = vm["bandwidth"].as<size_t>();

	if (vm.count("seeds-minimizer-count")) params.minimizerCount = vm["seeds-minimizer-count"].as<size_t>();
	if (vm.count("seeds-minimizer-length")) params.minimizerLength = vm["seeds-minimizer-length"].as<size_t>();
	if (vm.count("seeds-minimizer-windowsize")) params.minimizerWindowSize = vm["seeds-minimizer-windowsize"].as<size_t>();
	if (vm.count("seeds-minimizer-chunksize")) params.minimizerChunkSize = vm["seeds-minimizer-chunksize"].as<size_t>();
	if (vm.count("seeds-file")) params.seedFiles = vm["seeds-file"].as<std::vector<std::string>>();
	if (vm.count("seeds-mxm-length")) params.mxmLength = vm["seeds-mxm-length"].as<size_t>();
	if (vm.count("seeds-mem-count")) params.memCount = vm["seeds-mem-count"].as<size_t>();
	if (vm.count("seeds-mum-count")) params.mumCount = vm["seeds-mum-count"].as<size_t>();
	if (vm.count("seeds-mxm-cache-prefix")) params.seederCachePrefix = vm["seeds-mxm-cache-prefix"].as<std::string>();
	if (vm.count("seeds-first-full-rows")) params.dynamicRowStart = vm["seeds-first-full-rows"].as<int>();

	if (vm.count("ramp-bandwidth")) params.rampBandwidth = vm["ramp-bandwidth"].as<size_t>();
	if (vm.count("tangle-effort")) params.maxCellsPerSlice = vm["tangle-effort"].as<size_t>();
	if (vm.count("all-alignments"))
	{
		params.outputAllAlns = true;
		params.tryAllSeeds = true;
	}
	if (vm.count("verbose")) params.verboseMode = true;
	if (vm.count("try-all-seeds")) params.tryAllSeeds = true;
	if (vm.count("high-memory")) params.highMemory = true;
	if (vm.count("global-alignment")) params.forceGlobal = true;
	if (vm.count("precise-clipping")) params.preciseClipping = true;

	bool paramError = false;

	if (params.graphFile == "")
	{
		std::cerr << "graph file must be given" << std::endl;
		paramError = true;
	}
	if (params.fastqFiles.size() == 0)
	{
		std::cerr << "read file must be given" << std::endl;
		paramError = true;
	}
	if (params.outputAlignmentFile == "" && params.outputCorrectedFile == "" && params.outputCorrectedClippedFile == "")
	{
		std::cerr << "one of alignments-out, corrected-out or corrected-clipped-out must be given" << std::endl;
		paramError = true;
	}
	if (params.outputAlignmentFile != "" && params.outputAlignmentFile.substr(params.outputAlignmentFile.size()-4) != ".gam" && params.outputAlignmentFile.substr(params.outputAlignmentFile.size()-5) != ".json")
	{
		std::cerr << "unknown output alignment format, must be either .gam or .json" << std::endl;
		paramError = true;
	}
	if (params.outputCorrectedFile != "" && params.outputCorrectedFile.substr(params.outputCorrectedFile.size()-3) != ".fa" && params.outputCorrectedFile.substr(params.outputCorrectedFile.size()-6) != ".fasta" && params.outputCorrectedFile.substr(params.outputCorrectedFile.size()-6) != ".fa.gz" && params.outputCorrectedFile.substr(params.outputCorrectedFile.size()-9) != ".fasta.gz")
	{
		std::cerr << "unknown output corrected read format, must be .fa or .fa.gz" << std::endl;
		paramError = true;
	}
	if (params.outputCorrectedClippedFile != "" && params.outputCorrectedClippedFile.substr(params.outputCorrectedClippedFile.size()-3) != ".fa" && params.outputCorrectedClippedFile.substr(params.outputCorrectedClippedFile.size()-6) != ".fasta" && params.outputCorrectedClippedFile.substr(params.outputCorrectedClippedFile.size()-6) != ".fa.gz" && params.outputCorrectedClippedFile.substr(params.outputCorrectedClippedFile.size()-9) != ".fasta.gz")
	{
		std::cerr << "unknown output corrected read format, must be .fa or .fa.gz" << std::endl;
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
	if (params.initialBandwidth == 0 && params.rampBandwidth == 0 && params.maxCellsPerSlice == std::numeric_limits<decltype(params.maxCellsPerSlice)>::max())
	{
		//default extension parameters
		params.initialBandwidth = 5;
		params.rampBandwidth = 10;
		params.maxCellsPerSlice = 10000;
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
	if (params.mxmLength < 2)
	{
		std::cerr << "mum/mem minimum length must be >= 2" << std::endl;
		paramError = true;
	}
	if (params.minimizerLength >= sizeof(size_t)*8/2)
	{
		std::cerr << "Maximum minimizer length is " << (sizeof(size_t)*8/2)-1 << std::endl;
		paramError = true;
	}
	int pickedSeedingMethods = ((params.dynamicRowStart != 0) ? 1 : 0) + ((params.seedFiles.size() > 0) ? 1 : 0) + ((params.mumCount != 0) ? 1 : 0) + ((params.memCount != 0) ? 1 : 0) + ((params.minimizerCount != 0) ? 1 : 0);
	if (pickedSeedingMethods == 0)
	{
		//use minimizers as the default seeding method
		params.minimizerCount = 5;
		params.minimizerLength = 19;
		params.minimizerWindowSize = 30;
		params.minimizerChunkSize = 100;
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

	if (params.outputAlignmentFile.size() >= 5 && params.outputAlignmentFile.substr(params.outputAlignmentFile.size()-5) == ".json")
	{
		params.outputJSON = true;
	}

	if (params.outputCorrectedFile.size() >= 3 && params.outputCorrectedFile.substr(params.outputCorrectedFile.size()-3) == ".gz")
	{
		params.compressCorrected = true;
	}

	if (params.outputCorrectedClippedFile.size() >= 3 && params.outputCorrectedClippedFile.substr(params.outputCorrectedClippedFile.size()-3) == ".gz")
	{
		params.compressClipped = true;
	}

	omp_set_num_threads(params.numThreads);

	alignReads(params);

	return 0;
}
