#include <boost/program_options.hpp>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <csignal>
#include "Aligner.h"
#include "stream.hpp"
#include "ThreadReadAssertion.h"
#include "EValue.h"

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
		("alignments-out,a", boost::program_options::value<std::vector<std::string>>(), "output alignment file (.gaf/.gam/.json)")
		("corrected-out", boost::program_options::value<std::string>(), "output corrected reads file (.fa/.fa.gz)")
		("corrected-clipped-out", boost::program_options::value<std::string>(), "output corrected clipped reads file (.fa/.fa.gz)")
	;
	boost::program_options::options_description presets("Preset parameters");
	presets.add_options()
		("preset,x", boost::program_options::value<std::string>(), 
			"Preset parameters\n"
			"dbg - Parameters optimized for de Bruijn graphs\n"
			"vg - Parameters optimized for variation graphs")
	;
	boost::program_options::options_description general("General parameters");
	general.add_options()
		("help,h", "help message")
		("version", "print version")
		("threads,t", boost::program_options::value<size_t>(), "number of threads (int) (default 1)")
		("verbose", "print progress messages")
		("E-cutoff", boost::program_options::value<double>(), "discard alignments with E-value > arg (double)")
		("min-alignment-score", boost::program_options::value<double>(), "discard alignments with alignment score < arg (double) (default 0)")
		("multimap-score-fraction", boost::program_options::value<double>(), "discard alignments whose alignment score is less than this fraction of the best overlapping alignment (double) (default 0.9)")
	;
	boost::program_options::options_description seeding("Seeding");
	seeding.add_options()
		("max-cluster-extend", boost::program_options::value<size_t>(), "extend up to arg seed clusters (int) (-1 for all) (default 10)")
		("seeds-clustersize", boost::program_options::value<size_t>(), "discard seed clusters with fewer than arg seeds (int)")
		("seeds-minimizer-length", boost::program_options::value<size_t>(), "k-mer length for minimizer seeding (int)")
		("seeds-minimizer-windowsize", boost::program_options::value<size_t>(), "window size for minimizer seeding (int)")
		("seeds-minimizer-density", boost::program_options::value<double>(), "keep approximately (arg * sequence length) least frequent minimizers (double) (-1 for all)")
		("seeds-minimizer-ignore-frequent", boost::program_options::value<double>(), "ignore arg most frequent fraction of minimizers (double)")
		("seeds-mum-count", boost::program_options::value<size_t>(), "arg longest maximal unique matches (int) (-1 for all)")
		("seeds-mem-count", boost::program_options::value<size_t>(), "arg longest maximal exact matches (int) (-1 for all)")
		("seeds-mxm-length", boost::program_options::value<size_t>(), "minimum length for maximal unique / exact matches (int)")
		("seeds-mxm-cache-prefix", boost::program_options::value<std::string>(), "store the mum/mem seeding index to the disk for reuse, or reuse it if it exists (filename prefix)")
	;
	boost::program_options::options_description alignment("Extension");
	alignment.add_options()
		("bandwidth,b", boost::program_options::value<size_t>(), "alignment bandwidth (int)")
		("tangle-effort,C", boost::program_options::value<size_t>(), "tangle effort limit (int) (-1 for unlimited)")
		("X-drop", boost::program_options::value<int>(), "X-drop alignment ending score cutoff (int)")
		("precise-clipping", boost::program_options::value<double>(), "clip the alignment ends with arg as the identity cutoff between correct / wrong alignments (double) (default 0.66)")
		("max-trace-count", boost::program_options::value<size_t>(), "backtrace from up to arg highest scoring local maxima per cluster (int) (-1 for all)")
	;
	boost::program_options::options_description hidden("hidden");
	hidden.add_options()
		("cigar-match-mismatch", "use M for matches and mismatches in the cigar string instead of = and X")
		("seeds-file,s", boost::program_options::value<std::vector<std::string>>()->multitoken(), "external seeds (.gam)")
		("seedless-DP", "no seeding, instead use DP alignment algorithm for the entire first row. VERY SLOW except on tiny graphs")
		("DP-restart-stride", boost::program_options::value<size_t>(), "if --seedless-DP doesn't span the entire read, restart after arg base pairs (int)")
		("hpc-collapse-reads", "Collapse homopolymer runs in input reads")
		("discard-cigar", "Don't include CIGAR string in gaf output")
		("clip-ambiguous-ends", boost::program_options::value<int>(), "clip ambiguous alignment ends with alignment score cutoff arg")
		("overlap-incompatible-cutoff", boost::program_options::value<double>(), "consider two partial alignments incompatible if they overlap by arg% of the length of the shorter one")
	;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(mandatory).add(general).add(seeding).add(alignment).add(hidden).add(presets);

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
		std::cerr << mandatory << std::endl << general << std::endl << seeding << std::endl << alignment << std::endl;
		std::cerr << presets << std::endl;
		std::exit(0);
	}
	if (vm.count("version"))
	{
		std::cout << "Version " << VERSION << std::endl;
		std::exit(0);
	}

	AlignerParams params;
	params.graphFile = "";
	params.outputGAMFile = "";
	params.outputJSONFile = "";
	params.outputGAFFile = "";
	params.outputCorrectedFile = "";
	params.outputCorrectedClippedFile = "";
	params.numThreads = 1;
	params.alignmentBandwidth = 0;
	params.dynamicRowStart = false;
	params.maxCellsPerSlice = std::numeric_limits<decltype(params.maxCellsPerSlice)>::max();
	params.verboseMode = false;
	params.mxmLength = 20;
	params.mumCount = 0;
	params.memCount = 0;
	params.seederCachePrefix = "";
	params.selectionECutoff = -1;
	params.compressCorrected = false;
	params.compressClipped = false;
	params.minimizerSeedDensity = 0;
	params.minimizerLength = 19;
	params.minimizerWindowSize = 30;
	params.seedClusterMinSize = 1;
	params.minimizerDiscardMostNumerousFraction = 0.0002;
	params.maxClusterExtend = 5;
	params.preciseClippingIdentityCutoff = 0.66;
	params.Xdropcutoff = 50;
	params.DPRestartStride = 0;
	params.multimapScoreFraction = 0.9;
	params.cigarMatchMismatchMerge = false;
	params.minAlignmentScore = 0;
	params.hpcCollapse = false;
	params.includeCigar = true;
	params.clipAmbiguousEnds = -1;
	params.maxTraceCount = 10;
	params.overlapIncompatibleCutoff = 0.3;

	std::vector<std::string> outputAlns;
	bool paramError = false;

	if (vm.count("preset"))
	{
		std::string preset = vm["preset"].as<std::string>();
		if (preset == "dbg")
		{
			params.minimizerSeedDensity = 5;
			params.minimizerLength = 19;
			params.minimizerWindowSize = 30;
			params.maxClusterExtend = 5;
			params.alignmentBandwidth = 5;
			params.maxCellsPerSlice = 10000;
		}
		else if (preset == "vg")
		{
			params.minimizerSeedDensity = 10;
			params.minimizerLength = 15;
			params.minimizerWindowSize = 20;
			params.maxClusterExtend = 5;
			params.minimizerDiscardMostNumerousFraction = 0.001;
			params.alignmentBandwidth = 10;
		}
		else
		{
			std::cerr << "unknown preset \"" << preset << "\"" << std::endl;
			paramError = true;
		}
	}

	if (vm.count("graph")) params.graphFile = vm["graph"].as<std::string>();
	if (vm.count("reads")) params.fastqFiles = vm["reads"].as<std::vector<std::string>>();
	if (vm.count("alignments-out")) outputAlns = vm["alignments-out"].as<std::vector<std::string>>();
	if (vm.count("corrected-out")) params.outputCorrectedFile = vm["corrected-out"].as<std::string>();
	if (vm.count("corrected-clipped-out")) params.outputCorrectedClippedFile = vm["corrected-clipped-out"].as<std::string>();
	if (vm.count("threads")) params.numThreads = vm["threads"].as<size_t>();
	if (vm.count("bandwidth")) params.alignmentBandwidth = vm["bandwidth"].as<size_t>();

	if (vm.count("max-cluster-extend")) params.maxClusterExtend = vm["max-cluster-extend"].as<size_t>();
	if (vm.count("seeds-minimizer-ignore-frequent")) params.minimizerDiscardMostNumerousFraction = vm["seeds-minimizer-ignore-frequent"].as<double>();
	if (vm.count("seeds-clustersize")) params.seedClusterMinSize = vm["seeds-clustersize"].as<size_t>();
	if (vm.count("seeds-minimizer-density")) params.minimizerSeedDensity = vm["seeds-minimizer-density"].as<double>();
	if (vm.count("seeds-minimizer-length")) params.minimizerLength = vm["seeds-minimizer-length"].as<size_t>();
	if (vm.count("seeds-minimizer-windowsize")) params.minimizerWindowSize = vm["seeds-minimizer-windowsize"].as<size_t>();
	if (vm.count("seeds-file")) params.seedFiles = vm["seeds-file"].as<std::vector<std::string>>();
	if (vm.count("seeds-mxm-length")) params.mxmLength = vm["seeds-mxm-length"].as<size_t>();
	if (vm.count("seeds-mem-count")) params.memCount = vm["seeds-mem-count"].as<size_t>();
	if (vm.count("seeds-mum-count")) params.mumCount = vm["seeds-mum-count"].as<size_t>();
	if (vm.count("seeds-mxm-cache-prefix")) params.seederCachePrefix = vm["seeds-mxm-cache-prefix"].as<std::string>();
	if (vm.count("seedless-DP")) params.dynamicRowStart = true;
	if (vm.count("DP-restart-stride")) params.DPRestartStride = vm["DP-restart-stride"].as<size_t>();
	if (vm.count("multimap-score-fraction")) params.multimapScoreFraction = vm["multimap-score-fraction"].as<double>();

	if (vm.count("tangle-effort")) params.maxCellsPerSlice = vm["tangle-effort"].as<size_t>();
	if (vm.count("verbose")) params.verboseMode = true;
	if (vm.count("cigar-match-mismatch")) params.cigarMatchMismatchMerge = true;
	if (vm.count("min-alignment-score")) params.minAlignmentScore = vm["min-alignment-score"].as<double>();
	if (vm.count("max-trace-count")) params.maxTraceCount = vm["max-trace-count"].as<size_t>();

	if (vm.count("verbose")) params.verboseMode = true;
	if (vm.count("precise-clipping")) params.preciseClippingIdentityCutoff = vm["precise-clipping"].as<double>();
	if (vm.count("hpc-collapse-reads")) params.hpcCollapse = true;
	if (vm.count("discard-cigar")) params.includeCigar = false;
	if (vm.count("clip-ambiguous-ends")) params.clipAmbiguousEnds = vm["clip-ambiguous-ends"].as<int>();
	if (vm.count("overlap-incompatible-cutoff")) params.overlapIncompatibleCutoff = vm["overlap-incompatible-cutoff"].as<double>();

	if (vm.count("X-drop"))
	{
		params.Xdropcutoff = vm["X-drop"].as<int>();
	}
	else
	{
		// by default pick x-drop so that a block of 50 mismatches breaks an alignment
		if (params.preciseClippingIdentityCutoff >= 0.501)
		{
			params.Xdropcutoff = std::max((double)params.Xdropcutoff, (double)50 * (100 * (params.preciseClippingIdentityCutoff / (1.0 - params.preciseClippingIdentityCutoff) + 1.0)));
		}
	}
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
	if (outputAlns.size() == 0 && params.outputCorrectedFile == "" && params.outputCorrectedClippedFile == "")
	{
		std::cerr << "one of alignments-out, corrected-out or corrected-clipped-out must be given" << std::endl;
		paramError = true;
	}
	for (std::string file : outputAlns)
	{
		if (file.size() >= 4 && file.substr(file.size()-4) == ".gam")
		{
			params.outputGAMFile = file;
		}
		else if (file.size() >= 5 && file.substr(file.size()-5) == ".json")
		{
			params.outputJSONFile = file;
		}
		else if (file.size() >= 4 && file.substr(file.size()-4) == ".gaf")
		{
			params.outputGAFFile = file;
		}
		else
		{
			std::cerr << "unknown output alignment format (" << file << "), must be either .gaf, .gam or .json" << std::endl;
			paramError = true;
		}
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
	if (params.numThreads < 1)
	{
		std::cerr << "number of threads must be >= 1" << std::endl;
		paramError = true;
	}
	if (params.alignmentBandwidth < 1)
	{
		std::cerr << "alignment bandwidth must be >= 1" << std::endl;
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
	if (params.minimizerDiscardMostNumerousFraction < 0 || params.minimizerDiscardMostNumerousFraction >= 1)
	{
		std::cerr << "Minimizer discard fraction must be 0 <= x < 1" << std::endl;
		paramError = true;
	}
	if (params.minimizerSeedDensity < 0 && params.minimizerSeedDensity != -1)
	{
		std::cerr << "Minimizer density can't be negative" << std::endl;
		paramError = true;
	}
	if (params.multimapScoreFraction < 0)
	{
		std::cerr << "--multimap-score-fraction cannot be less than 0" << std::endl;
		paramError = true;
	}
	if (params.multimapScoreFraction > 1)
	{
		std::cerr << "--multimap-score-fraction cannot be more than 1" << std::endl;
		paramError = true;
	}
	if (params.preciseClippingIdentityCutoff < 0.501 || params.preciseClippingIdentityCutoff > 0.999)
	{
		std::cerr << "precise clipping identity cutoff must be between 0.501 and 0.999" << std::endl;
		paramError = true;
	}
	if (params.Xdropcutoff < 1)
	{
		std::cerr << "X-drop score cutoff must be > 1" << std::endl;
		paramError = true;
	}
	if (params.maxClusterExtend == 0)
	{
		std::cerr << "--max-cluster-extend cannot be 0" << std::endl;
		paramError = true;
	}
	int pickedSeedingMethods = ((params.dynamicRowStart) ? 1 : 0) + ((params.seedFiles.size() > 0) ? 1 : 0) + ((params.mumCount != 0) ? 1 : 0) + ((params.memCount != 0) ? 1 : 0) + ((params.minimizerSeedDensity != 0) ? 1 : 0);
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

	if (params.outputCorrectedFile.size() >= 3 && params.outputCorrectedFile.substr(params.outputCorrectedFile.size()-3) == ".gz")
	{
		params.compressCorrected = true;
	}

	if (params.outputCorrectedClippedFile.size() >= 3 && params.outputCorrectedClippedFile.substr(params.outputCorrectedClippedFile.size()-3) == ".gz")
	{
		params.compressClipped = true;
	}

	alignReads(params);

	return 0;
}
