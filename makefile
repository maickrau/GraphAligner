CC=gcc
GPP=g++
CPPFLAGS=-Wall -std=c++14 -O3 -g `pkg-config --cflags mummer`

ODIR=obj
BINDIR=bin

LIBS=`pkg-config --libs mummer` -lm -lprotobuf -lz -lboost_serialization -lumdmummer -lboost_program_options -lsdsl -ldivsufsort -ldivsufsort64
JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`

DEPS = vg.pb.h fastqloader.h GraphAlignerWrapper.h vg.pb.h BigraphToDigraph.h stream.hpp Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h GfaGraph.h AlignmentCorrectnessEstimation.h ByteStuff.h MummerSeeder.h

_OBJ = Aligner.o AlignerMain.o vg.pb.o fastqloader.o BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o GraphAlignerWrapper.o GfaGraph.o AlignmentCorrectnessEstimation.o ByteStuff.o MummerSeeder.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ -fopenmp $(JEMALLOCFLAGS)

GITCOMMIT := $(shell git rev-parse HEAD)
GITBRANCH := $(shell git rev-parse --abbrev-ref HEAD)
GITDATE := $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(ODIR)/GraphAlignerWrapper.o: GraphAlignerWrapper.cpp GraphAligner.h NodeSlice.h WordSlice.h ArrayPriorityQueue.h GraphAlignerVGAlignment.h GraphAlignerBitvectorBanded.h GraphAlignerBitvectorCommon.h GraphAlignerCommon.h GraphAlignerSubgraphExtraction.h $(DEPS)

$(ODIR)/AlignerMain.o: AlignerMain.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DGITBRANCH=\"$(GITBRANCH)\" -DGITCOMMIT=\"$(GITCOMMIT)\" -DGITDATE="\"$(GITDATE)\""

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/Aligner: $(OBJ)
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/SimulateReads: SimulateReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ReverseReads: ReverseReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/SupportedSubgraph: SupportedSubgraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/MafToAlignment: MafToAlignment.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ExtractPathSequence: ExtractPathSequence.cpp $(ODIR)/CommonUtils.o $(ODIR)/GfaGraph.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ExtractPathSubgraphNeighbourhood: ExtractPathSubgraphNeighbourhood.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/VisualizeAlignment: VisualizeAlignment.cpp $(ODIR)/AlignmentCorrectnessEstimation.o $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/NodePosCsv: NodePosCsv.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ExtractExactPathSubgraph: ExtractExactPathSubgraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/EstimateRepeatCount: EstimateRepeatCount.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/PickMummerSeeds: PickMummerSeeds.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/SelectLongestAlignment: SelectLongestAlignment.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/Postprocess: Postprocess.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/AlignmentSubsequenceIdentity: AlignmentSubsequenceIdentity.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/BruteForceExactPrefixSeeds: BruteForceExactPrefixSeeds.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/UntipRelative: UntipRelative.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/PickAdjacentAlnPairs: PickAdjacentAlnPairs.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/TestMummerSeeder: TestMummerSeeder.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/MummerSeeder.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

all: $(BINDIR)/Aligner $(BINDIR)/SimulateReads $(BINDIR)/ReverseReads $(BINDIR)/SupportedSubgraph $(BINDIR)/MafToAlignment $(BINDIR)/ExtractPathSequence $(BINDIR)/ExtractPathSubgraphNeighbourhood $(BINDIR)/VisualizeAlignment $(BINDIR)/NodePosCsv $(BINDIR)/ExtractExactPathSubgraph $(BINDIR)/EstimateRepeatCount $(BINDIR)/PickMummerSeeds $(BINDIR)/SelectLongestAlignment $(BINDIR)/Postprocess $(BINDIR)/AlignmentSubsequenceIdentity $(BINDIR)/BruteForceExactPrefixSeeds $(BINDIR)/PickAdjacentAlnPairs

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
