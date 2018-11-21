CC=gcc
GPP=g++
CPPFLAGS=-Wall -Wextra -std=c++14 -O3 -g `pkg-config --cflags mummer` -fopenmp -Wno-unused-parameter

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=`pkg-config --libs mummer` -lm -lprotobuf -lz -lboost_serialization -lumdmummer -lboost_program_options -lsdsl -ldivsufsort -ldivsufsort64
JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`

_DEPS = vg.pb.h fastqloader.h GraphAlignerWrapper.h vg.pb.h BigraphToDigraph.h stream.hpp Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h GfaGraph.h AlignmentCorrectnessEstimation.h MummerSeeder.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = Aligner.o AlignerMain.o vg.pb.o fastqloader.o BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o GraphAlignerWrapper.o GfaGraph.o AlignmentCorrectnessEstimation.o MummerSeeder.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ $(JEMALLOCFLAGS)

GITCOMMIT := $(shell git rev-parse HEAD)
GITBRANCH := $(shell git rev-parse --abbrev-ref HEAD)
GITDATE := $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/Aligner: $(OBJ)
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/GraphAlignerWrapper.o: $(SRCDIR)/GraphAlignerWrapper.cpp $(SRCDIR)/GraphAligner.h $(SRCDIR)/NodeSlice.h $(SRCDIR)/WordSlice.h $(SRCDIR)/ArrayPriorityQueue.h $(SRCDIR)/GraphAlignerVGAlignment.h $(SRCDIR)/GraphAlignerBitvectorBanded.h $(SRCDIR)/GraphAlignerBitvectorCommon.h $(SRCDIR)/GraphAlignerCommon.h $(DEPS)

$(ODIR)/AlignerMain.o: $(SRCDIR)/AlignerMain.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DGITBRANCH=\"$(GITBRANCH)\" -DGITCOMMIT=\"$(GITCOMMIT)\" -DGITDATE="\"$(GITDATE)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/SimulateReads: $(SRCDIR)/SimulateReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ReverseReads: $(SRCDIR)/ReverseReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/SupportedSubgraph: $(SRCDIR)/SupportedSubgraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/MafToAlignment: $(SRCDIR)/MafToAlignment.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ExtractPathSequence: $(SRCDIR)/ExtractPathSequence.cpp $(ODIR)/CommonUtils.o $(ODIR)/GfaGraph.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ExtractPathSubgraphNeighbourhood: $(SRCDIR)/ExtractPathSubgraphNeighbourhood.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/VisualizeAlignment: $(SRCDIR)/VisualizeAlignment.cpp $(ODIR)/AlignmentCorrectnessEstimation.o $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/NodePosCsv: $(SRCDIR)/NodePosCsv.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ExtractExactPathSubgraph: $(SRCDIR)/ExtractExactPathSubgraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/EstimateRepeatCount: $(SRCDIR)/EstimateRepeatCount.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/PickMummerSeeds: $(SRCDIR)/PickMummerSeeds.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/SelectLongestAlignment: $(SRCDIR)/SelectLongestAlignment.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/Postprocess: $(SRCDIR)/Postprocess.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/AlignmentSubsequenceIdentity: $(SRCDIR)/AlignmentSubsequenceIdentity.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/BruteForceExactPrefixSeeds: $(SRCDIR)/BruteForceExactPrefixSeeds.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/UntipRelative: $(SRCDIR)/UntipRelative.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/PickAdjacentAlnPairs: $(SRCDIR)/PickAdjacentAlnPairs.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

all: $(BINDIR)/Aligner $(BINDIR)/SimulateReads $(BINDIR)/ReverseReads $(BINDIR)/SupportedSubgraph $(BINDIR)/MafToAlignment $(BINDIR)/ExtractPathSequence $(BINDIR)/ExtractPathSubgraphNeighbourhood $(BINDIR)/VisualizeAlignment $(BINDIR)/NodePosCsv $(BINDIR)/ExtractExactPathSubgraph $(BINDIR)/EstimateRepeatCount $(BINDIR)/PickMummerSeeds $(BINDIR)/SelectLongestAlignment $(BINDIR)/Postprocess $(BINDIR)/AlignmentSubsequenceIdentity $(BINDIR)/BruteForceExactPrefixSeeds $(BINDIR)/PickAdjacentAlnPairs

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
