GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Iconcurrentqueue -IBBHash -Izstr/src -Iparallel-hashmap/parallel_hashmap/ `pkg-config --cflags protobuf` `pkg-config --cflags libsparsehash` `pkg-config --cflags mummer` -fopenmp -Wno-unused-parameter

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=-lm -lz -lboost_serialization -lboost_program_options `pkg-config --libs mummer`  `pkg-config --libs protobuf` -lsdsl
JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`

_DEPS = vg.pb.h fastqloader.h GraphAlignerWrapper.h vg.pb.h BigraphToDigraph.h stream.hpp Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h GfaGraph.h AlignmentCorrectnessEstimation.h MummerSeeder.h ReadCorrection.h MinimizerSeeder.h AlignmentSelection.h EValue.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = Aligner.o vg.pb.o fastqloader.o BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o GraphAlignerWrapper.o GfaGraph.o AlignmentCorrectnessEstimation.o MummerSeeder.o ReadCorrection.o MinimizerSeeder.o AlignmentSelection.o EValue.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ $(JEMALLOCFLAGS) `pkg-config --libs libdivsufsort` `pkg-config --libs libdivsufsort64`

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/GraphAligner: $(ODIR)/AlignerMain.o $(OBJ)
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/GraphAlignerWrapper.o: $(SRCDIR)/GraphAlignerWrapper.cpp $(SRCDIR)/GraphAligner.h $(SRCDIR)/NodeSlice.h $(SRCDIR)/WordSlice.h $(SRCDIR)/ArrayPriorityQueue.h $(SRCDIR)/ComponentPriorityQueue.h $(SRCDIR)/GraphAlignerVGAlignment.h $(SRCDIR)/GraphAlignerGAFAlignment.h $(SRCDIR)/GraphAlignerBitvectorBanded.h $(SRCDIR)/GraphAlignerBitvectorCommon.h $(SRCDIR)/GraphAlignerCommon.h $(SRCDIR)/DijkstraQueue.h $(SRCDIR)/GraphAlignerBitvectorDijkstra.h $(DEPS)

$(ODIR)/AlignerMain.o: $(SRCDIR)/AlignerMain.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(ODIR)/vg.pb.o: $(SRCDIR)/vg.pb.cc
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(SRCDIR)/%.pb.cc $(SRCDIR)/%.pb.h: $(SRCDIR)/%.proto
	protoc -I=$(SRCDIR) --cpp_out=$(SRCDIR) $<

$(BINDIR)/FusionFinder: $(SRCDIR)/FusionFinder.cpp $(OBJ)
	$(GPP) -o $@ $^ $(LINKFLAGS) -DVERSION="\"$(VERSION)\""

$(BINDIR)/ExtractPathSequence: $(SRCDIR)/ExtractPathSequence.cpp $(ODIR)/CommonUtils.o $(ODIR)/GfaGraph.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/SelectLongestAlignment: $(SRCDIR)/SelectLongestAlignment.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/Postprocess: $(SRCDIR)/Postprocess.cpp $(ODIR)/AlignmentSelection.o $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/AlignmentSubsequenceIdentity: $(SRCDIR)/AlignmentSubsequenceIdentity.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/UntipRelative: $(SRCDIR)/UntipRelative.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/PickAdjacentAlnPairs: $(SRCDIR)/PickAdjacentAlnPairs.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ExtractCorrectedReads: $(SRCDIR)/ExtractCorrectedReads.cpp $(ODIR)/ReadCorrection.o $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/GfaGraph.o $(ODIR)/fastqloader.o $(ODIR)/ThreadReadAssertion.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

all: $(BINDIR)/GraphAligner $(BINDIR)/ExtractPathSequence $(BINDIR)/SelectLongestAlignment $(BINDIR)/AlignmentSubsequenceIdentity $(BINDIR)/PickAdjacentAlnPairs $(BINDIR)/ExtractCorrectedReads $(BINDIR)/UntipRelative

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	rm -f $(SRCDIR)/vg.pb.cc
	rm -f $(SRCDIR)/vg.pb.h
