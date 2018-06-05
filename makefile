CC=gcc
GPP=g++
CPPFLAGS=-Wall -std=c++14 -O3 -g

ODIR=obj
BINDIR=bin

LIBS=-lm -lprotobuf -lz
LINKFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`

DEPS = vg.pb.h fastqloader.h GraphAlignerWrapper.h vg.pb.h BigraphToDigraph.h stream.hpp Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h GfaGraph.h AlignmentCorrectnessEstimation.h ByteStuff.h

_OBJ = Aligner.o AlignerMain.o vg.pb.o fastqloader.o BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o GraphAlignerWrapper.o GfaGraph.o AlignmentCorrectnessEstimation.o ByteStuff.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(ODIR)/GraphAlignerWrapper.o: GraphAlignerWrapper.cpp GraphAligner.h NodeSlice.h WordSlice.h ArrayPriorityQueue.h GraphAlignerVGAlignment.h GraphAlignerBitvectorBanded.h GraphAlignerBitvectorCommon.h GraphAlignerCommon.h $(DEPS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/Aligner: $(OBJ)
	$(GPP) -o $@ $^ $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ $(LINKFLAGS)

$(BINDIR)/ReadIndexToId: $(OBJ)
	$(GPP) -o $@ ReadIndexToId.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/CompareAlignments: $(OBJ)
	$(GPP) -o $@ CompareAlignments.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/MergeGraphs: $(OBJ)
	$(GPP) -o $@ MergeGraphs.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/SimulateReads: $(OBJ)
	$(GPP) -o $@ SimulateReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/ReverseReads: $(OBJ)
	$(GPP) -o $@ ReverseReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/PickSeedHits: $(OBJ)
	$(GPP) -o $@ PickSeedHits.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ $(LINKFLAGS)

$(BINDIR)/AlignmentSequenceInserter: $(OBJ)
	$(GPP) -o $@ AlignmentSequenceInserter.cpp $(ODIR)/CommonUtils.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/SupportedSubgraph: $(OBJ)
	$(GPP) -o $@ SupportedSubgraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/MafToAlignment: $(OBJ)
	$(GPP) -o $@ MafToAlignment.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/ExtractPathSequence: $(OBJ)
	$(GPP) -o $@ ExtractPathSequence.cpp $(ODIR)/CommonUtils.o $(ODIR)/GfaGraph.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/AlignmentOverlap: $(OBJ)
	$(GPP) -o $@ AlignmentOverlap.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/Bluntify: $(OBJ)
	$(GPP) -o $@ Bluntify.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/ExtractPathSubgraphNeighbourhood: $(OBJ)
	$(GPP) -o $@ ExtractPathSubgraphNeighbourhood.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/MergeGfas: $(OBJ)
	$(GPP) -o $@ MergeGfas.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/VisualizeAlignment: $(OBJ)
	$(GPP) -o $@ VisualizeAlignment.cpp $(ODIR)/AlignmentCorrectnessEstimation.o $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/SelectPartials: $(OBJ)
	$(GPP) -o $@ SelectPartials.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/NodePosCsv: $(OBJ)
	$(GPP) -o $@ NodePosCsv.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/ExtractExactPathSubgraph: $(OBJ)
	$(GPP) -o $@ ExtractExactPathSubgraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/EstimateRepeatCount: $(OBJ)
	$(GPP) -o $@ EstimateRepeatCount.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/GfaGraph.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/SelectFullAlignments: $(OBJ)
	$(GPP) -o $@ SelectFullAlignments.cpp $(ODIR)/CommonUtils.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/AddNodeNames: $(OBJ)
	$(GPP) -o $@ AddNodeNames.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

all: $(BINDIR)/Aligner $(BINDIR)/ReadIndexToId $(BINDIR)/CompareAlignments $(BINDIR)/SimulateReads $(BINDIR)/ReverseReads $(BINDIR)/PickSeedHits $(BINDIR)/AlignmentSequenceInserter $(BINDIR)/MergeGraphs $(BINDIR)/SupportedSubgraph $(BINDIR)/MafToAlignment $(BINDIR)/ExtractPathSequence $(BINDIR)/AlignmentOverlap $(BINDIR)/Bluntify $(BINDIR)/ExtractPathSubgraphNeighbourhood $(BINDIR)/MergeGfas $(BINDIR)/VisualizeAlignment $(BINDIR)/SelectPartials $(BINDIR)/NodePosCsv $(BINDIR)/ExtractExactPathSubgraph $(BINDIR)/EstimateRepeatCount $(BINDIR)/SelectFullAlignments $(BINDIR)/AddNodeNames

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
