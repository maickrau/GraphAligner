CC=gcc
GPP=g++
CPPFLAGS=-Wall -std=c++14 -O3 -g

ODIR=obj
BINDIR=bin

LIBS=-lm -lprotobuf -lz -lboost_serialization

DEPS = vg.pb.h fastqloader.h GraphAligner.h SubgraphFromSeed.h TopologicalSort.h vg.pb.h BigraphToDigraph.h mfvs_graph.h mfvs_utils.h stream.hpp 2dArray.h SparseBoolMatrix.h SparseMatrix.h ssw_cpp.h Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h CycleCutCalculation.h

_OBJ = Aligner.o GsswWrapper.o vg.pb.o fastqloader.o TopologicalSort.o SubgraphFromSeed.o mfvs_graph.o mfvs_utils.o BigraphToDigraph.o ssw_cpp.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o CycleCutCalculation.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(ODIR)/ssw.o:
	$(CC) -c -o $@ ssw.c $(CPPFLAGS)

$(BINDIR)/wrapper: $(OBJ) $(ODIR)/ssw.o
	$(GPP) -o $@ $^ $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(BINDIR)/ReadIndexToId: $(OBJ)
	$(GPP) -o $@ ReadIndexToId.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/CompareAlignments: $(OBJ)
	$(GPP) -o $@ CompareAlignments.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/MergeGraphs: $(OBJ)
	$(GPP) -o $@ MergeGraphs.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/SimulateReads: $(OBJ)
	$(GPP) -o $@ SimulateReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/ReverseReads: $(OBJ)
	$(GPP) -o $@ ReverseReads.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/PickSeedHits: $(OBJ)
	$(GPP) -o $@ PickSeedHits.cpp $(ODIR)/CommonUtils.o $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/AlignmentSequenceInserter: $(OBJ)
	$(GPP) -o $@ AlignmentSequenceInserter.cpp $(ODIR)/CommonUtils.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/OrderProtobufs: $(OBJ)
	$(GPP) -o $@ OrderProtobufs.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(ODIR)/TopologicalSort.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/SupportedSubgraph: $(OBJ)
	$(GPP) -o $@ SupportedSubgraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(ODIR)/TopologicalSort.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/MafToAlignment: $(OBJ)
	$(GPP) -o $@ MafToAlignment.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(ODIR)/TopologicalSort.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/ExtractPathSequence: $(OBJ)
	$(GPP) -o $@ ExtractPathSequence.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/PreprocessGraph: $(OBJ)
	$(GPP) -o $@ PreprocessGraph.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/BigraphToDigraph.o $(ODIR)/vg.pb.o $(ODIR)/ssw_cpp.o $(ODIR)/ssw.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(BINDIR)/AlignmentOverlap: $(OBJ)
	$(GPP) -o $@ AlignmentOverlap.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/Bluntify: $(OBJ)
	$(GPP) -o $@ Bluntify.cpp $(ODIR)/CommonUtils.o $(ODIR)/ThreadReadAssertion.o $(ODIR)/fastqloader.o $(ODIR)/vg.pb.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

all: $(BINDIR)/wrapper $(BINDIR)/ReadIndexToId $(BINDIR)/CompareAlignments $(BINDIR)/SimulateReads $(BINDIR)/ReverseReads $(BINDIR)/PickSeedHits $(BINDIR)/AlignmentSequenceInserter $(BINDIR)/MergeGraphs $(BINDIR)/OrderProtobufs $(BINDIR)/SupportedSubgraph $(BINDIR)/MafToAlignment $(BINDIR)/ExtractPathSequence $(BINDIR)/PreprocessGraph $(BINDIR)/AlignmentOverlap $(BINDIR)/Bluntify

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
