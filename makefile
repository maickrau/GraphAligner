CC=gcc
GPP=g++
CPPFLAGS=-Wall -std=c++14 -O3 -g -pg

ODIR=obj
BINDIR=bin

LIBS=-lm -lprotobuf -lz 

DEPS = vg.pb.h fastqloader.h GraphAligner.h SubgraphFromSeed.h TopologicalSort.h vg.pb.h BigraphToDigraph.h mfvs_graph.h mfvs_utils.h stream.hpp 2dArray.h SparseBoolMatrix.h SparseMatrix.h

_OBJ = GsswWrapper.o vg.pb.o fastqloader.o TopologicalSort.o SubgraphFromSeed.o mfvs_graph.o mfvs_utils.o BigraphToDigraph.o SparseBoolMatrix.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/wrapper: $(OBJ)
	$(GPP) -o $@ $^ $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(BINDIR)/ReadIndexToId: $(OBJ)
	$(GPP) -o $@ ReadIndexToId.cpp $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

$(BINDIR)/CompareAlignments: $(OBJ)
	$(GPP) -o $@ CompareAlignments.cpp $(ODIR)/vg.pb.o $(ODIR)/fastqloader.o $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

all: $(BINDIR)/wrapper $(BINDIR)/ReadIndexToId $(BINDIR)/CompareAlignments

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*