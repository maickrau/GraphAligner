PLATFORM=$(shell uname -s)
GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Iconcurrentqueue -IBBHash -Izstr/src -Iparallel-hashmap/parallel_hashmap/ `pkg-config --cflags protobuf` `pkg-config --cflags libsparsehash` -Wno-unused-parameter -IMEMfinder/src -I`jemalloc-config --includedir`

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=-lm -lz -lboost_program_options `pkg-config --libs protobuf` -lsdsl
JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`

_DEPS = vg.pb.h fastqloader.h GraphAlignerWrapper.h vg.pb.h BigraphToDigraph.h stream.hpp Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h GfaGraph.h ReadCorrection.h MinimizerSeeder.h AlignmentSelection.h EValue.h MEMSeeder.h DNAString.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = Aligner.o vg.pb.o fastqloader.o BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o GraphAlignerWrapper.o GfaGraph.o ReadCorrection.o MinimizerSeeder.o AlignmentSelection.o EValue.o MEMSeeder.o DNAString.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

ifeq ($(PLATFORM),Linux)
   JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`
   LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ $(JEMALLOCFLAGS) `pkg-config --libs libdivsufsort` `pkg-config --libs libdivsufsort64`
else
   CPPFLAGS += -D_LIBCPP_DISABLE_AVAILABILITY
   JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -ljemalloc `jemalloc-config --libs`
   LINKFLAGS = $(CPPFLAGS) $(LIBS) -lpthread -pthread -static-libstdc++ $(JEMALLOCFLAGS) `pkg-config --libs libdivsufsort` `pkg-config --libs libdivsufsort64`
endif

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/GraphAligner: $(ODIR)/AlignerMain.o $(OBJ) MEMfinder/lib/memfinder.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/GraphAlignerWrapper.o: $(SRCDIR)/GraphAlignerWrapper.cpp $(SRCDIR)/GraphAligner.h $(SRCDIR)/NodeSlice.h $(SRCDIR)/WordSlice.h $(SRCDIR)/ArrayPriorityQueue.h $(SRCDIR)/ComponentPriorityQueue.h $(SRCDIR)/GraphAlignerVGAlignment.h $(SRCDIR)/GraphAlignerGAFAlignment.h $(SRCDIR)/GraphAlignerBitvectorBanded.h $(SRCDIR)/GraphAlignerBitvectorCommon.h $(SRCDIR)/GraphAlignerCommon.h $(DEPS)

$(ODIR)/AlignerMain.o: $(SRCDIR)/AlignerMain.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(ODIR)/vg.pb.o: $(SRCDIR)/vg.pb.cc
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(SRCDIR)/%.pb.cc $(SRCDIR)/%.pb.h: $(SRCDIR)/%.proto
	protoc -I=$(SRCDIR) --cpp_out=$(SRCDIR) $<

MEMfinder/lib/memfinder.a:
	$(MAKE) -C MEMfinder lib

all: $(BINDIR)/GraphAligner

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	rm -f $(SRCDIR)/vg.pb.cc
	rm -f $(SRCDIR)/vg.pb.h
	$(MAKE) -C MEMfinder clean
