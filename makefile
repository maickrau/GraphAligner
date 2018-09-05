CC=gcc
GPP=g++
CPPFLAGS=-Wall -std=c++14 -O3 -g

ODIR=obj
BINDIR=bin

LIBS=-lm -lz

DEPS = fastqloader.h BigraphToDigraph.h Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h GfaGraph.h ByteStuff.h

_OBJ = Aligner.o AlignerMain.o fastqloader.o BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o GfaGraph.o ByteStuff.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -static-libstdc++

$(shell mkdir -p bin)
$(shell mkdir -p obj)

all: $(BINDIR)/Aligner

$(ODIR)/%.o: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/Aligner: $(OBJ)
	$(GPP) -o $@ $^ $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ $(LINKFLAGS)

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
