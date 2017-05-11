CC=gcc
GPP=g++
CFLAGS=-Wall -msse4.1 -O2 -g -pg -DNDEBUG
CPPFLAGS=-Wall -std=c++11 -msse4.1 -O3 -g -pg -DNDEBUG

ODIR=obj
BINDIR=bin

LIBS=-lm -lprotobuf -lz 

DEPS = gssw.h vg.pb.h

_OBJ = gssw.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))
_OBJPP = GsswWrapper.opp vg.pb.opp fastqloader.opp TopologicalSort.opp SubgraphFromSeed.opp
OBJPP = $(patsubst %, $(ODIR)/%, $(_OBJPP))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.opp: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/wrapper: $(OBJ) $(OBJPP)
	$(GPP) -o $@ $^ $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(BINDIR)/ReadIndexToId: $(OBJ) $(OBJPP)
	$(GPP) -o $@ ReadIndexToId.cpp $(ODIR)/vg.pb.opp $(ODIR)/fastqloader.opp $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

all: $(BINDIR)/wrapper $(BINDIR)/ReadIndexToId

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*