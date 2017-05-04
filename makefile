IDIR=/mnt/c/koulujutut/vg/vg/src
CC=gcc
GPP=g++
CFLAGS=-I$(IDIR) -Wall -msse4.1 -g -pg -O3 -DNDEBUG
CPPFLAGS=-I$(IDIR) -Wall -std=c++11 -msse4.1 -g -pg -O3 -DNDEBUG

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
	$(GPP) -o $@ $^ $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed

all: $(BINDIR)/wrapper

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*