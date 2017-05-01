IDIR=../include/
CC=gcc
GPP=g++
CFLAGS=-I$(IDIR) -Wall -msse4.1
CPPFLAGS=-I$(IDIR) -Wall

ODIR=obj
BINDIR=bin

LIBS=-lm -lz -lprotobuf

DEPS = gssw.h vg.pb.h

_OBJ = gssw.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))
_OBJPP = GsswWrapper.opp vg.pb.opp
OBJPP = $(patsubst %, $(ODIR)/%, $(_OBJPP))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.opp: %.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

$(BINDIR)/wrapper: $(OBJ) $(OBJPP)
	$(GPP) -o $@ $^ $(CPPFLAGS) $(LIBS)

all: $(BINDIR)/wrapper

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*