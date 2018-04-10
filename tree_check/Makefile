LIBS = $(shell root-config --evelibs)
INCS = -I$(shell root-config --incdir) -I.

WARN = -Wall -Wshadow

OBJ = compare_ntuple.o
EXE = compare_ntuple

all: $(EXE) Makefile
	@echo "done"

%.o:%.cpp
	g++ -o $@ $+ -c -O2 $(INCS) $(WARN) -std=c++11

$(EXE): $(OBJ)
	g++ -O2 -o $@ $+ $$CMSSW_BASE/lib/slc6_amd64_gcc530/pluginTreeWriterTreeWriterAuto.so $(LIBS) $(WARN)

clean:
	@rm -f $(OBJ)


