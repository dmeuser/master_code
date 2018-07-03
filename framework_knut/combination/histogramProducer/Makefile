LIBS = $(shell root-config --evelibs)
INCS = -I$(shell root-config --incdir) -I.
OWNLIBS = $(shell echo $CMSSW_BASE/lib/$SCRAM_ARCH)

WARN = -Wall -Wshadow

PARTICLELIB = /home/home4/institut_1b/dmeuser/CMSSW_TreeWriter/CMSSW_8_0_26_patch1/lib/slc7_amd64_gcc530/pluginTreeWriterTreeWriterAuto.so
EXE = CombinationHistogramProducer
OBJ = CombinationHistogramProducer.o

all: $(EXE) Makefile

%.o: %.cc
	g++ -o $@ $+ -c -O2 $(INCS) $(WARN) -std=c++11

$(EXE): $(OBJ)
	g++ -O2 -o $@ $(PARTICLELIB) $+ $(LIBS) $(WARN)

clean:
	@rm -f *.cxx lib*h # remove rootcints generated source and header files
	@rm -f AutoDict_* *_h.d *_cc.d *_C.d *_cc_ACLiC_dict_rdict.pcm *_cc.so# cint generated files
	@rm -f $(OBJ)


