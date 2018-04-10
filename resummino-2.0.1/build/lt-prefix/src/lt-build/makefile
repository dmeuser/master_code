# --- variables defined by configure ---

SRC = /home/home4/institut_1b/dmeuser/master_code/resummino-2.0.1/lib/LoopTools-2.13/src
PREFIX = /home/home4/institut_1b/dmeuser/master_code/resummino-2.0.1/build/lt-prefix
LIBDIRSUFFIX = 

EXE = 
DEF = -D
NOUNDERSCORE = 0

FC = gfortran
FFLAGS = -O3 -march=native -g -ffixed-line-length-none -fno-range-check $(DEF)QUAD=1 $(DEF)U77EXT=0

CC = clang
CFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64 -DNOUNDERSCORE=$(NOUNDERSCORE) -DBIGENDIAN=0

CXX = g++
CXXFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64 

ML = 1
MCC = mcc
MCFLAGS = -st -b64
MCLIBS = -lpthread

LDFLAGS = -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin/../libexec/gcc/x86_64-pc-linux-gnu/5.3.0/liblto_plugin.so -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin/../lib/gcc/x86_64-pc-linux-gnu/5.3.0 -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin/../lib/gcc -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin/../lib/gcc/x86_64-pc-linux-gnu/5.3.0/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin/../lib/gcc/x86_64-pc-linux-gnu/5.3.0/../../../../x86_64-pc-linux-gnu/lib -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin/../lib/gcc/x86_64-pc-linux-gnu/5.3.0/../../.. -lgfortran -lm -lgcc_s -lgcc -lquadmath -lm -lgcc_s -lgcc -lgcc_s -lgcc -m64

AR = ar
RANLIB = ranlib
NM = nm
DLLTOOL = 

PARALLEL = -j 8

# --- end defs by configure ---


BLD = build$(QUADSUFFIX)

LIBDIR = $(PREFIX)/lib$(LIBDIRSUFFIX)
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

LIB = libooptools$(QUADSUFFIX).a
FE = lt$(QUADSUFFIX)$(EXE)
MFE = LoopTools$(QUADSUFFIX)$(EXE)
INCLUDE = $(BLD)/looptools.h $(BLD)/clooptools.h

ARGS = $(PARALLEL) \
  LIB="$(LIB)" \
  FE="$(FE)" \
  MFE="$(MFE)" \
  EXE="$(EXE)" \
  DEF="$(DEF)" \
  NOUNDERSCORE="$(NOUNDERSCORE)" \
  XFC="$(FC$(QUADSUFFIX)) $(FFLAGS) $(FFLAGS-quad) -I." \
  F90="$(F90)" \
  CC="$(CC)" \
  CFLAGS="$(CFLAGS) $(CFLAGS-quad)" \
  CXX="$(CXX)" \
  CXXFLAGS="$(CXXFLAGS)" \
  ML="$(ML)" \
  MCC="$(MCC)" \
  MCFLAGS="$(MCFLAGS)" \
  MCLIBS="$(MCLIBS)" \
  AR="$(AR)" \
  RANLIB="$(RANLIB)" \
  NM="$(NM)" \
  DLLTOOL="$(DLLTOOL)" \
  LDFLAGS="$(LDFLAGS)" \
  LIBPATH="$(LIBPATH)" \
  OBJS-quad="$(OBJS-quad)"


default all lib frontend mma: force
	cd $(BLD) && $(MAKE) $(ARGS) $@

install: lib frontend
	-mkdir $(PREFIX)
	-mkdir $(LIBDIR) $(BINDIR) $(INCLUDEDIR)
	cp -p $(BLD)/$(LIB) $(LIBDIR)
	cp -p $(INCLUDE$(QUADSUFFIX)) $(INCLUDEDIR)
	strip $(BLD)/$(FE)
	cp -p $(BLD)/fcc $(BLD)/$(FE) $(BINDIR)
	rm -f $(BINDIR)/f++
	ln -s fcc $(BINDIR)/f++
	test ! -f $(BLD)/$(MFE) || { strip $(BLD)/$(MFE) ; cp -p $(BLD)/$(MFE) $(BINDIR); }

force: $(BLD)/timestamp

$(BLD)/timestamp:
	-mkdir $(BLD)
	find "`cd $(SRC) ; pwd`" -follow -exec ln -sf {} $(BLD) \;
	touch $(BLD)/timestamp

clean:
	rm -fr $(BLD)

