# --- variables defined by configure ---

SRC = ./src
PREFIX = x86_64-Linux
LIBDIRSUFFIX = 64

EXE = 
DEF = -D
NOUNDERSCORE = 0

FC = gfortran
FFLAGS = -O3 -march=native -g -ffixed-line-length-none -fno-range-check $(DEF)QUAD=0 $(DEF)U77EXT=0

CC = gcc
CFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64 -DNOUNDERSCORE=$(NOUNDERSCORE) -DBIGENDIAN=0

CXX = g++
CXXFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -Wall -m64 -stdlib=libstdc++

ML = 0
MCC = mcc
MCFLAGS = -st -b64
MCLIBS = -lpthread

LDFLAGS = -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s -lgcc -lgcc_s -lgcc -m64

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

