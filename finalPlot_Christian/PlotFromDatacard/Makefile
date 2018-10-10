TARGET = Main

C=g++
LD=g++

#O2 for optimization, g for debugging, pg for profiling
SPECIALFLAGS= -g #-g -Wall -pg#-O2
ROOTAUXCFLAGS=$(shell root-config --auxcflags)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) #-lMath # -lMinuit
CFLAGS = $(SPECIALFLAGS) -Wall $(ROOTAUXCFLAGS)
LFLAGS = $(SPECIALFLAGS) #-lz 

RCXX=$(SPECIALFLAGS) -Wno-deprecated -Wall $(ROOTCFLAGS) -I$(HOME)/include
RLXX=$(LFLAGS) $(ROOTLIBS)  -L$(HOME)/lib 

SRCEXT = cc
SRCS = $(shell find ./ -type f -name '*.$(SRCEXT)')

%.o: %.cc
		$(C) $(RCXX) -c $<

all: ${TARGET}

clean:
	@rm -f *.o 
	@rm -f ${TARGET}

$(TARGET): $(SRCS:.cc=.o) 
		$(LD)  $(SRCS:.cc=.o) $(RLXX) -o $@
		@echo '-> executable $@ created.'
