RED='\e[1;31m'
GREEN='\e[1;32m'
YELLOW='\e[1;33m'
BLUE='\e[1;34m'
WHITE='\e[0m'
VIOLET='\e[1;35m'
CYAN='\e[1;36m'
red='\e[0;31m'
green='\e[0;32m'
yellow='\e[0;33m'
blue='\e[0;34m'
violet='\e[0;35m'
cyan='\e[0;36m'

Reset=`tput sgr0`

#CC = g++ `root-config --cflags`
CC = g++ 
CFLAGS = -fPIC $(shell root-config --cflags)

ROOTLIBS = $(shell root-config --libs --glibs)
					

LSDIR = $(shell ls)

OBJS1 = Decaymain.C ODEsolver.C element.C ODEdecay.C MCdecay.C  Decayfunc.C
OBJS2 = MyMainFrame.C graph.C ODEsolver2.C ./atalho/cFCgraphics.C 


parte1: rDecay.exe
parte2:	rPlanets.exe  

rDecay.exe: $(OBJS1)
	@echo list of dependencies $^
	@echo -e creating executable... $(GREEN) $@ $(WHITE)
	@g++ $(CFLAGS) -o $@ $^ $(ROOTLIBS)

rPlanets.exe: $(OBJS2)
	rootcint -f graph.C -c MyMainFrame.h 
	@echo list of dependencies $^
	@echo -e creating executable... $(GREEN) $@ $(WHITE)
	@g++ $(CFLAGS) -o $@ -I ./atalho $^ $(ROOTLIBS)

%.o: %.C
	@echo -e compiling $(GREEN) $<  $(WHITE) to $(YELLOW) $@ $(WHITE)
	@g++ $(CFLAGS) -c $< -O3

test: 
	@echo Hello! $@
	@echo CC = $(CC)
	@echo ROOTLIBS = $(ROOTLIB)
	@echo LSDIR = $(LSDIR)
	@echo CFLAGS = $(CFLAGS)

clean: 
	@echo cleaning...
	rm -f *.o
	rm -f *.exe
