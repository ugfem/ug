##############################################################################
#																			 #
#	Makefile for all modules of ug version 3								 #
#																			 #
#	created 16 December 1994												 #
#																			 #
##############################################################################

# load architecture dependent makefile	
include ug.arch

# architecture-dependent makefile entries
include machines/mk.$(ARCHDIR)

# include switches
include ug.conf

# the following list may be extended
MODULES = LOW DEV DOM GM NUMERICS GRAPH GRAPE UI GG $(PMODULES)
UGMODULES = LOW GM NUMERICS GRAPH GRAPE UI GG $(PMODULES)

# dimension dependent targets
version = $(DIM)Dversion

# local C compiler flags
LCFLAGS = -Ilow -Idddif -Idev -Idom -Idom/$(DOM_MODULE) -Igm -Igraph -Iui -Inumerics -Igg

# object files for both dimensions
OBJECTS = initug.o

# make all
all: $(MODULES) $(OBJECTS)
	ar $(ARFLAGS) lib/libug$(LIBSUFFIX).a $(OBJECTS)
	echo "libug, libdom and libdev compiled"

uglib: $(UGMODULES) $(OBJECTS) 
	ar $(ARFLAGS) lib/libug$(LIBSUFFIX).a $(OBJECTS)
	echo "libug compiled"

LOW: include
	cd low; make -f Makefile.low; cd ..;

DDDIF: include
	cd dddif; make -f Makefile.dddif; cd ..;

DEV: include
	cd dev; make -f Makefile.dev; cd ..;

DOM: include
	cd dom; make -f Makefile.dom; cd ..;

GM: include
	cd gm; make -f Makefile.gm $(version); cd ..;

NUMERICS: include
	cd numerics; make -f Makefile.numerics $(version); cd ..;

GRAPH: include
	cd graph; make -f Makefile.graph $(version); cd ..;

GRAPE: include
	cd grape; make -f Makefile.grape $(version)GRAPE$(GRAPE); cd ..;

UI: include
	cd ui; make -f Makefile.ui $(version); cd ..;

GG: include
	cd gg; make -f Makefile.gg $(version); cd ..;
	cd gg3d; make -f Makefile.gg3d $(version); cd ..;

include:
	ugmakelinks;

# default rule
.c.o:
	$(CC) $(UGDEFINES) $(CFLAGS) $(LCFLAGS) $<

clean:
	rm -f $(OBJECTS)
	cd low; make -f Makefile.low clean; cd ..;
	cd dddif; make -f Makefile.dddif clean; cd ..;
	cd dev; make -f Makefile.dev clean; cd ..;
	cd dom; make -f Makefile.dom clean; cd ..;
	cd gm; make -f Makefile.gm clean; cd ..;
	cd numerics; make -f Makefile.numerics clean; cd ..;
	cd graph; make -f Makefile.graph clean; cd ..;
	cd ui; make -f Makefile.ui clean; cd ..;
	cd gg; make -f Makefile.gg clean; cd ..;
	cd gg3d; make -f Makefile.gg3d clean; cd ..;
#	cd machines/$(ARCHDIR); make clean; cd ..;

ifdef:
	cd gm; make -f Makefile.gm clean; cd ..;
	cd numerics; make -f Makefile.numerics clean; cd ..;
	cd graph; make -f Makefile.graph clean; cd ..;
	cd dom; make -f Makefile.dom clean; cd ..;
	cd ui; rm commands.o ; cd ..;
	rm -f initug.o;










