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
MODULES = LOW DEV GM NUMERICS GRAPH UI 

# dimension dependent targets
version = $(DIM)Dversion

# local C compiler flags
LCFLAGS = -Ilow -Idev -Igm -Igraph -Iui -Inumerics

# object files for both dimensions
OBJECTS = initug.o

# make all
all: $(MODULES) $(OBJECTS)
	ar $(ARFLAGS) lib/libug$(LIBSUFFIX).a $(OBJECTS)
	echo "ug 3 compiled"

LOW:
	cd low; make -f Makefile.low; cd ..;

DEV:
	cd dev; make -f Makefile.dev; cd ..;

GM:
	cd gm; make -f Makefile.gm $(version); cd ..;
	
NUMERICS:
	cd numerics; make -f Makefile.numerics $(version); cd ..;

GRAPH:
	cd graph; make -f Makefile.graph $(version); cd ..;
	
UI:
	cd ui; make -f Makefile.ui $(version); cd ..;
	
# default rule
.c.o:
	$(CC) $(UGDEFINES) $(CFLAGS) $(LCFLAGS) $<

clean:
	rm $(OBJECTS)
	cd low; make -f Makefile.low clean; cd ..;
	cd dev; make -f Makefile.dev clean; cd ..;
	cd gm; make -f Makefile.gm clean; cd ..;
	cd numerics; make -f Makefile.numerics clean; cd ..;
	cd graph; make -f Makefile.graph clean; cd ..;
	cd ui; make -f Makefile.ui clean; cd ..;

ifdef:
	cd gm; make -f Makefile.gm clean; cd ..;
	cd numerics; make -f Makefile.numerics clean; cd ..;
	cd graph; make -f Makefile.graph clean; cd ..;

