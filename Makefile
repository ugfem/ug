##############################################################################
#																			 #
#	Makefile for all modules of ug version 3								 #
#																			 #
#	created 16 December 1994												 #
#																			 #
##############################################################################

# include configuration and all makefile macro definitions
include ug.conf

# the following list may be extended
MODULES = LOW DEV DOM GM NUMERICS GRAPHICS UI $(MODEL_TARGET)
UGMODULES = LOW GM NUMERICS GRAPHICS UI $(MODEL_TARGET)

# modules for ug server daemon
UGDMODULES = LOW DEV 


# local C compiler flags
LCFLAGS = -I./include

# object files for both dimensions
OBJECTS = initug.o


##############################################################################

# make all
all: include $(OBJECTS) $(MODULES)
	$(ARCH_AR) $(ARCH_ARFLAGS) lib/libug$(UG_LIBSUFFIX).a $(OBJECTS)
	echo "libug, libdom and libdev compiled"

uglib: include $(OBJECTS) $(UGMODULES)
	$(ARCH_AR) $(ARCH_ARFLAGS) lib/libug$(UG_LIBSUFFIX).a $(OBJECTS)
	echo "libug compiled"

UGD:  include $(UGDMODULES) ugd.o
	$(UG_LINK) -o bin/ugd $(ARCH_LFLAGS) ugd.o lib/libdev.a lib/libug$(UG_LIBSUFFIX).a $(UG_LFLAGS)
	echo "ugd compiled"


##############################################################################

LOW: include
	cd low; make -f Makefile.low; cd ..;

DEV: include
	cd dev; make -f Makefile.dev; cd ..;

DOM: include
	cd dom; make -f Makefile.dom; cd ..;

GM: include
	cd gm; make -f Makefile.gm; cd ..;

NUMERICS: include
	cd numerics; make -f Makefile.numerics; cd ..;

GRAPHICS: include
	cd graphics; make -f Makefile.graphics; cd ..;

UI: include
	cd ui; make -f Makefile.ui; cd ..;


##############################################################################

# targets for MODEL, exactly one will be chosen

SEQUENTIAL:

PARALLEL: include
	cd parallel; make -f Makefile.parallel; cd ..;


SEQUENTIAL_clean:

PARALLEL_clean:
	cd parallel; make -f Makefile.parallel clean; cd ..;


##############################################################################

include:
	ugmakelinks;

# default rule
.c.o: 
	$(ARCH_CC) $(UG_CFLAGS) $(LCFLAGS) $<


##############################################################################

clean: $(MODEL_TARGET)_clean
	rm -f $(OBJECTS)
	cd low; make -f Makefile.low clean; cd ..;
	cd dev; make -f Makefile.dev clean; cd ..;
	cd dom; make -f Makefile.dom clean; cd ..;
	cd gm; make -f Makefile.gm clean; cd ..;
	cd numerics; make -f Makefile.numerics clean; cd ..;
	cd graphics; make -f Makefile.graphics clean; cd ..;
	cd ui; make -f Makefile.ui clean; cd ..;

ifdef: $(MODEL_TARGET)_clean
	cd gm; make -f Makefile.gm clean; cd ..;
	cd numerics; make -f Makefile.numerics clean; cd ..;
	cd graphics; make -f Makefile.graphics clean; cd ..;
	cd dom; make -f Makefile.dom clean; cd ..;
	cd ui; rm commands.o ; cd ..;
	rm -f initug.o;
