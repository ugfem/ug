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
UGMODULES = LOW GM NP GRAPHICS UI $(MODEL_TARGET)
MODULES = DEV DOM $(UGMODULES)

# modules for ug server daemon
UGDMODULES = LOW DEV 

# object files for both dimensions
OBJECTS = initug.o

##############################################################################

# make all
all: include $(MODULES)
	make $(UG_LIB) 	
	echo "libug, libdom and libdev compiled"

uglib: include $(UGMODULES)
	make $(UG_LIB) 	
	echo "libug compiled"

UGD: include $(UGDMODULES) ugd.o
	make $(UG_LIB) 	
	$(UG_LINK) -o bin/ugd $(ARCH_LFLAGS) ugd.o lib/libdev.a \
                          $(UG_LIB) $(UG_LFLAGS)
	echo "ugd compiled"

$(UG_LIB): $(OBJECTS) 
	$(ARCH_AR) $(ARCH_ARFLAGS) $(UG_LIB) $(OBJECTS) 

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

NP: include
	cd np; make -f Makefile.np; cd ..;

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
	cd np; make -f Makefile.np clean; cd ..;
	cd graphics; make -f Makefile.graphics clean; cd ..;
	cd ui; make -f Makefile.ui clean; cd ..;

ifdef: $(MODEL_TARGET)_clean
	cd gm; make -f Makefile.gm clean; cd ..;
	cd numerics; make -f Makefile.numerics clean; cd ..;
	cd np; make -f Makefile.np clean; cd ..;
	cd graphics; make -f Makefile.graphics clean; cd ..;
	cd dom; make -f Makefile.dom clean; cd ..;
	cd ui; make -f Makefile.ui clean; cd ..;
	rm -f initug.o;

extract:
	$(ARCH_AR) $(ARCH_EXFLAGS) lib/libug$(UG_LIBSUFFIX).a $(OBJECTS)
	cd low; make -f Makefile.low extract; cd ..;
	cd gm; make -f Makefile.gm extract; cd ..;
	cd numerics; make -f Makefile.numerics extract; cd ..;
	cd graphics; make -f Makefile.graphics extract; cd ..;
	cd ui; make -f Makefile.ui extract; cd ..;

xmc:
	cd low; make -f Makefile.low xmc; cd ..;
	cd gm; make -f Makefile.gm xmc; cd ..;
	cd numerics; make -f Makefile.numerics xmc; cd ..;
	cd graphics; make -f Makefile.graphics xmc; cd ..;
	cd ui; make -f Makefile.ui xmc; cd ..;
