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
all: include $(MODULES) $(OBJECTS)
	make $(UG_LIB) 	
	$(ARCH_AR) $(ARCH_ARFLAGS) $(UG_LIB) $(OBJECTS) 
	echo "libug, libdom and libdev compiled"

uglib: include $(UGMODULES)
	make $(UG_LIB) 	
	echo "libug compiled"

UGD: include $(UGDMODULES) ugd.o
	make $(UG_LIB) 	
	$(UG_LINK) -o bin/ugd $(ARCH_LFLAGS) ugd.o lib/libdev$(IF).a \
                          $(UG_LIB) $(UG_LFLAGS)
	echo "ugd compiled"

init: $(OBJECTS)
	make $(UG_LIB)

$(UG_LIB): $(OBJECTS)
	$(ARCH_AR) $(ARCH_ARFLAGS) $(UG_LIB) $(OBJECTS) 

##############################################################################

LOW: include
	cd low && make -f Makefile.low

DEV: include
	cd dev && make -f Makefile.dev

DOM: include
	cd dom && make -f Makefile.dom

GM: include
	cd gm && make -f Makefile.gm

NUMERICS: include
	cd numerics && make -f Makefile.numerics

NP: include
	cd np && make -f Makefile.np

GRAPHICS: include
	cd graphics && make -f Makefile.graphics

UI: include
	cd ui && make -f Makefile.ui


##############################################################################

# targets for MODEL, exactly one will be chosen

SEQUENTIAL:

PARALLEL: include
	cd parallel; $(ARCH_MAKE) -f Makefile.parallel


SEQUENTIAL_clean:

PARALLEL_clean:
	cd parallel; $(ARCH_MAKE) -f Makefile.parallel clean


##############################################################################

include:
	ugmakelinks;

.c.o:
	$(ARCH_CC) $(UG_CFLAGS) $(LCFLAGS) $<


##############################################################################

clean: $(MODEL_TARGET)_clean
	rm -f $(OBJECTS)
	cd low && make -f Makefile.low clean
	cd dev && make -f Makefile.dev clean
	cd dom && make -f Makefile.dom clean
	cd gm && make -f Makefile.gm clean
	cd np && make -f Makefile.np clean
	cd graphics && make -f Makefile.graphics clean
	cd ui && make -f Makefile.ui clean

ar: 
	cd gm && make -f Makefile.gm ar
	cd graphics && make -f Makefile.graphics ar
	cd ui && make -f Makefile.ui ar
	cd np && make -f Makefile.np ar
	$(ARCH_AR) $(ARCH_ARFLAGS) $(UG_LIB) $(OBJECTS);
	make;

ifdef: $(MODEL_TARGET)_clean
	cd gm && make -f Makefile.gm clean
	cd np && make -f Makefile.np clean
	cd graphics && make -f Makefile.graphics clean
	cd dom && make -f Makefile.dom clean
	cd ui && make -f Makefile.ui clean
	rm -f initug.o;

extract:
	$(ARCH_AR) $(ARCH_EXFLAGS) lib/libug$(UG_LIBSUFFIX).a $(OBJECTS)
	cd low && make -f Makefile.low extract
	cd dev && make -f Makefile.dev extract
	cd dom && make -f Makefile.dom extract
	cd gm && make -f Makefile.gm extract
	cd np && make -f Makefile.np extract
	cd graphics && make -f Makefile.graphics extract
	cd ui && make -f Makefile.ui extract

xmc:
	cd low && make -f Makefile.low xmc
	cd dev && make -f Makefile.dev xmc
	cd dom && make -f Makefile.dom xmc
	cd gm && make -f Makefile.gm xmc
	cd np && make -f Makefile.np xmc
	cd graphics && make -f Makefile.graphics xmc
	cd ui && make -f Makefile.ui xmc
