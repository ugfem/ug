##############################################################################
#                                                                            #
#   Makefile for famg              (UNIX version)                            #
#                                                                            #
#   created Aug 1997                                                         #
#                                                                            #
##############################################################################

FAMGLIB = libfamg.a
FAMGLINKLIB = -L. -lfamg 

APPL = famg

# SGI
#CC = aCC 
#CFLAGS = -c -n32 -mips4 -r10000 -O -fullwarn
#LN=$(CC) 
#LNFLAGS = -n32 -mips4 -r10000 -O -fullwarn
#LIBS= -lm -lc

# HP 10.20
#CC = aCC is more accurate but much slower
CC = aCC 
CFLAGS = -c -g
LN=$(CC) 
LNFLAGS = -g
LIBS= -lm -lc

ARCH_AR       = ar
ARCH_ARFLAGS  = ru
#RANLIB = ranlib
RANLIB = echo

OBJECTS = famg_heap.o famg_misc.o famg_matrix.o famg_transfer.o \
		  famg_decomp.o famg_graph.o famg_approx.o famg_construct.o \
		  famg_grid.o famg_multigrid.o famg_system.o famg_interface.o

APPLOBJ = famg.o

$(APPL)	: $(APPLOBJ) $(OBJECTS) $(FAMGLIB) makefile
	$(LN) $(LNFLAGS) $(APPLOBJ) $(FAMGLIB) $(LIBS) -o $(APPL)

$(FAMGLIB) : $(OBJECTS) makefile
	$(ARCH_AR) $(ARCH_ARFLAGS) $(FAMGLIB) $(OBJECTS) 
	$(RANLIB) $(FAMGLIB)

.C.o:		  
	$(CC) $(CFLAGS) $<

clean:
	rm -f *.o
	rm -f $(APPL)

new:
	make clean
	make
