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
CC = CC 
CFLAGS = -c -w -O
LN=$(CC) 
LNFLAGS = -O
LIBS= -lm -lc

ARCH_AR       = ar
ARCH_ARFLAGS  = ru
#RANLIB = ranlib
RANLIB = echo

OBJECTS = heap.o misc.o matrix.o transfer.o decomp.o graph.o approx.o \
construct.o grid.o multigrid.o system.o interface.o

APPLOBJ = famg.o

$(APPL)	: $(APPLOBJ) $(OBJECTS) $(FAMGLIB) Makefile
	$(LN) $(LNFLAGS) $(APPLOBJ) $(FAMGLIB) $(LIBS) -o $(APPL)

$(FAMGLIB) : $(OBJECTS) Makefile
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
