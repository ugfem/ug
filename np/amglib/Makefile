##############################################################################
#                                                                            #
#   Makefile for the amg library (UNIX version)                              #
#   (stand alone version)                                                    #
#                                                                            #
#   created  08 OKT 1997                                                     #
#                                                                            #
##############################################################################

COMPILE   = cc
CC        = $(COMPILE)
LINK      = cc

CFLAGS    = -D__SGI__ -g -mips2 -c
LLDFLAGS  = -mips2
LIBS      = -lm

# set application name
APPL        = test

OBJECTS  = amg_low.o amg_sp.o amg_blas.o amg_iter.o amg_coarsen.o amg_solve.o \
           test.o


$(APPL) : $(OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(APPL) $(OBJECTS) $(LIBS)
 
.c.o:
	$(COMPILE) $(CFLAGS) $<

clean:
	rm *.o
	rm $(APPL)
	
new:
	rm $(APPL)
	make
