##############################################################################
#																			 #
#	Makefile for parallel ug version										 #
#																			 #
#	created 17 August 1992													 #
#																			 #
##############################################################################


# ARCHDIR has to be one of: INDIGO, PARAGON, PARIX, PVM, KSR
include ../ug.arch

# architecture-dependent makefile entries
include ../$(ARCHDIR)/mk.$(ARCHDIR)

OBJECTS = GenerateRules.o simplex.o 

../lib/libug.a : $(OBJECTS) Makefile
	$(LINK) $(OBJECTS) $(LLDFLAGS) -o gr $(LIBS)

.c.o:
	$(COMPILE) $(LCFLAGS) -c $<

clean:
	rm gr

