##############################################################################
#                                                                            #
#   Makefile for the coupling numproc library                                #
#                                                                            #
#   created Sep 10, 2002                                                     #
#                                                                            #
##############################################################################

# include configuration and all makefile macro definitions
include $(UGROOT)/ug.conf

.SUFFIXES: .cc

LIBRARY	= libgcoupl$(UG_LIBSUFFIX).a

OBJECTS = partdisc.o globdisc.o sparsepattern.o multivd.o coupleddiscinit.o

HEADERS = partdisc.hh globdisc.h sparsepattern.h multivd.h coupleddiscinit.h

$(LIBRARY): $(OBJECTS)
	$(ARCH_AR) $(ARCH_ARFLAGS) $(LIBRARY) $(OBJECTS)

.c.o:
	$(ARCH_CC) $(UG_CFLAGS) $<

.cc.o:
	$(ARCH_C++) $(UG_CFLAGS) $< 

$(OBJECTS): $(HEADERS)

clean:
	rm -f *.o $(LIBRARY)

# End of File
