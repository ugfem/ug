##############################################################################
#																			 #
#	Makefile for netgen								 #
#																			 #
#       created by Joachim Schoeberl												 #
#																			 #
##############################################################################

# load architecture dependent makefile	
include ../../../../ug.conf

default: $(objects) 
	$(ARCH_AR) r  ../../../../lib/libnetgen.a $(objects)

.PRECIOUS:   .cc
# default rule
.cc.o:
	$(ARCH_C++) $(C++FLAGS) -I../general -I.. $<

# clean up
clean:
	rm -f *.o

