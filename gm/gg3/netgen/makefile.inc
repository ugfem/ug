##############################################################################
#																			 #
#	Makefile for netgen								 #
#																			 #
#       created by Joachim Schoeberl												 #
#																			 #
##############################################################################

# load architecture dependent makefile	
include ../../../ug.arch

# architecture-dependent makefile entries
include ../../../machines/mk.$(ARCHDIR)

# load configuration switches
include ../../../ug.conf


default: $(objects) 
	ar r ../../../lib/libnetgen.a $(objects)

# default rule
.cc.o:
	$(C++) $(C++FLAGS) -I../general -I.. $<

# clean up
clean:
	rm *.o

