##############################################################################
#                                                                            #
#   Makefile for parallel ug version                                         #
#                                                                            #
#   created 10 Jan 1994                                                      #
#                                                                            #
##############################################################################


# ARCHDIR has to be one of: INDIGO, PARAGON, PARIX, PVM, KSR
include ../../../ug.conf

# main/main.c is not compiled and linked */
# util/smalloc.c is not compiled and linked */
# util/nrerror.c util/randomize.c util/smalloc.c util/scadd.c 
FILES.c= checkbp.c inits2d.c inits3d.c \
		genvals2d.c genvals3d.c map2d.c \
		map3d.c mergesort.c movevtxs.c \
		sorts2d.c sorts3d.c

FILES.o=	$(FILES.c:.c=.o) 

# local C compiler flags
LCFLAGS = -I../../include

dummy:	${FILES.o} Makefile

lint:
		lint ${FILES.c} -lm

alint:
		alint ${FILES.c} -lm

clean:
		rm -f *.o

.c.o:
		${ARCH_CC} $(UG_CFLAGS) ${LCFLAGS} -c $<

