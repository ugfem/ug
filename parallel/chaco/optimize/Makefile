##############################################################################
#                                                                            #
#   Makefile for parallel ug version                                         #
#                                                                            #
#   created 10 Jan 1994                                                      #
#                                                                            #
##############################################################################


# ARCHDIR has to be one of: INDIGO, PARAGON, PARIX, PVM, KSR
include ../../../ug.conf

FILES.c= determinant.c func2d.c \
		func3d.c opt2d.c opt3d.c

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

