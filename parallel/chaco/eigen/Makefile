##############################################################################
#                                                                            #
#   Makefile for parallel ug version                                         #
#                                                                            #
#   created 10 Jan 1994                                                      #
#                                                                            #
##############################################################################


# ARCHDIR has to be one of: INDIGO, PARAGON, PARIX, PVM, KSR
include ../../../ug.conf

FILES.c= bisect.c checkeig.c checkorth.c  \
		eigensolve.c get_ritzvals.c lanczos_IO.c \
		lanczos_FO.c lanczos_SO.c lanczos_SO_n.c \
		lanpause.c makeorthlnk.c mkscanlist.c \
		orthogonaliz.c orthogvec.c orthog1.c \
		ql.c scanmax.c scanmin.c \
		solistout.c sorthog.c splarax.c \
		sturmcnt.c Tevec.c

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

