##############################################################################
#                                                                            #
#   Makefile for parallel ug version                                         #
#                                                                            #
#   created 10 Jan 1994                                                      #
#                                                                            #
##############################################################################


# ARCHDIR has to be one of: INDIGO, PARAGON, PARIX, PVM, KSR
include ../../../ug.conf

FILES.c= affirm.c array_alloc2.c bit_reverse.c \
		checkpnt.c cpvec.c dot.c doubleout.c \
		input_int.c gray.c machine_para.c \
		makevwsqrt.c mkvec.c norm.c normalize.c \
		nrerror.c randomize.c smalloc.c scadd.c \
		scale.c seconds.c setvec.c shell_sort.c \
		tri_prod.c true_or_fals.c update.c \
		vecout.c vecran.c vecscale.c

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

