# makefile for sparse supernodal LU, implemented in ANSI C
include ../make.inc

#######################################################################
#  This is the makefile to create a library for SuperLU.
#  The files are organized as follows:
#
#       ALLAUX -- Auxiliary routines called from all precisions
#       SCLAUX -- Auxiliary routines called from both real and complex
#       DZLAUX -- Auxiliary routines called from both double precision
#                 and complex*16
#	SLUSRC -- Single precision real SuperLU routines
#       DLUSRC -- Double precision real SuperLU routines
#	CLUSRC -- Single precision complex SuperLU routines
#       ZLUSRC -- Double precision complex SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single double
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################

ALLAUX = superlu_timer.o lsame.o util.o memory.o get_perm_c.o mmd.o \
	 sp_coletree.o sp_preorder.o sp_ienv.o relax_snode.o \
	 xerbla.o colamd.o

DZLAUX = dlamch.o

DLUSRC = \
	dgssv.o dgssvx.o \
	dsp_blas2.o dsp_blas3.o dgscon.o dlacon.o \
	dlangs.o dgsequ.o dlaqgs.o dpivotgrowth.o  \
	dgsrfs.o dgstrf.o dgstrs.o dcopy_to_ucol.o \
	dsnode_dfs.o dsnode_bmod.o \
	dpanel_dfs.o dpanel_bmod.o dreadhb.o \
	dcolumn_dfs.o dcolumn_bmod.o dpivotL.o dpruneL.o \
	dmemory.o dutil.o

all:    double 

double: $(DLUSRC) $(ALLAUX) $(DZLAUX)
	$(ARCH) $(ARCHFLAGS) $(SUPERLULIB) $(DLUSRC) $(ALLAUX) $(DZLAUX)
	$(RANLIB) $(SUPERLULIB)

##################################
# Do not optimize this routine   #
##################################
dlamch.o: dlamch.c ; $(CC) -c $(NOOPTS) $<

superlu_timer.o:  superlu_timer.c ; $(CC) -c $(NOOPTS) $<

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) -c $< $(VERBOSE)


