# $Id$

# !!! missing sourcefiles...
SOURCES = idamax.c dasum.c daxpy.c dcopy.c ddot.c dnrm2.c drot.c dscal.c dgemv.c dsymv.c dtrsv.c dger.c dsyr2.c dmyblas2.c

if BUILD2D
  LIB_2D = libcblas2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libcblas3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libcblas2_la_SOURCES = $(SOURCES)
libcblas2_la_CPPFLAGS = $(UG_2DFLAGS) -DAdd_

libcblas3_la_SOURCES = $(SOURCES)
libcblas3_la_CPPFLAGS = $(UG_3DFLAGS) -DAdd_
