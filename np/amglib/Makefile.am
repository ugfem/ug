# $Id$

SOURCES = amg_blas.c amg_coarsen.c amg_iter.c amg_low.c amg_solve.c amg_sp.c amg_ug.c amg_blas.h amg_coarsen.h amg_header.h amg_iter.h amg_low.h amg_solve.h amg_sp.h amg_ug.h

if BUILD2D
  LIB_2D = libamg2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libamg3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libamg2_la_SOURCES = $(SOURCES)
libamg2_la_CPPFLAGS = $(UG_2DFLAGS)

libamg3_la_SOURCES = $(SOURCES)
libamg3_la_CPPFLAGS = $(UG_3DFLAGS)
