# $Id$

SOURCES = gen.c gen_parallel.c gen.h

if BUILD2D
  LIB_2D = libdomG2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libdomG3.la
else
  LIB_3D =
endif

lib_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libdomG2_la_SOURCES = $(SOURCES)
libdomG2_la_CPPFLAGS = $(UG_2DFLAGS)

libdomG3_la_SOURCES = $(SOURCES)
libdomG3_la_CPPFLAGS = $(UG_3DFLAGS)

