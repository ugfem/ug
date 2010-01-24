# $Id$

SOURCES = std_domain.c std_parallel.c std_domain.h std_internal.h

if BUILD2D
  LIB_2D = libdomS2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libdomS3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

include_HEADERS = std_domain.h

libdomS2_la_SOURCES = $(SOURCES) domains2d.c
libdomS2_la_CPPFLAGS = $(UG_2DFLAGS)

libdomS3_la_SOURCES = $(SOURCES) domains3d.c
libdomS3_la_CPPFLAGS = $(UG_3DFLAGS)
