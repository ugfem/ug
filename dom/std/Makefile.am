# $Id$

SOURCES = std_domain.c std_parallel.c std_domain.h std_internal.h

EXTRALIBS = ../../low/liblow.la ../../dev/libdevices.la
if BUILD2D
  LIB_2D = libdomS2.la
if PARALLEL
  EXTRALIBS += ../../parallel/libparallel2.la ../../gm/libug_gm2.la
endif
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libdomS3.la
if PARALLEL
  EXTRALIBS += ../../parallel/libparallel3.la ../../gm/libug_gm3.la
endif
else
  LIB_3D =
endif

lib_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

include_HEADERS = std_domain.h

libdomS2_la_SOURCES = $(SOURCES) domains2d.c
libdomS2_la_CPPFLAGS = $(UG_2DFLAGS)
libdomS2_la_LIBADD = $(EXTRALIBS)

libdomS3_la_SOURCES = $(SOURCES) repair.c domains3d.c
libdomS3_la_CPPFLAGS = $(UG_3DFLAGS)
libdomS3_la_LIBADD = $(EXTRALIBS)
