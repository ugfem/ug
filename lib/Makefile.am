# $Id$

SUBDIRS = ugdata

if BUILD2D
  LIB_2D = libugS2.la libugL2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libugS3.la libugL3.la
else
  LIB_3D =
endif

2DLIBS = ../gm/libug_gm2.la \
         ../np/libnp2.la \
         ../graphics/libgraphics2.la \
         ../ui/libugui2.la \
         ../dev/libdevices.la \
         ../low/liblow2.la

3DLIBS = ../gm/libug_gm3.la \
         ../np/libnp3.la \
         ../graphics/libgraphics3.la \
         ../ui/libugui3.la \
         ../dev/libdevices.la \
         ../low/liblow3.la

lib_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libugS2_la_SOURCES = ../initug.c
libugS2_la_CPPFLAGS = $(UG_2DFLAGS)
libugS2_la_LIBADD = $(2DLIBS) ../dom/std/libdomS2.la
if PARALLEL
  libugS2_la_LIBADD += ../parallel/libparallel2.la
endif

libugL2_la_SOURCES = ../initug.c
libugL2_la_CPPFLAGS = $(UG_2DFLAGS)
libugL2_la_LIBADD = $(2DLIBS) ../dom/lgm/libdomL2.la

libugS3_la_SOURCES = ../initug.c
libugS3_la_CPPFLAGS = $(UG_3DFLAGS)
libugS3_la_LIBADD = $(3DLIBS) ../dom/std/libdomS3.la
if PARALLEL
  libugS3_la_LIBADD += ../parallel/libparallel3.la
endif

libugL3_la_SOURCES = ../initug.c
libugL3_la_CPPFLAGS = $(UG_3DFLAGS)
libugL3_la_LIBADD = $(3DLIBS) ../dom/lgm/libdomL3.la
