# $Id$

if BUILD2D
  LIB_2D = libug2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libug3.la
else
  LIB_3D =
endif

if PARALLEL
  PLIBS_2D = ../parallel/libparallel2.la
  PLIBS_3D = ../parallel/libparallel3.la
else
  PLIBS_2D =
  PLIBS_3D =
endif

lib_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

# seems to miss...
#                   ../parallel/util/libparutil.la

libug2_la_SOURCES = ../initug.c
libug2_la_CPPFLAGS = $(UG_2DFLAGS)
libug2_la_LIBADD = ../gm/libug_gm2.la \
                   ../np/libnp2.la \
                   ../graphics/libgraphics2.la \
                   ../ui/libugui2.la \
                   ../dev/libdevices.la \
                   ../low/liblow2.la \
                   ../low/liblow.la $(PLIBS_2D)

libug3_la_SOURCES = ../initug.c
libug3_la_CPPFLAGS = $(UG_3DFLAGS)
libug3_la_LIBADD = ../gm/libug_gm3.la \
                   ../np/libnp3.la \
                   ../graphics/libgraphics3.la \
                   ../ui/libugui3.la \
                   ../dev/libdevices.la \
                   ../low/liblow3.la \
                   ../low/liblow.la $(PLIBS_3D)
