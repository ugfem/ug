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
  PLIBS = ../parallel/ppif/libppif.la ../parallel/ddd/libddd.la
else
  PLIBS =
endif

lib_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

# seems to miss...
#                   ../parallel/util/libparutil.la

libug2_la_SOURCES = ../initug.c
libug2_la_CPPFLAGS = $(UG_2DFLAGS)
libug2_la_LIBADD = ../gm/libug_gm2.la \
                   ../np/libnp2.la \
                   ../graphics/libgraphics2.la \
                   ../low/libtree2.la \
                   ../ui/libugui2.la \
                   ../dev/libdevices.la \
                   ../low/liblow.la $(PLIBS)

libug3_la_SOURCES = ../initug.c
libug3_la_CPPFLAGS = $(UG_3DFLAGS)
libug3_la_LIBADD = ../gm/libug_gm3.la \
                   ../np/libnp3.la \
                   ../graphics/libgraphics3.la \
                   ../low/libtree3.la \
                   ../ui/libugui3.la \
                   ../dev/libdevices.la \
                   ../low/liblow.la $(PLIBS)
