# $Id$

SUBDIRS = cblas src

SOURCES = slu.c slu.h

if BUILD2D
  LIB_2D = libslu2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libslu3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libslu2_la_SOURCES = $(SOURCES)
libslu2_la_CPPFLAGS = $(UG_2DFLAGS) -DAdd_
libslu2_la_LIBADD = cblas/libcblas2.la \
                    src/librealslu2.la

libslu3_la_SOURCES = $(SOURCES)
libslu3_la_CPPFLAGS = $(UG_3DFLAGS) -DAdd_
libslu3_la_LIBADD = cblas/libcblas3.la \
                    src/librealslu3.la
