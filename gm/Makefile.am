# $Id$

SUBDIRS = gg2 gg3

SOURCES = algebra.c enrol.c evm.c mgio.c ugio.c ugm.c cw.c initgm.c elements.c shapes.c evalproc.c rm.c refine.c smooth.c dlmgr.c gmcheck.c er.c mgheapmgr.c algebra.h dlmgr.h enrol.h er.h evm.h GenerateRules.h initgm.h mgheapmgr.h mgio.h shapes.h ugio.h

# put this "template" into distribution-tarball as well
EXTRA_DIST = dlmgr.t

gmincludedir = $(pkgincludedir)
gminclude_HEADERS = cw.h elements.h gm.h pargm.h evm.h shapes.h ugm.h dlmgr.h algebra.h \
  rm.h refine.h



if BUILD2D
  LIB_2D = libug_gm2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libug_gm3.la
else
  LIB_3D =
endif

# this library should only exist during the build (noinst)
noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libug_gm2_la_SOURCES = $(SOURCES)
libug_gm2_la_CPPFLAGS = $(UG_2DFLAGS)
libug_gm2_la_LIBADD = gg2/libgg2.la

libug_gm3_la_SOURCES = $(SOURCES)
libug_gm3_la_CPPFLAGS = $(UG_3DFLAGS)
libug_gm3_la_LIBADD = gg3/libgg3.la
