# $Id$

SUBDIRS = gg2 gg3

SOURCES = algebra.c enrol.c evm.c mgio.c ugio.c ugm.c cw.c initgm.c elements.c shapes.c evalproc.c rm.c refine.c smooth.c dlmgr.c gmcheck.c er.c mgheapmgr.c algebra.h cw.h dlmgr.h elements.h enrol.h er.h evm.h GenerateRules.h gm.h initgm.h mgheapmgr.h mgio.h pargm.h refine.h rm.h shapes.h ugio.h ugm.h

# put this "template" into distribution-tarball as well
EXTRA_DIST = dlmgr.t

gmincludedir = $(pkgincludedir)
gminclude_HEADERS = gm.h pargm.h evm.h shapes.h dlmgr.h algebra.h \
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
