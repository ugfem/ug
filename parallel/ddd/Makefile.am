# $Id$

SUBDIRS = analyser basic ctrl f77 ident if join mgr prio xfer

if PARALLEL

if BUILD2D
  LIB_2D = libddd2.la
else
  LIB_2D =
endif
if BUILD3D
  LIB_3D = libddd3.la
else
  LIB_3D =
endif

else
  LIB_2D =
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

DDDFLAGS = -DDDD_FRONTEND_C

DDD = ddd.c dddi.h

# header files that are not explicitly used but need to be distributed
pkginclude_HEADERS = include/ddd.h

EXTRA_DIST = include/dddaddon.h include/ddd.h include/memmgr.h include/dddio.h \
             include/fddd.f include/dddresources.h \
             dddstr.h

libddd2_la_SOURCES = $(DDD)
libddd2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)
libddd2_la_LIBADD = analyser/libanalyser2.la basic/libbasic2.la ctrl/libctrl2.la f77/libf772.la ident/libident2.la if/libif2.la join/libjoin2.la mgr/libmgr2.la prio/libprio2.la xfer/libxfer2.la

libddd3_la_SOURCES = $(DDD)
libddd3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
libddd3_la_LIBADD = analyser/libanalyser3.la basic/libbasic3.la ctrl/libctrl3.la f77/libf773.la ident/libident3.la if/libif3.la join/libjoin3.la mgr/libmgr3.la prio/libprio3.la xfer/libxfer3.la

