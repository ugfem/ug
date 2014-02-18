# $Id$

SUBDIRS = basic ctrl if mgr xfer

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

# set options from Makefile.dddobj
DDDFLAGS = -DDDD_FRONTEND_C

# reproduce the style of Makefile.ddd

IDENT = ident/ident.c

PRIO = prio/pcmds.c

JOIN = join/join.c join/jcmds.c

MISC = analyser/inv.c

DDD = ddd.c dddi.h

# header files that are not explicitly used but need to be distributed
pkginclude_HEADERS = include/ddd.h

EXTRA_DIST = include/dddaddon.h include/ddd.h include/memmgr.h include/dddio.h \
             include/fddd.f include/dddresources.h join/join.h \
             xif/if.h dddstr.h

libddd2_la_SOURCES = $(IDENT) $(PRIO) $(JOIN) $(MISC) $(DDD) f77/f2c.c
libddd2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)
libddd2_la_LIBADD = basic/libbasic2.la ctrl/libctrl2.la if/libif2.la mgr/libmgr2.la xfer/libxfer2.la

libddd3_la_SOURCES = $(IDENT) $(PRIO) $(JOIN) $(MISC) $(DDD) f77/f2c.c
libddd3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
libddd3_la_LIBADD = basic/libbasic3.la ctrl/libctrl3.la if/libif3.la mgr/libmgr3.la xfer/libxfer3.la

