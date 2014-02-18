# $Id$

SUBDIRS = basic mgr xfer

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

# The first one is only used in ctrl
DDDFLAGS = -DConsMemFromHeap -DDDD_FRONTEND_C

# reproduce the style of Makefile.ddd

CTRL = ctrl/debug.c ctrl/cons.c ctrl/stat.c

IDENT = ident/ident.c

IF = if/ifcreate.c if/ifuse.c if/ifcmds.c if/ifcheck.c if/ifobjsc.c

PRIO = prio/pcmds.c

JOIN = join/join.c join/jcmds.c

MISC = analyser/inv.c

DDD = ddd.c dddi.h

# header files that are not explicitly used but need to be distributed
pkginclude_HEADERS = include/ddd.h

EXTRA_DIST = include/dddaddon.h include/ddd.h include/memmgr.h include/dddio.h \
             include/fddd.f include/dddresources.h join/join.h ctrl/stat.h \
             xif/if.h if/ifstd.ct if/ifcmd.ct \
             dddstr.h

libddd2_la_SOURCES = $(CTRL) $(IDENT) $(IF) $(XFER) \
                    $(PRIO) $(JOIN) $(MISC) $(DDD) f77/f2c.c
libddd2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)
libddd2_la_LIBADD = basic/libbasic2.la mgr/libmgr2.la xfer/libxfer2.la

libddd3_la_SOURCES = $(CTRL) $(IDENT) $(IF) $(XFER) \
                    $(PRIO) $(JOIN) $(MISC) $(DDD) f77/f2c.c
libddd3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
libddd3_la_LIBADD = basic/libbasic3.la mgr/libmgr3.la xfer/libxfer3.la

