# $Id$

if PARALLEL
  LIB = libddd.la
else
  LIB =
endif

noinst_LTLIBRARIES = $(LIB)

# set options from Makefile.dddobj

DDDFLAGS = -DXferMemFromHeap -DConsMemFromHeap -DDDD_FRONTEND_C

# reproduce the style of Makefile.ddd

BASIC = basic/io.c basic/topo.c basic/notify.c basic/reduct.c basic/lowcomm.c

CTRL = ctrl/debug.c ctrl/cons.c ctrl/stat.c

MGR = mgr/typemgr.c mgr/objmgr.c mgr/cplmgr.c mgr/prio.c mgr/attr.c

IDENT = ident/ident.c

IF = if/ifcreate.c if/ifuse.c if/ifcmds.c if/ifcheck.c if/ifobjsc.c

XFER = xfer/xfer.c xfer/cmds.c xfer/pack.c xfer/unpack.c xfer/ctrl.c \
       xfer/cplmsg.c xfer/cmdmsg.c xfer/supp.c

PRIO = prio/pcmds.c

JOIN = join/join.c join/jcmds.c

MISC = analyser/inv.c

DDD = ddd.c dddi.h

INCLUDES = include/dddaddon.h include/ddd.h include/memmgr.h include/dddio.h \
           include/fddd.f include/dddresources.h

libddd_la_SOURCES = $(BASIC) $(CTRL) $(MGR) $(IDENT) $(IF) $(XFER) \
                    $(PRIO) $(JOIN) $(MISC) $(DDD) f77/f2c.c
libddd_la_CPPFLAGS = $(UG_CPPFLAGS) $(DDDFLAGS)

EXTRA_DIST = $(INCLUDES)
