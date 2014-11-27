# $Id$

SUBDIRS =

if PARALLEL
  PLIB = libxfer2.la libxfer3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS = sll.h xfer.h

EXTRA_DIST = sll.cc sll.h

SOURCES = cmds.cc \
          cmdmsg.cc \
          cplmsg.cc \
          ctrl.cc \
          pack.cc \
          supp.cc \
          unpack.cc \
          xfer.ccc

libxfer2_la_SOURCES = $(SOURCES)
libxfer2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS) -DXferMemFromHeap

libxfer3_la_SOURCES = $(SOURCES)
libxfer3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS) -DXferMemFromHeap
