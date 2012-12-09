# $Id$

SUBDIRS = test

if PARALLEL
  PLIB = libbasic2.la libbasic3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

# header files that are not explicitly used but need to be distributed
EXTRA_DIST = basic/notify.h basic/lowcomm.h basic/oopp.h basic/ooppcc.h basic/memutil.h

SOURCES = io.c lowcomm.c notify.c reduct.c topo.c

libbasic2_la_SOURCES = $(SOURCES)
libbasic2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libbasic3_la_SOURCES = $(SOURCES)
libbasic3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
