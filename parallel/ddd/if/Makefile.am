# $Id$

SUBDIRS =

if PARALLEL
  PLIB = libif2.la libif3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS = if.h

EXTRA_DIST = ifstd.ct ifcmd.ct

SOURCES = ifcreate.c ifuse.c ifcmds.c ifcheck.c ifobjsc.c

libif2_la_SOURCES = $(SOURCES)
libif2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libif3_la_SOURCES = $(SOURCES)
libif3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
