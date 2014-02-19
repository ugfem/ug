# $Id$

SUBDIRS =

if PARALLEL
  PLIB = libanalyser2.la libanalyser3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS =

SOURCES = inv.c

libanalyser2_la_SOURCES = $(SOURCES)
libanalyser2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libanalyser3_la_SOURCES = $(SOURCES)
libanalyser3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
