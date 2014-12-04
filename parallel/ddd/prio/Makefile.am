# $Id$

SUBDIRS =

if PARALLEL
  PLIB = libprio2.la libprio3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS =

SOURCES = pcmds.cc

libprio2_la_SOURCES = $(SOURCES)
libprio2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libprio3_la_SOURCES = $(SOURCES)
libprio3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
