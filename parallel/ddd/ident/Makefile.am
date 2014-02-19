# $Id$

SUBDIRS =

if PARALLEL
  PLIB = libident2.la libident3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS =

SOURCES = ident.c

libident2_la_SOURCES = $(SOURCES)
libident2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libident3_la_SOURCES = $(SOURCES)
libident3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
