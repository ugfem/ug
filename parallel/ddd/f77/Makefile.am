# $Id$

SUBDIRS =

if PARALLEL
  PLIB = libf772.la libf773.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS =

SOURCES = f2c.c

libf772_la_SOURCES = $(SOURCES)
libf772_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libf773_la_SOURCES = $(SOURCES)
libf773_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
