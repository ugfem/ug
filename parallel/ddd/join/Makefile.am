SUBDIRS =

if PARALLEL
  PLIB = libjoin2.la libjoin3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS = join.h

SOURCES = jcmds.cc join.cc

libjoin2_la_SOURCES = $(SOURCES)
libjoin2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libjoin3_la_SOURCES = $(SOURCES)
libjoin3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
