SUBDIRS =

if PARALLEL
  PLIB = libctrl2.la libctrl3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS = stat.h

SOURCES = cons.c debug.c stat.c

libctrl2_la_SOURCES = $(SOURCES)
libctrl2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS) -DConsMemFromHeap

libctrl3_la_SOURCES = $(SOURCES)
libctrl3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS) -DConsMemFromHeap
