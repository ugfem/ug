SUBDIRS = test

if PARALLEL
  PLIB = libbasic2.la libbasic3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

pkginclude_HEADERS = notify.h lowcomm.h oopp.h ooppcc.h memutil.h

SOURCES = io.cc lowcomm.cc notify.cc reduct.cc topo.cc

libbasic2_la_SOURCES = $(SOURCES)
libbasic2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libbasic3_la_SOURCES = $(SOURCES)
libbasic3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
