noinst_LTLIBRARIES = libdddif2.la libdddif3.la

SOURCES = initddd.cc memmgr.cc handler.cc lb.cc debugger.cc lbrcb.cc gridcons.cc support.cc \
          identify.cc identify.h trans.cc priority.cc pgmcheck.cc overlap.cc partition.cc pamg.cc debugger.h

dddifincludedir = $(pkgincludedir)
dddifinclude_HEADERS = parallel.h

libdddif2_la_SOURCES = $(SOURCES)
libdddif2_la_CPPFLAGS = $(UG_2DFLAGS)

libdddif3_la_SOURCES = $(SOURCES)
libdddif3_la_CPPFLAGS = $(UG_3DFLAGS)
