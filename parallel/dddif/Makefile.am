noinst_LTLIBRARIES = libdddif2.la libdddif3.la

SOURCES = initddd.c memmgr.c handler.c lb.c debugger.c lbrcb.c gridcons.c support.c \
          identify.c identify.h trans.c priority.c pgmcheck.c overlap.c partition.c pamg.c debugger.h

dddifincludedir = $(pkgincludedir)
dddifinclude_HEADERS = parallel.h

libdddif2_la_SOURCES = $(SOURCES)
libdddif2_la_CPPFLAGS = $(UG_2DFLAGS)

libdddif3_la_SOURCES = $(SOURCES)
libdddif3_la_CPPFLAGS = $(UG_3DFLAGS)
