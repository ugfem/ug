SUBDIRS = dddif ddd util ppif

if PARALLEL
  PLIB = libparallel2.la libparallel3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

SOURCES = initparallel.cc initparallel.h

SUBLIBS = ppif/MPI/libppifmpi.la util/libparutil.la

libparallel2_la_SOURCES = $(SOURCES)
libparallel2_la_CPPFLAGS = $(UG_2DFLAGS)
libparallel2_la_LIBADD = $(SUBLIBS) ddd/libddd2.la dddif/libdddif2.la

libparallel3_la_SOURCES = $(SOURCES)
libparallel3_la_CPPFLAGS = $(UG_3DFLAGS)
libparallel3_la_LIBADD = $(SUBLIBS) ddd/libddd3.la dddif/libdddif3.la
