# $Id$

# Build the following subdirectories.  The chaco subdir
# is not built because it is not available in the anonymous
# Heidelberg CVS and not in the set of patches we provide
# either.  AFAIK it doesn't work anyways.
SUBDIRS = dddif ddd util ppif #chaco

if PARALLEL
  PLIB = libparallel2.la libparallel3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

SOURCES = initparallel.c initparallel.h

SUBLIBS = ddd/libddd.la ppif/MPI/libppifmpi.la

libparallel2_la_SOURCES = $(SOURCES)
libparallel2_la_CPPFLAGS = $(UG_2DFLAGS)
libparallel2_la_LIBADD = $(SUBLIBS) dddif/libdddif2.la

libparallel3_la_SOURCES = $(SOURCES)
libparallel3_la_CPPFLAGS = $(UG_3DFLAGS)
libparallel3_la_LIBADD = $(SUBLIBS) dddif/libdddif3.la
