if MPI
  MLIB = libppifmpi.la
else
  MLIB =
endif

noinst_LTLIBRARIES = $(MLIB)

libppifmpi_la_SOURCES = ppif.cc
libppifmpi_la_CPPFLAGS = $(MPI_CPPFLAGS) $(UG_EXTRAFLAGS)
libppifmpi_la_LDFLAGS = $(MPI_LDFLAGS)
