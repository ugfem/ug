# $Id$

if MPI
  MLIB = libppifmpi.la
else
  MLIB =
endif

noinst_LTLIBRARIES = $(MLIB)

libppifmpi_la_SOURCES = ppif.c
libppifmpi_la_CPPFLAGS = $(MPI_CPPFLAGS)
libppifmpi_la_LDFLAGS = $(MPI_LDFLAGS)
