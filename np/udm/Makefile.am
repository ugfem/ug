SOURCES = \
  pcr.cc \
  disctools.cc \
  npscan.cc \
  udm.cc \
  formats.cc \
  numproc.cc \
  dio.cc \
  data_io.cc \
  data_io.h \
  dio.h \
  disctools.h \
  formats.h \
  npscan.h \
  numproc.h \
  pcr.h \
  udm.h

if BUILD2D
  LIB_2D = libudm2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libudm3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

npincludedir = $(pkgincludedir)
npinclude_HEADERS = udm.h numproc.h formats.h npscan.h disctools.h pcr.h

libudm2_la_SOURCES = $(SOURCES)
libudm2_la_CPPFLAGS = $(UG_2DFLAGS)

libudm3_la_SOURCES = $(SOURCES)
libudm3_la_CPPFLAGS = $(UG_3DFLAGS)
