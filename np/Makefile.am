# $Id$

SUBDIRS = algebra amglib field procs udm slu

# famglib ist auch im alten Makefile nicht drin

include_HEADERS = np.h

SOURCES = initnp.c initnp.h initnumerics.h np.h num.h

if BUILD2D
  LIB_2D = libnp2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libnp3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libnp2_la_SOURCES = $(SOURCES)
libnp2_la_CPPFLAGS = $(UG_2DFLAGS)
libnp2_la_LIBADD = algebra/libalgebra2.la \
                   amglib/libamg2.la \
                   field/libfield2.la \
                   procs/libprocs2.la \
                   slu/libslu2.la \
                   udm/libudm2.la \
                   ../low/liblow2.la

libnp3_la_SOURCES = $(SOURCES)
libnp3_la_CPPFLAGS = $(UG_3DFLAGS)
libnp3_la_LIBADD = algebra/libalgebra3.la \
                   amglib/libamg3.la \
                   field/libfield3.la \
                   procs/libprocs3.la \
                   slu/libslu3.la \
                   udm/libudm3.la \
                   ../low/liblow3.la

