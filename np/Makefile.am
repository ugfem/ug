SUBDIRS = algebra amglib field procs udm

# famglib ist auch im alten Makefile nicht drin

npincludedir = $(pkgincludedir)
npinclude_HEADERS = np.h

SOURCES = initnp.cc initnp.h initnumerics.h np.h num.h

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
                   udm/libudm2.la

libnp3_la_SOURCES = $(SOURCES)
libnp3_la_CPPFLAGS = $(UG_3DFLAGS)
libnp3_la_LIBADD = algebra/libalgebra3.la \
                   amglib/libamg3.la \
                   field/libfield3.la \
                   procs/libprocs3.la \
                   udm/libudm3.la

