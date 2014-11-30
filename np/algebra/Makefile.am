SOURCES = \
  ugblas.cc \
  ugiter.cc \
  transgrid.cc \
  block.cc \
  quadrature.cc \
  fvgeom.cc \
  fegeom.cc \
  ff_gen.cc \
  ff.cc \
  amgtools.cc \
  npcheck.cc \
  sm.cc \
  blasm.cc \
  ugeblas.cc \
  amgtools.h \
  blasm.h \
  block.h \
  fegeom.h \
  ff_gen.h \
  ff.h \
  fvgeom.h \
  np.h \
  eck.h \
  quadrature.h \
  sm.h \
  transgrid.h \
  ugblas.h

if BUILD2D
  LIB_2D = libalgebra2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libalgebra3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

# include templates in distribution tarball
EXTRA_DIST = \
  matfunc.ct \
  matloop.ct \
  matmode.ct \
  vecfunc.ct \
  vecloop.ct \
  vecmode.ct

npincludedir = $(pkgincludedir)
npinclude_HEADERS = \
  block.h \
  fegeom.h \
  fvgeom.h \
  quadrature.h \
  sm.h

libalgebra2_la_SOURCES = $(SOURCES)
libalgebra2_la_CPPFLAGS = $(UG_2DFLAGS)

libalgebra3_la_SOURCES = $(SOURCES)
libalgebra3_la_CPPFLAGS = $(UG_3DFLAGS)
