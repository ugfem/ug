# $Id$

SOURCES = ugblas.c ugiter.c transgrid.c block.c quadrature.c fvgeom.c fegeom.c ff_gen.c ff.c amgtools.c npcheck.c sm.c blasm.c ugeblas.c amgtools.h blasm.h block.h fegeom.h ff_gen.h ff.h fvgeom.h npcheck.h quadrature.h sm.h transgrid.h ugblas.h

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
EXTRA_DIST = matfunc.ct matloop.ct matmode.ct vecfunc.ct vecloop.ct vecmode.ct

include_HEADERS = sm.h fvgeom.h fegeom.h quadrature.h block.h

libalgebra2_la_SOURCES = $(SOURCES)
libalgebra2_la_CPPFLAGS = $(UG_2DFLAGS)

libalgebra3_la_SOURCES = $(SOURCES)
libalgebra3_la_CPPFLAGS = $(UG_3DFLAGS)
