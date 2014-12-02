if BUILD2D
  LIB_2D = libdomL2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libdomL3.la
else
  LIB_3D =
endif

# !!! properly integrate netgen
SUBDIRS = ngin ngin2d

# !!! not good...
EXTRA_DIST = netgen/lgm_gginterface.h

SOURCES = lgm_domain.cc lgm_load.cc lgm_parallel.cc ansys2lgm.h lgm_domain.h lgm_load.h lgm_macros.h lgm_transfer.h

domincludedir = $(pkgincludedir)
dominclude_HEADERS = lgm_domain.h

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libdomL2_la_SOURCES = $(SOURCES) lgm_transfer2d.cc lgm_domain2d.cc
libdomL2_la_CPPFLAGS = $(UG_2DFLAGS)
libdomL2_la_LIBADD = ngin2d/libngin2d.la

libdomL3_la_SOURCES = $(SOURCES) lgm_transfer3d.cc lgm_domain3d.cc ansys2lgm.cc lgm_accel.cc lgm_accel.h
libdomL3_la_CPPFLAGS = $(UG_3DFLAGS)
libdomL3_la_LIBADD = ngin/libngin.la
