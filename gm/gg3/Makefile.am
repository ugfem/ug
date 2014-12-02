if BUILD3D
  LIB_3D = libgg3.la
else
  LIB_3D =
endif

# this lib should be directly installed
noinst_LTLIBRARIES = libgg3.la

libgg3_la_SOURCES = gg3d.cc gg3d.h gginterface.h
libgg3_la_CPPFLAGS = $(UG_3DFLAGS)
