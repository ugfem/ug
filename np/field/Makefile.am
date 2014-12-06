SOURCES = stoch.c field.h

if BUILD2D
  LIB_2D = libfield2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libfield3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libfield2_la_SOURCES = $(SOURCES)
libfield2_la_CPPFLAGS = $(UG_2DFLAGS)

libfield3_la_SOURCES = $(SOURCES)
libfield3_la_CPPFLAGS = $(UG_3DFLAGS)
