if BUILD2D
  LIB_2D = libuggrape2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libuggrape3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

# !!! sollte mit "if GRAPE" die richtigen Sourcen nehmen
SOURCES = connectuggrapeOFF.cc connectuggrape.h defs.h

libuggrape2_la_SOURCES = $(SOURCES)
libuggrape2_la_CPPFLAGS = $(UG_2DFLAGS)

libuggrape3_la_SOURCES = $(SOURCES)
libuggrape3_la_CPPFLAGS = $(UG_3DFLAGS)
