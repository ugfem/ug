SUBDIRS = uggraph grape

SOURCES = graphics.c graphics.h

if BUILD2D
  LIB_2D = libgraphics2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libgraphics3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libgraphics2_la_SOURCES = $(SOURCES)
libgraphics2_la_CPPFLAGS = $(UG_2DFLAGS)
libgraphics2_la_LIBADD = uggraph/libuggraph2.la \
                         grape/libuggrape2.la

libgraphics3_la_SOURCES = $(SOURCES)
libgraphics3_la_CPPFLAGS = $(UG_3DFLAGS)
libgraphics3_la_LIBADD = uggraph/libuggraph3.la \
                         grape/libuggrape3.la

graphicsincludedir = $(pkgincludedir)
graphicsinclude_HEADERS = graphics.h

