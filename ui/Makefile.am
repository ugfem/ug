SOURCES = uginterface.cc initui.cc cmdint.cc cmdline.cc \
    dataexplorer.cc commands.cc helpmsg.cc avs.cc \
    tecplot.cc mmio.cc fieldio.cc \
    avs.hh cmdint.hh cmdline.hh commands.hh dataexplorer.hh \
    fieldio.hh helpmsg.hh initui.hh mmio.hh tecplot.hh uginterface.hh

uiincludedir = $(pkgincludedir)
uiinclude_HEADERS = cmdint.hh cmdline.hh commands.hh

if BUILD2D
  LIB_2D = libugui2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libugui3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libugui2_la_SOURCES = $(SOURCES)
libugui2_la_CPPFLAGS = $(UG_2DFLAGS)

libugui3_la_SOURCES = $(SOURCES)
libugui3_la_CPPFLAGS = $(UG_3DFLAGS)
