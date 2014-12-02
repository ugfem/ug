if BUILD2D
  LIB_2D = libgg2.la
else
  LIB_2D =
endif

# should be installed
noinst_LTLIBRARIES = $(LIB_2D)

libgg2_la_SOURCES = ggmain.cc ggaccel.cc ggm.cc ggaccel.h ggmain.h ggm.h
libgg2_la_CPPFLAGS = $(UG_2DFLAGS)
