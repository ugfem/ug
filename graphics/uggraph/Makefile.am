# $Id$

if BUILD2D
  LIB_2D = libuggraph2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libuggraph3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

SOURCES = wpm.c graph.c wop.c initgraph.c plotproc.c placer.c bullet.c bullet.h graph.h initgraph.h plotproc.h wop.h wpm.h

libuggraph2_la_SOURCES = $(SOURCES)
libuggraph2_la_CPPFLAGS = $(UG_2DFLAGS)

libuggraph3_la_SOURCES = $(SOURCES)
libuggraph3_la_CPPFLAGS = $(UG_3DFLAGS)
