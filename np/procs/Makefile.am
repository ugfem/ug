# $Id$

SOURCES = nls.c ls.c assemble.c basics.c iter.c transfer.c error.c newton.c ts.c bdf.c ew.c ewn.c amgtransfer.c freebnd.c db.c fas.c nliter.c project.c order.c tstep.c enewton.c pstep.c reinit.c els.c eiter.c iter_2.c amgtransfer.h assemble.h basics.h bdf.h db.h eiter.h els.h enewton.h error.h ew.h fas.h freebnd.h iter.h ls.h newton.h nliter.h nls.h order.h project.h pstep.h reinit.h transfer.h ts.h tstep.h

include_HEADERS = assemble.h transfer.h error.h ts.h nls.h ls.h

if BUILD2D
  LIB_2D = libprocs2.la
else
  LIB_2D =
endif

if BUILD3D
  LIB_3D = libprocs3.la
else
  LIB_3D =
endif

noinst_LTLIBRARIES = $(LIB_2D) $(LIB_3D)

libprocs2_la_SOURCES = $(SOURCES)
libprocs2_la_CPPFLAGS = $(UG_2DFLAGS)

libprocs3_la_SOURCES = $(SOURCES)
libprocs3_la_CPPFLAGS = $(UG_3DFLAGS)
