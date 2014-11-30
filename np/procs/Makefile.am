SOURCES =\
  nls.cc \
  ls.cc \
  assemble.cc \
  basics.cc \
  blocking.cc \
  iter.cc \
  transfer.cc \
  error.cc \
  newton.cc \
  ts.cc \
  bdf.cc \
  ew.cc \
  ewn.cc \
  amgtransfer.cc \
  freebnd.cc \
  db.cc \
  fas.cc \
  nliter.cc \
  project.cc \
  order.cc \
  tstep.cc \
  enewton.cc \
  pstep.cc \
  reinit.cc \
  els.cc \
  eiter.cc \
  iter_2.cc \
  amgtransfer.h \
  assemble.h \
  basics.h \
  bdf.h \
  blocking.h \
  db.h \
  eiter.h \
  els.h \
  enewton.h \
  error.h \
  ew.h \
  fas.h \
  freebnd.h \
  iter.h \
  ls.h \
  newton.h \
  nliter.h \
  nls.h \
  order.h \
  project.h \
  pstep.h \
  reinit.h \
  transfer.h \
  ts.h \
  tstep.h

npincludedir = $(pkgincludedir)
npinclude_HEADERS = assemble.h transfer.h ts.h nls.h ls.h

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
