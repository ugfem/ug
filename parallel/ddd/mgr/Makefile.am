SUBDIRS =

if PARALLEL
  PLIB = libmgr2.la libmgr3.la
else
  PLIB =
endif

noinst_LTLIBRARIES = $(PLIB)

EXTRA_DIST = handler.ct

SOURCES = typemgr.cc objmgr.cc cplmgr.cc prio.cc attr.cc

libmgr2_la_SOURCES = $(SOURCES)
libmgr2_la_CPPFLAGS = $(UG_2DFLAGS) $(DDDFLAGS)

libmgr3_la_SOURCES = $(SOURCES)
libmgr3_la_CPPFLAGS = $(UG_3DFLAGS) $(DDDFLAGS)
