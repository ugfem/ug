noinst_LTLIBRARIES = liblow.la liblow2.la liblow3.la

# dimension-independent functions
liblow_la_SOURCES = ugenv.cc heaps.cc fifo.cc misc.cc \
                  defaults.cc initlow.cc fileopen.cc \
                  ugstruct.cc debug.cc bio.cc \
                  scan.cc ugtimer.cc tree.cc \
architecture.h bio.h boxtree.h  debug.h  defaults.h  dimension.h  fifo.h \
fileopen.h general.h heaps.h initlow.h misc.h pfile.h \
scan.h smallvecs.h tree.h ugenv.h ugfile.h ugmemory.h ugstrings.h \
ugstruct.h ugtimer.h ugtypes.h ugtime.h namespace.h bbtree.cc bbtree.h
# the header in SOURCES tell automake to distribute them in the tarball

liblow_la_CPPFLAGS = $(UG_CPPFLAGS)

# normal headers
lowincludedir = $(pkgincludedir)
lowinclude_HEADERS = heaps.h ugenv.h ugstruct.h misc.h debug.h ugtypes.h \
                  general.h dimension.h smallvecs.h \
                  tree.h boxtree.h ugtimer.h scan.h defaults.h fileopen.h \
                  namespace.h pfile.h

# special library for the parts of low that need a dimension
DIMDEPSOURCES = boxtree.cc pfile.cc

liblow2_la_SOURCES = $(DIMDEPSOURCES)
liblow2_la_CPPFLAGS = $(UG_2DFLAGS)
liblow2_la_LIBADD = liblow.la

liblow3_la_SOURCES = $(DIMDEPSOURCES)
liblow3_la_CPPFLAGS = $(UG_3DFLAGS)
liblow3_la_LIBADD = liblow.la
