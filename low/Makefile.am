# $Id$

noinst_LTLIBRARIES = liblow.la libtree2.la libtree3.la

# dimension-independent functions
liblow_la_SOURCES = ugenv.c heaps.c fifo.c misc.c \
                  defaults.c initlow.c fileopen.c \
                  ugstruct.c debug.c bio.c pfile.c \
                  scan.c ugtimer.c \
architecture.h bio.h boxtree.h  debug.h  defaults.h  dimension.h  fifo.h \
fileopen.h general.h heaps.h initlow.h misc.h pfile.h \
scan.h smallvecs.h tree.h ugenv.h ugfile.h ugmemory.h ugstrings.h \
ugstruct.h ugtimer.h ugtypes.h ugtime.h namespace.h bbtree.c bbtree.h
# the header in SOURCES tell automake to distribute them in the tarball

liblow_la_CPPFLAGS = $(UG_CPPFLAGS)

# normal headers
include_HEADERS = heaps.h ugenv.h ugstruct.h misc.h debug.h ugtypes.h \
                  general.h dimension.h smallvecs.h \
                  tree.h boxtree.h ugtimer.h scan.h defaults.h fileopen.h \
                  namespace.h

# tree and boxtree are dimension dependent
TREESOURCES = tree.c boxtree.c

libtree2_la_SOURCES = $(TREESOURCES)
libtree2_la_CPPFLAGS = $(UG_2DFLAGS)

libtree3_la_SOURCES = $(TREESOURCES)
libtree3_la_CPPFLAGS = $(UG_3DFLAGS)
