# $Id$

AM_CFLAGS = $(UG_CPPFLAGS)

# for plugins later:
#
# !!! dynamic options should be global or so... or in the lib?
# AM_LDFLAGS = -dlopen self

if X11
 XUGV = xugv
else
 XUGV =
endif

TOOLS = ugnetgen ugstl $(XUGV) ugmesh m2ps

if BUILD2D
  TOOLS_2D = ugmgs2 ugpfm2 ugmarc2
else
  TOOLS_2D =
endif

if BUILD3D
  TOOLS_3D = ugmgs3 ugpfm3 ugmarc3
else
  TOOLS_3D =
endif

bin_PROGRAMS = $(TOOLS) $(TOOLS_2D) $(TOOLS_3D)

xugv_LDFLAGS = $(UG_XLIBS)

# !!! std-domain / device SIF ok?
UG_2DLIBS = $(top_builddir)/lib/libug2.la $(top_builddir)/dom/std/libdomS2.la \
  $(top_builddir)/dev/sif/libdevS.la
UG_3DLIBS = $(top_builddir)/lib/libug3.la $(top_builddir)/dom/std/libdomS3.la \
  $(top_builddir)/dev/sif/libdevS.la

ugpfm2_SOURCES = ugpfm.c
ugpfm2_CPPFLAGS = $(UG_2DFLAGS)
ugpfm2_LDFLAGS = $(UG_LDFLAGS)
ugpfm2_LDADD = $(UG_2DLIBS)

ugpfm3_SOURCES = ugpfm.c
ugpfm3_CFLAGS = $(UG_3DFLAGS)
ugpfm3_LDFLAGS = $(UG_3DLIBS)

ugmgs2_SOURCES = ugmgs.c
ugmgs2_CPPFLAGS = $(UG_2DFLAGS)
ugmgs2_LDFLAGS = $(UG_LDFLAGS)
ugmgs2_LDADD = $(UG_2DLIBS)

ugmgs3_SOURCES = ugmgs.c
ugmgs3_CPPFLAGS = $(UG_2DFLAGS)
ugmgs3_LDFLAGS = $(UG_LDFLAGS)
ugmgs3_LDADD = $(UG_2DLIBS)

ugmarc2_SOURCES = ugmarc.c
ugmarc2_CPPFLAGS = $(UG_2DFLAGS)
ugmarc2_LDFLAGS = $(UG_LDFLAGS)
ugmarc2_LDADD = $(UG_2DLIBS)

ugmarc3_SOURCES = ugmarc.c
ugmarc3_CPPFLAGS = $(UG_3DFLAGS)
ugmarc3_LDFLAGS = $(UG_LDFLAGS)
ugmarc3_LDADD = $(UG_3DLIBS)
