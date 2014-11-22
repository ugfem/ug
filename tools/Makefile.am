AM_CFLAGS = $(UG_CPPFLAGS)

# for plugins later:
#
# !!! dynamic options should be global or so... or in the lib?
# AM_LDFLAGS = -dlopen self

# xugv needs X11 and XAW, the latter includes the first
if XAW
 XUGV = xugv
else
 XUGV =
endif

TOOLS = ugnetgen ugstl $(XUGV) ugmesh m2ps

if BUILD2D
  TOOLS_2D = ugmgs2 ugpfm2
else
  TOOLS_2D =
endif

if BUILD3D
  TOOLS_3D = ugmgs3 ugpfm3
else
  TOOLS_3D =
endif

bin_PROGRAMS = $(TOOLS) $(TOOLS_2D) $(TOOLS_3D)

xugv_LDFLAGS = $(UG_XLIBS)

# !!! std-domain / device SIF ok?
UG_2DLIBS = $(top_builddir)/lib/libugS2.la \
  $(top_builddir)/dev/sif/libdevS.la
UG_3DLIBS = $(top_builddir)/lib/libugS3.la \
  $(top_builddir)/dev/sif/libdevS.la

ugnetgen_SOURCES = ugnetgen.cc
ugnetgen_CPPFLAGS = $(UG_2DFLAGS)
ugnetgen_LDFLAGS = $(UG_LDFLAGS)
ugnetgen_LDADD = $(UG_2DLIBS)

ugstl_SOURCES = ugstl.cc
ugstl_CPPFLAGS = $(UG_2DFLAGS)
ugstl_LDFLAGS = $(UG_LDFLAGS)
ugstl_LDADD = $(UG_2DLIBS)

xugv_SOURCES = xugv.cc
xugv_CPPFLAGS = $(UG_2DFLAGS)
xugv_LDFLAGS = $(UG_LDFLAGS)
xugv_LDADD = $(UG_2DLIBS)

ugmesh_SOURCES = ugmesh.cc
ugmesh_CPPFLAGS = $(UG_2DFLAGS)
ugmesh_LDFLAGS = $(UG_LDFLAGS)
ugmesh_LDADD = $(UG_2DLIBS)

m2ps_SOURCES = m2ps.cc
m2ps_CPPFLAGS = $(UG_2DFLAGS)
m2ps_LDFLAGS = $(UG_LDFLAGS)
m2ps_LDADD = $(UG_2DLIBS)

ugpfm2_SOURCES = ugpfm.cc
ugpfm2_CPPFLAGS = $(UG_2DFLAGS)
ugpfm2_LDFLAGS = $(UG_LDFLAGS)
ugpfm2_LDADD = $(UG_2DLIBS)

ugpfm3_SOURCES = ugpfm.cc
ugpfm3_CPPFLAGS = $(UG_3DFLAGS)
ugpfm3_LDFLAGS = $(UG_LDFLAGS)
ugpfm3_LDADD = $(UG_3DLIBS)

ugmgs2_SOURCES = ugmgs.cc
ugmgs2_CPPFLAGS = $(UG_2DFLAGS)
ugmgs2_LDFLAGS = $(UG_LDFLAGS)
ugmgs2_LDADD = $(UG_2DLIBS)

ugmgs3_SOURCES = ugmgs.cc
ugmgs3_CPPFLAGS = $(UG_3DFLAGS)
ugmgs3_LDFLAGS = $(UG_LDFLAGS)
ugmgs3_LDADD = $(UG_3DLIBS)
