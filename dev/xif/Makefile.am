# $Id$

# commands.h will gm.h und der braucht Dimension
AM_CPPFLAGS = $(UG_2DFLAGS) -I$(top_srcdir)/low

lib_LTLIBRARIES = libdevX.la

libdevX_la_SOURCES = xmain.c xgraph.c xshell.c \
  shades.h xgraph.h xmain.h xshell.h
libdevX_la_LIBADD = ../libdevices.la

# put icons into tarball as well...
EXTRA_DIST = tool0 tool1 tool2 tool3 tool4 tool5 tool6 view-icon shell-icon

#plugin_LTLIBRARIES = device_x.la
#
#device_x_la_SOURCES = xmain.c xgraph.c xshell.c \
#  shades.h xgraph.h xmain.h xshell.h
#device_x_la_LDFLAGS = -module -no-undefined -avoid-version
#device_x_la_LIBADD = $(UG_XLIBS)

# !!! gui-Toolbox?
