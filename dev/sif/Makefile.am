# $Id$

AM_CPPFLAGS = $(UG_CPPFLAGS) -I$(top_srcdir)/low

lib_LTLIBRARIES = libdevS.la
libdevS_la_SOURCES = sif.c
libdevS_la_LDFLAGS = -release $(PACKAGE_VERSION)
libdevS_la_LIBADD = ../libdevices.la

# for the plugin-stuff later:
#
#plugin_LTLIBRARIES = device_s.la
#
#device_s_la_SOURCES = sif.c
#device_s_la_LDFLAGS = -module -no-undefined -avoid-version

