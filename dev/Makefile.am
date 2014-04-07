# $Id$

# !!! test and include them:
#
# mif (MacOS)

SUBDIRS = ppm meta ps . sif xif rif $(EXTRADIR)

if MACOS
  EXTRADIR = mif
else
  EXTRADIR =
endif

noinst_LTLIBRARIES = libdevices.la

# normal flags and the position of the installed plugins
AM_CPPFLAGS = $(UG_CPPFLAGS) $(UG_PLUGINDIR)

SOURCES = ugdevices.c sockcomm.c initdev.h sockcomm.h ugdevices.h

devincludedir = $(pkgincludedir)
devinclude_HEADERS = ugdevices.h

libdevices_la_SOURCES = $(SOURCES)
libdevices_la_LIBADD = ppm/libdevice_ppm.la \
                       meta/libdevice_meta.la \
                       ps/libdevice_ps.la

# for system with dynamic drivers:
#
# SOURCES += loader.c active_device.h devmodule.h passive_device.h
# libdevices_la_LIBADD = $(LIBLTDL)

## small tester for the device interface
#
#noinst_PROGRAMS = tester
#
#tester_SOURCES = tester.c
#tester_LDFLAGS = -dlopen self
# -dlopen sif/device_s.la -dlopen xif/device_x.la
#tester_LDADD = libdevices.la ../low/liblow.la

# TESTS = tester
