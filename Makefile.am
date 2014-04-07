# $Id$

# do everything in the subdirs (order matters!)
SUBDIRS = parallel low dev gm dom graphics np ui lib tools
#
# libltdl after low will be needed for plugins later

ugincludedir = $(pkgincludedir)
uginclude_HEADERS = initug.h

# install pkg-config information
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libug.pc

#install cmake package configuration
cmakeconfigdir = $(libdir)/cmake/ug
cmakeconfig_DATA = ug-config.cmake ug-config-version.cmake

#
## automake-options
#

# don't follow the full GNU-standard
AUTOMAKE_OPTIONS = foreign 1.7

ACLOCAL_AMFLAGS = -I m4

# possibly use the same compiler for distcheck
# DISTCHECK_CONFIGURE_FLAGS = CXX="$(CXX)" CC="$(CC)"
