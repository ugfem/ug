# $Id$

# do everything in the subdirs (order matters!)
SUBDIRS = arch low dev dom gm graphics np ui parallel lib tools
#
# libltdl after low will be needed for plugins later

include_HEADERS = initug.h

# install pkg-config information
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libug.pc


#
## automake-options
#

# don't follow the full GNU-standard
AUTOMAKE_OPTIONS = foreign 1.7
