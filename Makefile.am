# $Id$

# do everything in the subdirs (order matters!)
SUBDIRS = arch low dev dom gm graphics np ui parallel lib tools
#
# libltdl after low will be needed for plugins later

include_HEADERS = initug.h

# install pkg-config information
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libug.pc

# !!! crazy trick to get include-dir compatibility
install-data-hook:
	find $(top_srcdir) -name \*.h -exec cp {} $(includedir)/ \;
	rm -f $(includedir)/config.h
	touch $(includedir)/config.h

#
## automake-options
#

# don't follow the full GNU-standard
AUTOMAKE_OPTIONS = foreign 1.7
