# $Id$

# do everything in the subdirs (order matters!)
SUBDIRS = arch parallel low dev gm dom graphics np ui lib tools
#
# libltdl after low will be needed for plugins later

include_HEADERS = initug.h

# install pkg-config information
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libug.pc

# !!! crazy trick to get include-dir compatibility
# Crazy indeed.  Let's at least not install f2c.h (because it creates
# problems) and std_internal.h (because it is newer than the build system).
install-data-hook:
	find $(top_srcdir) -name \*.h -not -name "f2c.h" -not -name "std_internal.h" -exec $(INSTALL_DATA) {} $(DESTDIR)/$(includedir)/ \;
	rm -f $(DESTDIR)/$(includedir)/config.h

# distcheck checks uninstall as well: revert sick trick above
#uninstall-hook:
#	for HEADER in `find $(top_srcdir) -name \*.h` ; do \
#	  BASE=`basename $$HEADER` ; \
#	  rm -f $(includedir)/$$BASE ; \
#	done

#
## automake-options
#

# don't follow the full GNU-standard
AUTOMAKE_OPTIONS = foreign 1.7

# possibly use the same compiler for distcheck
# DISTCHECK_CONFIGURE_FLAGS = CXX="$(CXX)" CC="$(CC)"
