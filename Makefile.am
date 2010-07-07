# $Id$

# do everything in the subdirs (order matters!)
SUBDIRS = arch parallel low dev gm dom graphics np ui lib tools
#
# libltdl after low will be needed for plugins later

ugincludedir = $(pkgincludedir)
uginclude_HEADERS = initug.h

# install pkg-config information
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libug.pc

# !!! crazy trick to get include-dir compatibility
# Crazy indeed.  Let's at least not install f2c.h (because it creates
# problems) and std_internal.h (because we definitely know that it is an internal header).
install-data-hook:
	test -z "$(DESTDIR)$(pkgincludedir)" || $(MKDIR_P) "$(DESTDIR)$(pkgincludedir)"
	$(INSTALL_DATA) $(top_srcdir)/arch/compiler.h "$(DESTDIR)$(pkgincludedir)"
	find $(top_srcdir)/parallel -name \*.h -not -name "f2c.h" -exec $(INSTALL_DATA) {} "$(DESTDIR)$(pkgincludedir)" \;
	find $(top_srcdir)/low -name \*.h -exec $(INSTALL_DATA) {} "$(DESTDIR)$(pkgincludedir)" \;
	find $(top_srcdir)/dev -name \*.h -exec $(INSTALL_DATA) {} "$(DESTDIR)$(pkgincludedir)" \;
	find $(top_srcdir)/gm -name \*.h -exec $(INSTALL_DATA) {} "$(DESTDIR)$(pkgincludedir)" \;
	find $(top_srcdir)/dom -name \*.h -not -name "std_internal.h" -exec $(INSTALL_DATA) {} "$(DESTDIR)$(pkgincludedir)" \;
	find $(top_srcdir)/graphics -name \*.h -exec $(INSTALL_DATA) {} "$(DESTDIR)$(pkgincludedir)" \;
	find $(top_srcdir)/np -name \*.h -not -name "error.h" -exec $(INSTALL_DATA) {} "$(DESTDIR)$(pkgincludedir)" \;

# distcheck checks uninstall as well: revert sick trick above
uninstall-hook:
	rm -rf $(DESTDIR)/$(pkgincludedir)

#
## automake-options
#

# don't follow the full GNU-standard
AUTOMAKE_OPTIONS = foreign 1.7

# possibly use the same compiler for distcheck
# DISTCHECK_CONFIGURE_FLAGS = CXX="$(CXX)" CC="$(CC)"
