#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(UG, 3.9.1, UG-REPORT-ADDRESS)
AM_INIT_AUTOMAKE(UG, 3.9.1, UG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([initug.c])
AM_CONFIG_HEADER([config.h])

AM_MAINTAINER_MODE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_LN_S
# for dom/lgm/ngin
AM_PROG_LEX
AC_PROG_YACC

# libtool (lines stolen from libtool-docs)
dnl Enable building of the convenience library
dnl and set LIBLTDL accordingly
AC_LIBLTDL_CONVENIENCE
dnl Substitute LTDLINCL and LIBLTDL in the Makefiles
AC_SUBST(LTDLINCL)
AC_SUBST(LIBLTDL)
# was f�r libs wollen wir:
#AC_ENABLE_SHARED
#AC_DISABLE_STATIC
dnl Check for dlopen support
AC_LIBTOOL_DLOPEN
dnl Configure libtool
AC_PROG_LIBTOOL
dnl Configure libltdl
AC_CONFIG_SUBDIRS(libltdl)

# get Fortran-compiler-specs (for ddd-FORTRAN-bindings)
AC_PROG_F77
if test x$ac_cv_prog_f77_g != xno ; then
   AC_F77_MAIN
   AC_F77_WRAPPERS
else
   # tell about missing FORTRAN
   AC_DEFINE(NOF77, 1, [define if no FORTRAN-compiler is available])
fi

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lfamg':
AC_CHECK_LIB([famg], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])

# !!!
AM_CONDITIONAL(NETGEN, false)

# Checks for header files.
AC_PATH_XTRA
AC_FUNC_ALLOCA
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS([arpa/inet.h fcntl.h float.h limits.h malloc.h memory.h netdb.h netinet/in.h stddef.h stdlib.h string.h strings.h sys/param.h sys/socket.h unistd.h values.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])
AC_C_BIGENDIAN

# Checks for library functions.
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([bzero floor gethostbyaddr gethostbyname gethrtime memmove memset mkdir pow socket sqrt strchr strdup strrchr strstr strtol])

#
## parallel support
#

AC_ARG_ENABLE([parallel],AC_HELP_STRING([--enable-parallel],[enable parallel computing (default is NO)]),[parallel=$enableval],[parallel=no])

if test x$parallel = xyes ; then
  # set define
  UG_EXTRAFLAGS='-DModelP '
fi

# debugging

AC_ARG_ENABLE([debug],AC_HELP_STRING([--enable-debug],[enable UG debugging (default is NO)]),[debug=$enableval],[debug=no])

if test x$debug = xyes ; then
  # set define
  UG_EXTRAFLAGS='-DDebug '
fi

# the DUNE-test for MPI-package and parameters
DUNE_MPI

# set results as automake-conditionals
AM_CONDITIONAL(PARALLEL, test x$parallel = xyes)

#
## build options
#

# !!! reduce include pathes!!!
UG_CPPFLAGS="-I\$(top_srcdir) -I\$(top_srcdir)/arch -I\$(top_srcdir)/low -I\$(top_srcdir)/gm -I\$(top_srcdir)/dev -I\$(top_srcdir)/dom -I\$(top_srcdir)/np -I\$(top_srcdir)/ui -I\$(top_srcdir)/graphics/uggraph -I\$(top_srcdir)/np/algebra -I\$(top_srcdir)/np/udm -I\$(top_srcdir)/np/procs -I\$(top_srcdir)/parallel/ppif -I\$(top_srcdir)/parallel/ddd/include -I\$(top_srcdir)/parallel/dddif -DUGLIB $UG_EXTRAFLAGS"
AC_SUBST(UG_CPPFLAGS)
# export EXTRAFLAGS for libug.pc, CPPFLAGS would be quite wrong :)
AC_SUBST(UG_EXTRAFLAGS)
AC_SUBST(UG_2DFLAGS, '$(UG_CPPFLAGS) -D_2')
AC_SUBST(UG_3DFLAGS, '$(UG_CPPFLAGS) -D_3')

AC_ARG_ENABLE([2d],AC_HELP_STRING([--enable-2d],[build library for 2d (default is YES)]),[build2d=$enableval],[build2d=yes])
AM_CONDITIONAL(BUILD2D, test x$build2d = xyes)

AC_ARG_ENABLE([3d],AC_HELP_STRING([--enable-3d],[build library for 3d (default is YES)]),[build3d=$enableval],[build3d=yes])
AM_CONDITIONAL(BUILD3D, test x$build3d = xyes)

# lib-stuff
AC_SUBST(UG_LDFLAGS, '-L$(top_srcdir)/lib')

# where should the plugins go to?
plugindir=${libdir}/ug/
AC_SUBST(plugindir)
AC_SUBST(UG_PLUGINDIR, '-DPLUGINDIR=\"$(plugindir)\"')

# already expand variable (no obvious drawback and we can insert the
# value into libug.pc
UG_XLIBS="$X_PRE_LIBS $X_LIBS -lX11 $X_EXTRA_LIBS -lXt -lXaw"
AC_SUBST(UG_XLIBS)

# move architecture-parts into header to still allow crude tricks :(
AC_DEFINE_UNQUOTED(ARCH_CPU, "$build_cpu", [CPU type])
AC_DEFINE_UNQUOTED(ARCH_VENDOR, "$build_vendor", [OS Vendor])
AC_DEFINE_UNQUOTED(ARCH_OS, "$build_os", [OS type])
# !!! somehow backward compatible Symbol... (should be replaced by the
# defines above)
AC_DEFINE_UNQUOTED(ARCHNAME, "$build", [obsolete architecture string])

AC_DEFINE(DYNAMIC_MEMORY_ALLOCMODEL, 1, [UG memory allocation model])

# a special define to tell compiler.h to use the autoconf-stuff so
# that we can keep both build systems compatible for a while
AC_DEFINE(AUTOTOOLS_BUILD, 1, [switch headers to autoconf-support])

# tell automake what to build
AM_CONDITIONAL(MACOS, test x$build_vendor = xapple)

AC_CONFIG_FILES([Makefile
                 arch/Makefile
                 dev/Makefile
		 dev/meta/Makefile
		 dev/mif/Makefile
		 dev/ppm/Makefile
		 dev/ps/Makefile
		 dev/rif/Makefile
		 dev/sif/Makefile
		 dev/xif/Makefile
		 dom/Makefile
		 dom/std/Makefile
		 dom/lgm/Makefile
		 dom/lgm/ngin/Makefile
		 dom/lgm/ngin2d/Makefile
		 dom/gen/Makefile
                 gm/Makefile
		 gm/gg2/Makefile
		 gm/gg3/Makefile
		 graphics/Makefile
		 graphics/uggraph/Makefile
		 graphics/grape/Makefile
		 lib/Makefile
		 low/Makefile
		 np/Makefile
		 np/algebra/Makefile
		 np/amglib/Makefile
		 np/field/Makefile
		 np/procs/Makefile
		 np/slu/Makefile
		 np/slu/cblas/Makefile
		 np/slu/src/Makefile
		 np/udm/Makefile
		 tools/Makefile
                 ui/Makefile
                 parallel/Makefile
                 parallel/util/Makefile
                 libug.pc
])
AC_OUTPUT

dnl Still missing:
dnl
dnl                  dev/mif/Makefile
dnl                  gm/gg3/netgen/general/makefile
dnl                  gm/gg3/netgen/geom/makefile
dnl                  gm/gg3/netgen/linalg/makefile
dnl                  gm/gg3/netgen/meshing/makefile
dnl                  gm/gg3/netgen/opti/makefile
dnl                  parallel/chaco/assign/Makefile
dnl                  parallel/chaco/bpmatch/Makefile
dnl                  parallel/chaco/coarsen/Makefile
dnl                  parallel/chaco/connect/Makefile
dnl                  parallel/chaco/eigen/Makefile
dnl                  parallel/chaco/graph/Makefile
dnl                  parallel/chaco/inertial/Makefile
dnl                  parallel/chaco/input/Makefile
dnl                  parallel/chaco/klspiff/Makefile
dnl                  parallel/chaco/main/Makefile
dnl                  parallel/chaco/misc/Makefile
dnl                  parallel/chaco/optimize/Makefile
dnl                  parallel/chaco/symmlq/Makefile
dnl                  parallel/chaco/util/Makefile