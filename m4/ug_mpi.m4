# $Id$

# wrapper for the autoconf-archive check. Note: compiling MPI-stuff sucks!

# Explanation:
# ============
#

# compilation of MPI-programs is normally done by a
# mpicc/mpiCC-wrapper that adds all options needed. Thus, it may seem
# possible to just replace the compiler call by the wrapper and
# everything works. Unfortunately that's not the case: automake and
# libtool both show strange behaviour.
#
# In detail: replacing the compiler globally via ./configure CXX=mpiCC
# should work (at least I've found reports claiming this) but that is
# not what we want: mainly, it just adds a level of possible errors
# (mpiCC from MPICH does _nothing_ if "mpicc -c dummy.cc" is called!)
# and might introduce nice library-clashes.
#
# The next approach would be to include
#       if MPI
#         CXX = $(MPICXX)
#       endif
# in the Makefile.am where MPI is needed. First, this will change
# compilations of all binaries in this directory and secondly the
# dependency-tracking seems to break: the first compilation worked but
# the second failed with the compiler complaining about mismatching
# flags... There is no 'program_CXX = ...' in automake but even if
# there were it would break as well
#
# Thus, the best solution is to extract the flags needed for
# compilation and linking. Unfortunately, the parameters and behaviour
# of mpicc is not at all consistent over different
# implementations. For MPICH the parameters -compile_info and
# -link_info exist (albeit not being documented in the manpage, only
# in -help), for LAM dummy-calls of compilation and linking together
# with a -showme parameter (which is called -show in MPICH...) have to
# be used. Obviously, we have to identify the type of package... this
# is done via mpiCC-calls for now, I wouldn't be surprised if ths
# breaks often. Bad luck. Blame the MPI folks for this mess. And blame
# them a lot. [Thimo 26.8.2004]

# TODO:
#
# - add --disable-mpi

AC_DEFUN([UG_MPI],[
  AC_PREREQ(2.50) dnl for AC_LANG_CASE

  AS_IF([test -n "$MPICXX" -a -z "$MPICC"],[
    MPICC=$MPICXX
    AC_MSG_WARN(["You specified the MPICXX variable, but not MPICC. UG extracts the MPI flags using the C interface!])
  ])
  
  AC_LANG_PUSH([C])

  # enable/disable parallel features
  AC_ARG_ENABLE(parallel,
    AS_HELP_STRING([--enable-parallel],
      [Enable the parallel features of UG. If enabled
       configure will try to determine your MPI automatically. You can
       overwrite this setting by specifying the MPICC variable]))
  AC_SUBST(ENABLE_PARALLEL, "$enable_parallel")

  ## do nothing if --disable-parallel is used
  AS_IF([test "x$enable_parallel" = "xyes"],[
    ACX_MPI([
      MPICOMP="$MPICC"

      MPI_CONFIG()
      MPI_CPPFLAGS="$DUNEMPICPPFLAGS $MPI_NOCXXFLAGS -DENABLE_MPI=1"

      with_mpi="yes ($dune_MPI_VERSION)"
    ],[
      # ACX_MPI didn't find anything
      with_mpi="no"
    ])])

  # if an MPI implementation was found..
  AS_IF([test "x$with_mpi" != "xno"],[
    ### do a sanity check: can we compile and link a trivial MPI program?
    AC_MSG_CHECKING([whether compiling with $dune_MPI_VERSION works])

    # store old values
    ac_save_LIBS="$LIBS"
    ac_save_LDFLAGS="$LDFLAGS"
    ac_save_CPPFLAGS="$CPPFLAGS"
    
    # looks weird but as the -l... are contained in the MPI_LDFLAGS these
    # parameters have to be last on the commandline: with LIBS this is true
    LIBS="$DUNEMPILIBS $LIBS"
    LDFLAGS="$LDFLAGS $DUNEMPILDFLAGS"
    CPPFLAGS="$CPPFLAGS $DUNEMPICPPFLAGS"

    # try to create MPI program
    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE(
      AC_LANG_SOURCE(
        [ #include <mpi.h>
          int main (int argc, char** argv) { 
          MPI_Init(&argc, &argv); 
          MPI_Finalize(); }]),
        [ AC_MSG_RESULT([yes]) ],
        [ AC_MSG_RESULT([no])
          AC_MSG_ERROR([could not compile MPI testprogram!
          See config.log for details])
          with_mpi=no]
    )
    AC_LANG_POP

    # Check for MPI-2 Standard
    # We have to provide a dummy lib here as we do not know what the name
    # of the mpi is. -lm should be save.
    AC_CHECK_LIB(m,[MPI_Finalized], [AC_DEFINE(MPI_2, 1, [Define to 1 MPI supports MPI-2])])

    # restore variables
    LIBS="$ac_save_LIBS"
    CPPFLAGS="$ac_save_CPPFLAGS"
  ])
    
  # set flags
  AS_IF([test "x$with_mpi" != "xno"],[
    MPI_LDFLAGS="$DUNEMPILDFLAGS $DUNEMPILIBS"
    AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])
    AC_SUBST(MPI_CPPFLAGS, $DUNEMPICPPFLAGS)
    AC_SUBST(MPI_LDFLAGS, $MPI_LDFLAGS)
  ],[
    AC_SUBST(MPI_CPPFLAGS, "")
    AC_SUBST(MPI_LDFLAGS, "")
  ])
  AC_SUBST(MPI_VERSION, $dune_MPI_VERSION)
  AM_CONDITIONAL(MPI, [test "x$with_mpi" != "xno"])

  AC_LANG_POP
])
