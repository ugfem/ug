############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################

# include configuration and all makefile macro definitions
include $(UGROOT)/ug.conf

#  The name of the libraries to be created/linked to
SUPERLULIB   = $(UG_LIB)
BLASLIB      = $(UG_LIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = $(ARCH_AR)
ARCHFLAGS    = $(ARCH_ARFLAGS)
RANLIB       = echo

CC           = $(ARCH_CC)
CFLAGS       = $(UG_CFLAGS)
CDEFS        = -DAdd_










