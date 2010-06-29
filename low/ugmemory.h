// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* $Id$ */

/*

   portable memory-allocation, makes the following functions usable:

   malloc
   free
   alloca
   realloc

 */

#ifndef UG_UGMEMORY_H
#define UG_UGMEMORY_H

#ifndef UGLIB
#error Internal UG-lib header, must not be used in applications!
#endif

/*  try ISO header first */
#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif

/* !!! malloc.h seems to define ??? */

/* !!! missing:

   C_ALLOCA
   HAVE_ALLOCA
   HAVE_MALLOC
   HAVE_REALLOC

 */

#ifndef HAVE_ALLOCA
#error need replacement-alloca.. :(
#endif

/* -- treat alloca like the autoconf-docs suggest -- */

/* AIX requires this to be the first thing in the file.  */
#ifndef __GNUC__
# if HAVE_ALLOCA_H
#  include <alloca.h>
# else
#  ifdef _AIX
 #pragma alloca
#  else
#   ifndef alloca /* predefined by HP cc +Olibcalls */
char *alloca ();
#   endif
#  endif
# endif
#endif

#endif
