// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* $Id$ */

/*

  portable string-functions

*/

#ifndef UG_UGSTRINGS_H
#define UG_UGSTRINGS_H

#ifndef UGLIB
#error Internal UG-lib header, must not be used in applications!
#endif

#include "config.h"

/*
  Try ISO C header <string.h>

  http://www.opengroup.org/onlinepubs/7908799/xsh/string.h.html

  if it's not availabe try the <strings.h> which doesn't seem to have
  any standard associated
*/

/* stolen from autoconf-docs */
#if STDC_HEADERS
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  ifdef HAVE_STRINGS_H
#   include <strings.h>
#  endif
# endif
#else
# if !HAVE_STRCHR
#  define strchr index
#  define strrchr rindex
# endif
char *strchr (), *strrchr ();
# if !HAVE_MEMCPY
#  define memcpy(d, s, n) bcopy ((s), (d), (n))
#  define memmove(d, s, n) bcopy ((s), (d), (n))
# endif
#endif

/* --- replacement functions for those that may be missing --- */

/*

  we include the implementation directly for two reasons:

  1. even when no 'static' is used the function can never appear in
     the libs symbols (no clashes with other libs which may possibly
     include replacements themselves
     
  2. if the compiler supports it   
     
 */

/* !!! check for:

  HAVE_MEMSET
  HAVE_MEMMOVE
  HAVE_MEMORY_H
  HAVE_BZERO

  HAVE_STRRCHR
  HAVE_STRSTR
  HAVE_STRTOL

*/

#ifndef HAVE_STRDUP
/* need malloc */
#include "ugmemory.h"
static inline char *strdup(const char *s) {    
  return strcpy(malloc((strlen(s)+1)*sizeof(char), s);
}
#endif

#endif
