// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
   This nothingness makes several, say n, broadcasts into
   one. Should be n times faster for broadcasting n small
   variables not consecutive in memory than doing n single
   broadcasts.
 */

/*
   $Header$
 */

#include "config.h"
#include <stddef.h>
#include <stdarg.h>
#include <string.h>

#include "ppif.h"

USING_PPIF_NAMESPACE

#define XBCMAX 32768

static char buffer[XBCMAX];

void XBroadcast(int n, void *p, size_t s, ...)
{
  char *b;
  void *pp;
  int i;
  size_t ss;
  va_list a, aa;

  b = buffer;
  memcpy(b, p, s); b += s;

  pp = p; ss = s;
  va_start(a, s);

#if defined va_copy
  va_copy(aa, a);            /* C99 */
#elif defined __va_copy
  __va_copy (aa, a);         /* draft C99 */
#else
  aa = a;                    /* hope */
#endif

  for (i = 1; i < n; i++) {
    p = va_arg(a, void*);
    s = va_arg(a, size_t);
    memcpy(b, p, s); b += s;
  }

  Broadcast(buffer, b-buffer);

  b = buffer;
  memcpy(pp, b, ss); b += ss;

  for (i = 1; i < n; i++) {
    p = va_arg(aa, void*);
    s = va_arg(aa, size_t);
    memcpy(p, b, s); b += s;
  }

  va_end(a); va_end(aa);
}
