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

#ifndef __XBC__
#define __XBC__

void XBroadcast(int n, void *p, int s, ...);

#endif
