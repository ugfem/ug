// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>
#include        <string.h>
#include        "../main/defs.h"
#include        "../main/params.h"

/* Robust routine to read a non-negative integer */
int input_int()
{
  char line[LINE_LENGTH];       /* space to read input line */
  int done;                     /* flag for end of integer */
  int val;                      /* value returned */
  int i;                        /* loop counter */
  int isdigit();

  i = 0;
  done = FALSE;
  while (!done) {
    line[i] = getchar();
    if (isdigit(line[i])) i++;
    else if (i != 0) done = TRUE;
  }

  sscanf(line, "%d", &val);
  return(val);
}
