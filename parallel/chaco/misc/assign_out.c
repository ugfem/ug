// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>

void assign_out(nvtxs, sets, outname)
int nvtxs;                      /* number of vertices to output */
short *sets;                    /* values to be printed */
char *outname;                  /* name of output file */
{
  FILE *fout;                   /* output file */
  int i;                        /* loop counter */

  fout = fopen(outname, "w");
  for (i=1; i<=nvtxs; i++) fprintf(fout, "%d\n", sets[i]);
  fclose(fout);
}
