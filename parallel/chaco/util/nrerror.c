// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>

/* Error handler for Numerical Recipies. */
void nrerror(error_text)
char *error_text;
{
  void exit();

  (void)fprintf(stderr,"Numerical Recipes run-time error...\n");
  (void)fprintf(stderr,"%s\n",error_text);
  (void)fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
