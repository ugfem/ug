// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdio.h>

/* Debug break point. */
void checkpnt(tag)
char *tag;
{
  int affirm();
  void exit();

  {char buf[150]; sprintf(buf,"%s: ", tag);UserWrite(buf);}
  if (!affirm("continue")) {exit(0);}
}
