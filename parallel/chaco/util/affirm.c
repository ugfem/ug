// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include "../main/defs.h"

/* Record a return TRUE if answer is yes, FALSE if no. */
int affirm(prompt)
char *prompt;
{
  char response;                /* character typed in */
  int done;                     /* loop control */
  void exit();

  done = FALSE;
  while (!done) {
    {char buf[150]; sprintf(buf,"%s? ", prompt);UserWrite(buf);}
    response = (char) getchar();
    while (response == ' ' || response == '\n') {
      response = (char) getchar();
    }
    if (!strncmp((char *) &response,"y",1)) {return(TRUE);}
    if (!strncmp((char *) &response,"Y",1)) {return(TRUE);}
    if (!strncmp((char *) &response,"n",1)) {return(FALSE);}
    if (!strncmp((char *) &response,"N",1)) {return(FALSE);}
    if (!strncmp((char *) &response,"q",1)) {exit(0);}
    if (!strncmp((char *) &response,"Q",1)) {exit(0);}
    if (!strncmp((char *) &response,"x",1)) {exit(0);}
    if (!strncmp((char *) &response,"X",1)) {exit(0);}
    {char buf[150]; sprintf(buf,"Valid responses begin with: y Y n N q Q x X\n");UserWrite(buf);}
  }

  /* Won't ever get this far, but alint really wants a return value. */
  return(FALSE);
}
