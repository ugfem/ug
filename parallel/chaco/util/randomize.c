// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include        <stdlib.h>
#include        "../main/defs.h"

int RAND_MAXIMUM = 32767;       /* largest value returnable from rand() */

/* Randomly permute elements of an array. */
void randomize(array, n)
int *array;                     /* array of integer values */
int n;                          /* number of values */
{
  double value;                 /* random value */
  int index;                    /* array index to swap with */
  int temp;                     /* holds value being swapped */
  int i;                        /* loop counter */
  double drandom();

  for (i=1; i<=n; i++) {
    value = drandom();
    index = n*value + 1;
    temp = array[i];
    array[i] = array[index];
    array[index] = temp;
  }
}



double drandom()
{
  extern int RAND_MAXIMUM;      /* Largest value rand can return */
  int val;

  val = rand();
  while (val > RAND_MAXIMUM) RAND_MAXIMUM = 2*RAND_MAXIMUM+1;

  return(((double) val) / (1.0+RAND_MAXIMUM));
}


int irandom()
{
  extern int RAND_MAXIMUM;      /* Largest value rand can return */
  int val;

  val = rand();
  while (val > RAND_MAXIMUM) RAND_MAXIMUM = 2*RAND_MAXIMUM+1;

  return(val);
}


void setrandom(seed)
long seed;
{
  int iseed;

  iseed = seed;

  srand((unsigned) iseed);
}
