// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* This software was developed by Bruce Hendrickson and Robert Leland   *
* at Sandia National Laboratories under US Department of Energy        *
* contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <math.h>
#include "../main/defs.h"
#include "../main/structs.h"

/* Check orthogonality of vector set */
void checkorth(mat,n,dim)
double **mat;
int n;
int dim;
{
  int i,j;                              /* loop idices */
  double measure;                       /* Froebenius norm */
  double prod;                          /* value of dot product */
  double worst;                         /* greatest off-diagonal dot product */
  int lim;                              /* index of last vec to check against */
  int screenlim;                        /* value of lim that will fit on screen */
  int option;                           /* which option to use */

  double dot();                         /* standard dot product routine */

  /* The T/F argument in the conditionals is just a convenient option: */

  screenlim = 20;
  option = 3;

  /* Check orthogonality over whole set. */
  if (option == 1) {
    {char buf[150]; sprintf(buf,"Orthogonality check:\n");UserWrite(buf);}
    for (i=1; i<=dim; i++)  {
      {char buf[150]; sprintf(buf,"%2d)",i);UserWrite(buf);}
      for (j=1; j<=i; j++)  {
        prod =  dot(mat[i],1,n,mat[j]);
        /*{char buf[150]; sprintf(buf," %g ",prod);UserWrite(buf);} */
        /*{char buf[150]; sprintf(buf," %4.2e ",prod);UserWrite(buf);} */
        /*{char buf[150]; sprintf(buf," %4.2e ",fabs(prod));UserWrite(buf);} */
        {char buf[150]; sprintf(buf," %2d",-(int)log10(prod));UserWrite(buf);}
      }
      {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
    }
  }

  if (option == 2) {
    {char buf[150]; sprintf(buf,"Frobenius orthogonality measure:");UserWrite(buf);}
    measure = 0;
    for (i=1; i<=dim; i++) {
      for (j=i; j<=dim; j++) {
        prod =  dot(mat[i],1,n,mat[j]);
        if (i==j) {
          measure += fabs(1.0 - prod);
        }
        else {
          measure += 2.0 * fabs(prod);
        }
      }
    }
    {char buf[150]; sprintf(buf,"%g \n",measure);UserWrite(buf);}
  }

  /* Check orthogonality against last vector. Allows you
     to build up orthogonality matrix much faster if previous
     columns stay the same when add a new column, but may interact
     with other debug output to give a confusing presentation. */
  if (option == 3) {
    {char buf[150]; sprintf(buf,"%3d) ",dim);UserWrite(buf);}
    lim = min(dim,screenlim);
    worst = 0;
    for (i=1; i<=dim; i++)  {
      prod =  dot(mat[i],1,n,mat[dim]);
      if (i <= lim) {
        {char buf[150]; sprintf(buf," %2d",-(int)log10(fabs(prod)));UserWrite(buf);}
      }
      if ( (i!=dim) && (fabs(prod) > fabs(worst)) ) {
        worst = prod;
      }
    }
    {char buf[150]; sprintf(buf," worst %4.2e\n",worst);UserWrite(buf);}
  }
}
