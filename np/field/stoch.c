// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*																		*/
/* File:	stoch.c														*/
/*																		*/
/* Purpose:	Interface between the fieldgenerator and ug					*/
/*																		*/
/* Author:	Carsten Schwarz												*/
/*		Institut f"ur Hydromechanik und Wasserwirtschaft				*/
/*		ETH H"onggerberg												*/
/*		8093 Z"urich													*/
/*		Schweiz															*/
/*		internet: carsten@teverone.ethz.ch                                                              */
/*																		*/
/* History:   Feb. 1997 begin, ug3-version								*/
/*																		*/
/* Remarks:                                                                                                                     */
/*																		*/
/************************************************************************/

/************************************************************************/
/*																		*/
/* include files														*/
/*	system include files												*/
/*	application include files											*/
/*																		*/
/************************************************************************/

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* low module */
#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "misc.h"
#include "general.h"

/* dev module */
#include "devices.h"
#include "debug.h"

/* np module */
#include "udm.h"
#include "numproc.h"
#include "np.h"
#include "scan.h"
#include "field.h"

/************************************************************************/
/*																		*/
/* defines in the following order										*/
/*																		*/
/*	compile time constants defining static data size (i.e. arrays)		*/
/*	other constants														*/
/*	macros																*/
/*																		*/
/************************************************************************/

#define NumOfOptions    2

#ifdef __THREEDIM__
    #ifdef __TWODIM__
        #error Define only one of __TWODIM__ or __THREEDIM__
    #endif
    #define REALPART(N,Feld,i,j,k) (Feld)[((k)*(N)[1]+(j))*(N)[0]+(i)]
    #define IMAPART(N,Feld,i,j,k)  (Feld)[(((N)[2]+(k))*(N)[1]+(j))*(N)[0]+(i)]
    #define INTEXPDM(i) ((i)*(i)*(i))
#else
    #ifdef __TWODIM__
        #define REALPART(N,Feld,i,j) (Feld)[(j)*(N)[0]+(i)]
        #define IMAPART(N,Feld,i,j)  (Feld)[((N)[1]+(j))*(N)[0]+(i)]
        #define INTEXPDM(i) ((i)*(i))
    #else
        #error Define __TWODIM__ or __THREEDIM__
    #endif
#endif

#ifndef PI
    #define PI (3.1415926535898)
#endif

#define M1      259200
#define M2      134456
#define M3      243000
#define RM1     (1.0 / (DOUBLE) M1)
#define RM2     (1.0 / (DOUBLE) M2)
#define RM3     (1.0 / (DOUBLE) M3)
#define IA1     7141
#define IA2     8121
#define IA3     4561
#define IC1     54773
#define IC2     28411
#define IC3     51349

#define SGN(i) (((i)>0) ? 1 : (((i)==0) ? 0 : -1))

static INT theSeed = 0;
static INT initialize=0;
static DOUBLE *H;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/************************************************************************/
/*																		*/
/* data structures used in this source file (exported data structures	*/
/*	are in the corresponding include file!)								*/
/*																		*/
/************************************************************************/

/************************************************************************/
/*																		*/
/* definition of variables global to this source file only (static!)	*/
/*																		*/
/************************************************************************/

REP_ERR_FILE;

/************************************************************************/
/*																		*/
/* forward declarations of functions used before they are defined		*/
/*																		*/
/************************************************************************/

/************************************************************************/
/*
   getLinRnd - returns a homogen distributed random number

   SYNOPSIS:
   DOUBLE getLinRnd();

   PARAMETERS:
   .  none

   DESCRIPTION:
   getLinRnd calculates a isodistributed random value in the
   intervall [0,1].

   SEE ALSO:
   getGaussRnd

   RETURN VALUE:
   .  DOUBLE random number
 */
/************************************************************************/

static DOUBLE getLinRnd(void)
{
  static DOUBLE RND[97];
  static INT i, ix1, ix2, ix3;

  DOUBLE RetValue;
  if (!initialize)
  {
    ix1 = abs((IC1 - theSeed) % M1);

    ix1 = (IA1 * ix1 + IC1) % M1;
    ix2 =        ix1        % M2;

    ix1 = (IA1 * ix1 + IC1) % M1;
    ix3 =        ix1        % M3;

    for (i=0; i<97; i++)
    {
      ix1 = (IA1 * ix1 + IC1) % M1;
      ix2 = (IA2 * ix2 + IC2) % M2;
      RND[i] = ((DOUBLE) ix2 * RM2 + (DOUBLE) ix1) * RM1;
    }
    initialize = 1;
  }

  ix1 = (IA1 * ix1 + IC1) % M1;
  ix2 = (IA2 * ix2 + IC2) % M2;
  ix3 = (IA3 * ix3 + IC3) % M3;

  i = (97 * ix3) / M3;
  RetValue = RND[i];

  RND[i] = ((DOUBLE) ix2 * RM2 + (DOUBLE) ix1) * RM1;

  return(RetValue);
}
/************************************************************************/
/*
   getGaussRnd - returns a Gauss-distributed random number

   SYNOPSIS:
   DOUBLE getGaussRnd();

   PARAMETERS:
   .  none

   DESCRIPTION:
   getGaussRnd calls a low level random number generator for
   linear distributed random numbers and transforms this values to
   Gauss-distributed random values. One random double number is
   returned.

   SEE ALSO:
   getLinRnd

   RETURN VALUE:
   .  DOUBLE random number
 */
/************************************************************************/

static DOUBLE getGaussRnd(void)
{
  DOUBLE ReLinRnd, ImLinRnd, AbsLinRnd;

  do
  {
    ReLinRnd = 2.0 * getLinRnd() - 1.0;
    ImLinRnd = 2.0 * getLinRnd() - 1.0;
    AbsLinRnd = ReLinRnd*ReLinRnd + ImLinRnd*ImLinRnd;
  } while(AbsLinRnd>1);

  return( ReLinRnd * sqrt(-2.0*log(AbsLinRnd)/AbsLinRnd) );
}

static void getGaussRnd2(DOUBLE *out)
{
  DOUBLE ReLinRnd, ImLinRnd, AbsLinRnd;

  do
  {
    ReLinRnd = 2.0 * getLinRnd() - 1.0;
    ImLinRnd = 2.0 * getLinRnd() - 1.0;
    AbsLinRnd = ReLinRnd*ReLinRnd + ImLinRnd*ImLinRnd;
  } while(AbsLinRnd>1);

  out[0] = ReLinRnd * sqrt(-2.0*log(AbsLinRnd)/AbsLinRnd);
  out[1] = ImLinRnd * sqrt(-2.0*log(AbsLinRnd)/AbsLinRnd);
}

/************************************************************************/
/*
   StochModel - calculates fourier-coefficients for two
   crosscorrelated stochastic fields

   SYNOPSIS:
   INT StochModel(INT i, INT j, INT k, DOUBLE *out,
                       DOUBLE Var, DOUBLE Nug,
                       DOUBLE *Cor, INT actype, DOUBLE *F);

   PARAMETERS:
   .  i,j,k - frequence-vector
   .  out - array for two complex values as output
   .  Var, Nug - mean value, variance, nugget
   .  Cor, F - Array of correlation lengths, basic fourier frequencies
   .  actype - type of autocorrelation function

   DESCRIPTION:
   This function calculates randomly the fourier-coefficients of
   a stochastic field with given stochastic properties
   for one fourier-vector i,j,k.

   In two dimensions the integer k is neglected by the function.

   SEE ALSO:
   stochFourierFld

   RETURN VALUE:
   INT errorcode
   .n zero if none
 */
/************************************************************************/

static INT StochModel(INT i, INT j, INT k, DOUBLE *out, DOUBLE Var, DOUBLE Nug, DOUBLE *Cor, INT actype, DOUBLE *F)
{
  INT l, m, Counter;
  DOUBLE Frequence[DIM], NormOfWavelength, Weight;
  DOUBLE temp;
  DOUBLE Matrix[3];       /* lower left matrix with less than 3 entries, arranged column after column, down from
                             the diagonal */
  DOUBLE RandomNumber[2];

  Frequence[0] = ((DOUBLE) i) * F[0];
  Frequence[1] = ((DOUBLE) j) * F[1];
#ifdef __THREEDIM__
  Frequence[2] = ((DOUBLE) k) * F[2];
#endif

  NormOfWavelength = 0.0;
  for (m=0; m<DIM; m++)
  {
    temp = 2.0*PI*Frequence[m];
    NormOfWavelength += temp*temp*Cor[m]*Cor[m];
  }

  switch(actype)
  {
  case 1 :
    Weight = INTEXPDM(2.0*PI)*Var*Cor[DIM]*exp(-NormOfWavelength/4.0)
             /(((DOUBLE) INTEXPDM(2))*sqrt((DOUBLE) INTEXPDM(PI)) + Nug);
    break;
  case 2 :
    Weight = INTEXPDM(2.0*PI)*Var*Cor[DIM]
             /(pow(PI*(1.0+NormOfWavelength),((DOUBLE) DIM + 1.0)/2.0) + Nug);
    break;
  default :
    return(1);
  }

  Matrix[0] = sqrt(Weight);          /* Matrixentry 1 1 */
  Matrix[1] = 0.0;                      /* Matrixentry 2 1 */

  Matrix[2] = Matrix[0];                /*             2 2 */

  for (l=0; l<2; l++)
  {
    RandomNumber[l] = sqrt(0.5) * getGaussRnd();
    out[l] = 0.0;
  }


  Counter = 0;
  for (l=0; l<2; l++)
    for(m=l; m<2; m++)
      out[m] += Matrix[Counter++]*RandomNumber[l];                   /* think abaut the scaling, some factors
                                                                        may vanish due to the scaling of the fourier space */

  return(0);
}

/************************************************************************/
/*
   DoBetrag - Calculates a real number as fourier coefficient

   SYNOPSIS:
   INT DoBetrag(DOUBLE *Feld, INT i, INT j, INT k, INT *NOfN);
   INT DoBetrag(DOUBLE *Feld, INT i, INT j, INT *NOfN)

   PARAMETERS:
   .  Feld - array where the values are
   .  i,j,k - frequence-vector
   .  NOfN - array with dimensional information

   DESCRIPTION:
   Calculates a real number for a given frequence as fourier
   coefficient, the number should have the same mean absolut value
   and the same variance as the complex number originally calculated

   SEE ALSO:
   stochFourierFld

   RETURN VALUE:
   INT errorcode
   .n zero if none
 */
/************************************************************************/

#ifdef __THREEDIM__
static INT DoBetrag(DOUBLE *Feld, INT i, INT j, INT k, INT *NOfN)
{
  DOUBLE Value;

  Value = sqrt(2.0) * REALPART(NOfN, Feld, i, j, k);
  REALPART(NOfN, Feld, i, j, k) = Value; IMAPART(NOfN, Feld, i, j, k) = 0.0;

  return(0);
}
#else
static INT DoBetrag(DOUBLE *Feld, INT i, INT j, INT *NOfN)
{
  DOUBLE Value;

  Value = sqrt(2.0) * REALPART(NOfN, Feld, i, j);
  REALPART(NOfN, Feld, i, j) = Value; IMAPART(NOfN, Feld, i, j) = 0.0;

  return(0);
}
#endif

/************************************************************************/
/*
   stochFourierFld - calculates a stochastic field in
   fourier space

   SYNOPSIS:
   INT stochFourierFld(INT *NOfN, DOUBLE Var,
                       DOUBLE Nug, DOUBLE *Cor, INT actype, DOUBLE *F);

   PARAMETERS:
   .  N - array with dimensional information
   .  NField - Number of generated fields (one or two)
   .  Var, Nug, Cor, actype - stochastic properties,
      that should be generated
   .  F - highest fourier frequencies

   DESCRIPTION:
   This function gives the fourier coefficients of a stochastic field.
   The properties such as mean value, variance, nuggets (for the
   autocorrelation function), correlation length and type of
   autocorrelation function.

   This procedure loops over all fourier frequencies of the array of
   size NOfN[0]*NOfN[1]*NOfN[2]. For each frequence it invokes the
   function StochModel to get random-numbers with frequence-dependend
   distribution. These values are calculated such that the
   non-fourier-values will have the desired stochastical properties
   and are real values.

   The mean value is not incorperated in this function, it will be
   added to the backwards-fourier-transformed values afterwards.

   This function needs the mean value, the nugget, the variance,
   DIM correlation lengths and the autocorrelation type, as well as
   the highest fourier frequencies of the array.

   The autocorrelation types are:
   .  1 - bell-shaped anisotropic covariance
   .  2 - Exponential anisotropic covariance
   For both coherency functions the correlation length is the
   scaling vector of the coordinate axis.

   SEE ALSO:
   StochModel

   RETURN VALUE:
   INT errorcode
   .n zero if none
 */
/************************************************************************/

static INT stochFourierFld(INT *NOfN, DOUBLE Var,
                           DOUBLE Nug, DOUBLE *Cor, INT actype, DOUBLE *F)
{

  INT i,j;
#ifdef __THREEDIM__
  INT k;
#endif
  DOUBLE out[4];

#ifdef __THREEDIM__
  for (k = 1; k <= NOfN[2]/2; k++)
  {
    for(j = 1; j <= NOfN[1]/2; j++)
    {
      for(i = 1; i <= NOfN[0]/2; i++)
      {
        StochModel(i, j, k, out, Var, Nug, Cor, actype, F);
        REALPART(NOfN,H,         i,         j,         k) =  out[0];
        IMAPART(NOfN,H,         i,         j,         k) =  out[1];
        REALPART(NOfN,H, NOfN[0]-i, NOfN[1]-j, NOfN[2]-k) =  out[0];
        IMAPART(NOfN,H, NOfN[0]-i, NOfN[1]-j, NOfN[2]-k) = -out[1];

        StochModel(i, -j, k, out, Var, Nug, Cor, actype, F);
        REALPART(NOfN,H,        i , NOfN[1]-j,         k) =  out[0];
        IMAPART(NOfN,H,         i, NOfN[1]-j,         k) =  out[1];
        REALPART(NOfN,H, NOfN[0]-i,         j, NOfN[2]-k) =  out[0];
        IMAPART(NOfN,H, NOfN[0]-i,         j, NOfN[2]-k) = -out[1];

        StochModel(i, j, -k, out, Var, Nug, Cor, actype, F);
        REALPART(NOfN,H,         i,         j, NOfN[2]-k) =  out[0];
        IMAPART(NOfN,H,         i,         j, NOfN[2]-k) =  out[1];
        REALPART(NOfN,H, NOfN[0]-i, NOfN[1]-j,         k) =  out[0];
        IMAPART(NOfN,H, NOfN[0]-i, NOfN[1]-j,         k) = -out[1];

        StochModel(i, -j, -k, out, Var, Nug, Cor, actype, F);
        REALPART(NOfN,H,         i, NOfN[1]-j, NOfN[2]-k) =  out[0];
        IMAPART(NOfN,H,         i, NOfN[1]-j, NOfN[2]-k) =  out[1];
        REALPART(NOfN,H, NOfN[0]-i,         j,         k) =  out[0];
        IMAPART(NOfN,H, NOfN[0]-i,         j,         k) = -out[1];
      }

      StochModel(0, j, k, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H, 0,         j,         k) =  out[0];
      IMAPART(NOfN,H, 0,         j,         k) =  out[1];
      REALPART(NOfN,H, 0, NOfN[1]-j, NOfN[2]-k) =  out[0];
      IMAPART(NOfN,H, 0, NOfN[1]-j, NOfN[2]-k) = -out[1];

      StochModel(0, -j, k, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H, 0, NOfN[1]-j,         k) =  out[0];
      IMAPART(NOfN,H, 0, NOfN[1]-j,         k) =  out[1];
      REALPART(NOfN,H, 0,         j, NOfN[2]-k) =  out[0];
      IMAPART(NOfN,H, 0,         j, NOfN[2]-k) = -out[1];
    }

    for (i = 1; i <= NOfN[0]/2; i++)
    {
      StochModel( i, 0, k, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H,         i, 0,         k) =  out[0];
      IMAPART(NOfN,H,         i, 0,         k) =  out[1];
      REALPART(NOfN,H, NOfN[0]-i, 0, NOfN[2]-k) =  out[0];
      IMAPART(NOfN,H, NOfN[0]-i, 0, NOfN[2]-k) = -out[1];

      StochModel( i, 0, -k, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H,         i, 0, NOfN[2]-k) =  out[0];
      IMAPART(NOfN,H,         i, 0, NOfN[2]-k) =  out[1];
      REALPART(NOfN,H, NOfN[0]-i, 0,         k) =  out[0];
      IMAPART(NOfN,H, NOfN[0]-i, 0,         k) = -out[1];
    }

    StochModel( 0, 0, k, out, Var, Nug, Cor, actype, F);
    REALPART(NOfN,H, 0, 0,         k) =  out[0];
    IMAPART(NOfN,H, 0, 0,         k) =  out[1];
    REALPART(NOfN,H, 0, 0, NOfN[2]-k) =  out[0];
    IMAPART(NOfN,H, 0, 0, NOfN[2]-k) = -out[1];
  }

  for(j = 1; j <= NOfN[1]/2; j++)
  {
    for(i = 1; i <= NOfN[0]/2; i++)
    {
      StochModel( i, j, 0, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H,         i,         j, 0) =  out[0];
      IMAPART(NOfN,H,         i,         j, 0) =  out[1];
      REALPART(NOfN,H, NOfN[0]-i, NOfN[1]-j, 0) =  out[0];
      IMAPART(NOfN,H, NOfN[0]-i, NOfN[1]-j, 0) = -out[1];

      StochModel( i, -j, 0, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H,         i, NOfN[1]-j, 0) =  out[0];
      IMAPART(NOfN,H,         i, NOfN[1]-j, 0) =  out[1];
      REALPART(NOfN,H, NOfN[0]-i,         j, 0) =  out[0];
      IMAPART(NOfN,H, NOfN[0]-i,         j, 0) = -out[1];
    }

    StochModel( 0, j, 0, out, Var, Nug, Cor, actype, F);
    REALPART(NOfN,H, 0,         j, 0) =  out[0];
    IMAPART(NOfN,H, 0,         j, 0) =  out[1];
    REALPART(NOfN,H, 0, NOfN[1]-j, 0) =  out[0];
    IMAPART(NOfN,H, 0, NOfN[1]-j, 0) = -out[1];
  }

  for(i = 1; i <= NOfN[0]/2; i++)
  {
    StochModel( i, 0, 0, out, Var, Nug, Cor, actype, F);
    REALPART(NOfN,H,         i, 0, 0) =  out[0];
    IMAPART(NOfN,H,         i, 0, 0) =  out[1];
    REALPART(NOfN,H, NOfN[0]-i, 0, 0) =  out[0];
    IMAPART(NOfN,H, NOfN[0]-i, 0, 0) = -out[1];
  }

  REALPART(NOfN,H, 0, 0, 0) =  0.0;
  IMAPART(NOfN,H, 0, 0, 0) =  0.0;

  DoBetrag(H, NOfN[0]/2, NOfN[1]/2, NOfN[2]/2, NOfN);
  DoBetrag(H,         0, NOfN[1]/2, NOfN[2]/2, NOfN);
  DoBetrag(H, NOfN[0]/2,         0, NOfN[2]/2, NOfN);
  DoBetrag(H, NOfN[0]/2, NOfN[1]/2,         0, NOfN);
  DoBetrag(H,         0,         0, NOfN[2]/2, NOfN);
  DoBetrag(H,         0, NOfN[1]/2,         0, NOfN);
  DoBetrag(H, NOfN[0]/2,         0,         0, NOfN);

#else
  for(j = 1; j <= NOfN[1]/2; j++)
  {
    for(i = 1; i <= NOfN[0]/2; i++)
    {
      StochModel( i, j, 0, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H,         i,         j) =  out[0];
      IMAPART(NOfN,H,         i,         j) =  out[1];
      REALPART(NOfN,H, NOfN[0]-i, NOfN[1]-j) =  out[0];
      IMAPART(NOfN,H, NOfN[0]-i, NOfN[1]-j) = -out[1];

      StochModel( i, -j, 0, out, Var, Nug, Cor, actype, F);
      REALPART(NOfN,H,         i, NOfN[1]-j) =  out[0];
      IMAPART(NOfN,H,         i, NOfN[1]-j) =  out[1];
      REALPART(NOfN,H, NOfN[0]-i,         j) =  out[0];
      IMAPART(NOfN,H, NOfN[0]-i,         j) = -out[1];
    }

    StochModel( 0, j, 0, out, Var, Nug, Cor, actype, F);
    REALPART(NOfN,H, 0,         j) =  out[0];
    IMAPART(NOfN,H, 0,         j) =  out[1];
    REALPART(NOfN,H, 0, NOfN[1]-j) =  out[0];
    IMAPART(NOfN,H, 0, NOfN[1]-j) = -out[1];

  }

  for(i = 1; i <= NOfN[0]/2; i++)
  {
    StochModel( i, 0, 0, out, Var, Nug, Cor, actype, F);
    REALPART(NOfN,H,         i, 0) =  out[0];
    IMAPART(NOfN,H,         i, 0) =  out[1];
    REALPART(NOfN,H, NOfN[0]-i, 0) =  out[0];
    IMAPART(NOfN,H, NOfN[0]-i, 0) = -out[1];
  }

  REALPART(NOfN,H, 0, 0) =  0.0;
  IMAPART(NOfN,H, 0, 0) =  0.0;

  DoBetrag(H, NOfN[0]/2, NOfN[1]/2, NOfN);
  DoBetrag(H,         0, NOfN[1]/2, NOfN);
  DoBetrag(H, NOfN[0]/2,         0, NOfN);

#endif

  return(0);
}

/************************************************************************/
/*
   comAdd, comSub, comMult - some pseudo complex number operations

   SYNOPSIS:
   void comAdd(DOUBLE *Summand1, DOUBLE *Summand2, DOUBLE *Summe);
   void comSub(DOUBLE *Summand1, DOUBLE *Summand2, DOUBLE *Summe);
   void comMult(DOUBLE *Multiplikator, DOUBLE *Multiplikant, DOUBLE *Produkt);

   PARAMETERS:
   .  Summand1, Summand2, Multiplikator, Multiplikant - array of two DOUBLE variables
   describing one complex number as input values
   .  Summe, Produkt - array of two DOUBLE variables for output

   DESCRIPTION:
   This functions compute the complex sum, difference and product of two
   complex numbers and put that value into the third array of the commandline.
   The user of these functions has to control, that the pointer Summe or
   Produkt points to the memory for two DOUBLE numbers.

   SEE ALSO:
   comZuweisung, EhochiPhi

   RETURN VALUE:
   void
 */
/************************************************************************/

static void comAdd(DOUBLE *Summand1, DOUBLE *Summand2, DOUBLE *Summe)
{
  Summe[0] = Summand1[0] + Summand2[0];
  Summe[1] = Summand1[1] + Summand2[1];
}

static void comSub(DOUBLE *Summand1, DOUBLE *Summand2, DOUBLE *Summe)
{
  Summe[0] = Summand1[0] - Summand2[0];
  Summe[1] = Summand1[1] - Summand2[1];
}

static void comMult(DOUBLE *Multiplikator, DOUBLE *Multiplikant, DOUBLE *Produkt)
{
  Produkt[0] = Multiplikator[0] * Multiplikant[0] - Multiplikator[1] * Multiplikant[1];
  Produkt[1] = Multiplikator[0] * Multiplikant[1] + Multiplikator[1] * Multiplikant[0];
}

/************************************************************************/
/*
   EhochiPhi - the complex exponential function with imaginary arguments

   SYNOPSIS:
   void EhochiPhi(DOUBLE Phi, DOUBLE *out);

   PARAMETERS:
   .  Phi - the argument of the exponential function divided by i
   .  out - pointer to the target location of the complex number

   DESCRIPTION:
   This function computes exp(i*Phi) and returns that value in the
   variable out.

   SEE ALSO:
   comAdd, comSub, comMult, comZuweisung

   RETURN VALUE:
   void
 */
/************************************************************************/

static void EhochiPhi(DOUBLE Phi, DOUBLE *out)
{
  out[0] = cos(Phi);
  out[1] = sin(Phi);
}

/************************************************************************/
/*
   comZuweisung - copy one complex number to another variable

   SYNOPSIS:
   void comZuweisung(DOUBLE *in, DOUBLE *out);

   PARAMETERS:
   .  in - array of two DOUBLE variables
   describing one complex number as input value
   .  out - pointer to the target location of the complex number

   DESCRIPTION:
   This function just copies the entries in[0] and in[1] to the variables
   out[0] and out[1].

   SEE ALSO:
   comAdd, comSub, comMult, EhochiPhi

   RETURN VALUE:
   void
 */
/************************************************************************/

static void comZuweisung(DOUBLE *in, DOUBLE *out)
{
  out[0] = in[0];
  out[1] = in[1];
}

/************************************************************************/
/*
   fFT - fast fourier transform

   SYNOPSIS:
   INT fFT(DOUBLE *H, INT *N)

   PARAMETERS:
   .  H - array, at start fourier coefficients, after fFT retransformed values
   .  N - array with dimensional information

   DESCRIPTION:
   This procedure is the fourier backwards transform for the generator of
   stochastic fields. The normation in this procedure is adequate for the
   normation of the stochastic fourier components.

   The sidelengths of the transformed field have to be a power of 2.

   The indices of the calculated values are in bitwise invers order.
   The index calculation can be done with the function Rho.

   SEE ALSO:
   stochFourierFld, Rho

   RETURN VALUE:
   INT errorcode
   .n zero if none
 */
/************************************************************************/

static INT fFT(DOUBLE *Field, INT *N)
{
#ifdef __THREEDIM__
  INT Slice1, Slice2;
#else
  INT Slice;
#endif
  INT k,l,m,n;
  DOUBLE T[2],U[2],V[2],W[2],X[2],Z[2];

#ifdef __THREEDIM__
  for (Slice1=0; Slice1<N[1]; Slice1++)
    for (Slice2=0; Slice2<N[2]; Slice2++)
#else
  for (Slice=0; Slice<N[1]; Slice++)
#endif
    {
      EhochiPhi(PI * 2.0 / ((DOUBLE) (N[0])),T);
      for(n=N[0]; n>1; n=m)
      {
        m = n/2;
        for(k=0; k<N[0]; k+=n)
        {
          W[0] = 1.0;
          W[1] = 0.0;
          for(l=0; l<m; l++)
          {
#ifdef __THREEDIM__
            U[0] = REALPART(N, Field, l+k  , Slice1, Slice2);
            U[1] =  IMAPART(N, Field, l+k  , Slice1, Slice2);
            V[0] = REALPART(N, Field, l+k+m, Slice1, Slice2);
            V[1] =  IMAPART(N, Field, l+k+m, Slice1, Slice2);
#else
            U[0] = REALPART(N, Field, l+k  , Slice);
            U[1] =  IMAPART(N, Field, l+k  , Slice);
            V[0] = REALPART(N, Field, l+k+m, Slice);
            V[1] =  IMAPART(N, Field, l+k+m, Slice);
#endif
            comAdd(U,V,Z);
#ifdef __THREEDIM__
            REALPART(N, Field, l+k  , Slice1, Slice2) = Z[0];
            IMAPART(N, Field, l+k  , Slice1, Slice2) = Z[1];
#else
            REALPART(N, Field, l+k  , Slice) = Z[0];
            IMAPART(N, Field, l+k  , Slice) = Z[1];
#endif
            comSub(U,V,Z);
            comMult(Z,W,X);
#ifdef __THREEDIM__
            REALPART(N, Field, l+k+m, Slice1, Slice2) = X[0];
            IMAPART(N, Field, l+k+m, Slice1, Slice2) = X[1];
#else
            REALPART(N, Field, l+k+m, Slice) = X[0];
            IMAPART(N, Field, l+k+m, Slice) = X[1];
#endif
            comMult(W,T,Z);
            comZuweisung(Z,W);
          }
        }
        comMult(T,T,Z);
        comZuweisung(Z,T);
      }
    }

#ifdef __THREEDIM__
  for (Slice1=0; Slice1<N[0]; Slice1++)
    for (Slice2=0; Slice2<N[2]; Slice2++)
#else
  for (Slice=0; Slice<N[0]; Slice++)
#endif
    {
      EhochiPhi(PI * 2.0 / ((DOUBLE) (N[1])),T);
      for(n=N[1]; n>1; n=m)
      {
        m = n/2;
        for(k=0; k<N[1]; k+=n)
        {
          W[0] = 1.0;
          W[1] = 0.0;
          for(l=0; l<m; l++)
          {
#ifdef __THREEDIM__
            U[0] = REALPART(N, Field, Slice1, l+k  , Slice2);
            U[1] =  IMAPART(N, Field, Slice1, l+k  , Slice2);
            V[0] = REALPART(N, Field, Slice1, l+k+m, Slice2);
            V[1] =  IMAPART(N, Field, Slice1, l+k+m, Slice2);
#else
            U[0] = REALPART(N, Field, Slice, l+k  );
            U[1] =  IMAPART(N, Field, Slice, l+k  );
            V[0] = REALPART(N, Field, Slice, l+k+m);
            V[1] =  IMAPART(N, Field, Slice, l+k+m);
#endif
            comAdd(U,V,Z);
#ifdef __THREEDIM__
            REALPART(N, Field, Slice1, l+k  , Slice2) = Z[0];
            IMAPART(N, Field, Slice1, l+k  , Slice2) = Z[1];
#else
            REALPART(N, Field, Slice, l+k  ) = Z[0];
            IMAPART(N, Field, Slice, l+k  ) = Z[1];
#endif
            comSub(U,V,Z);
            comMult(Z,W,X);
#ifdef __THREEDIM__
            REALPART(N, Field, Slice1, l+k+m, Slice2) = X[0];
            IMAPART(N, Field, Slice1, l+k+m, Slice2) = X[1];
#else
            REALPART(N, Field, Slice, l+k+m) = X[0];
            IMAPART(N, Field, Slice, l+k+m) = X[1];
#endif
            comMult(W,T,Z);
            comZuweisung(Z,W);
          }
        }
        comMult(T,T,Z);
        comZuweisung(Z,T);
      }
    }
#ifdef __THREEDIM__
  for (Slice1=0; Slice1<N[0]; Slice1++)
    for (Slice2=0; Slice2<N[1]; Slice2++)
    {
      EhochiPhi(PI * 2.0 / ((DOUBLE) (N[2])),T);
      for(n=N[2]; n>1; n=m)
      {
        m = n/2;
        for(k=0; k<N[2]; k+=n)
        {
          W[0] = 1.0;
          W[1] = 0.0;
          for(l=0; l<m; l++)
          {
            U[0] = REALPART(N, Field, Slice1, Slice2, l+k  );
            U[1] =  IMAPART(N, Field, Slice1, Slice2, l+k  );
            V[0] = REALPART(N, Field, Slice1, Slice2, l+k+m);
            V[1] =  IMAPART(N, Field, Slice1, Slice2, l+k+m);
            comAdd(U,V,Z);
            REALPART(N, Field, Slice1, Slice2, l+k  ) = Z[0];
            IMAPART(N, Field, Slice1, Slice2, l+k  ) = Z[1];
            comSub(U,V,Z);
            comMult(Z,W,X);
            REALPART(N, Field, Slice1, Slice2, l+k+m) = X[0];
            IMAPART(N, Field, Slice1, Slice2, l+k+m) = X[1];
            comMult(W,T,Z);
            comZuweisung(Z,W);
          }
        }
        comMult(T,T,Z);
        comZuweisung(Z,T);
      }
    }
#endif

  return(0);
}

/************************************************************************/
/*
   Rho - change the bitorder of an integer number

   SYNOPSIS:
   INT Rho(INT Index, INT len);

   PARAMETERS:
   .  Index - an integer number
   .  len - length of bitfield in bit

   DESCRIPTION:
   This function computes a number with the bits in Index read in invers order
   in the lowest len bits. This function is used for index-arithmetics for the
   fast-fourier-transformation.

   SEE ALSO:
   fFT

   RETURN VALUE:
   INT
   .n The invers order bitfield.
 */
/************************************************************************/

static INT Rho(INT Index, INT len, INT n)
{
  INT Value;
  INT left, right;
  INT shift;

  shift = len;
  Value = 0;
  for (left=1, right=n/2, shift--; shift>=0; left<<=1, right>>=1, shift-=2)
  {
    Value |= (Index & left)<<shift;
    Value |= (Index & right)>>shift;
  }
  return(Value);
}

/************************************************************************/
/*
   CopyTo - Copy the realpart of the input field to the output field

   SYNOPSIS:
   void CopyTo(DOUBLE *inField, INT *NOfN, DOUBLE *outField, INT NField);

   PARAMETERS:
   .  inField - a complex field, filled by the stochastic generator as input
   .  NOfN - information about array size
   .  outField - pointer to the output field

   DESCRIPTION:
   This procedure copies the realparts of the value of the input field
   inField to the output field outField. During this operation the
   index correction for the fast fourier transform is included.

   SEE ALSO:
   fFT, Field_genStochField

   RETURN VALUE:
   void
 */
/************************************************************************/

static void CopyTo(DOUBLE *inField, INT *NOfN, DOUBLE *outField)
{
  INT log2[DIM];
#ifdef __THREEDIM__
  INT i,j,k;
#else
  INT i,j;
#endif

  for (i=0; i<DIM; i++)
    for(log2[i]=0, j=1; j<NOfN[i]; log2[i] +=1)
      j*=2;

  for(i=0; i<NOfN[0]; i++)
    for(j=0; j<NOfN[1]; j++)
#ifdef __THREEDIM__
      for (k = 0; k < NOfN[2]; k++)
        outField[(k*NOfN[1]+j)*NOfN[0]+i]
          = REALPART(NOfN, inField, Rho(i,log2[0],NOfN[0]), Rho(j,log2[1],NOfN[1]), Rho(k,log2[2],NOfN[2]));
#else
      outField[j*NOfN[0]+i]
        = REALPART(NOfN, inField, Rho(i,log2[0],NOfN[0]), Rho(j,log2[1],NOfN[1]));
#endif
}

/************************************************************************/
/*
   add - add a value to each entry of an array

   SYNOPSIS:
   void add(DOUBLE *Feld, INT *N, DOUBLE value);

   PARAMETERS:
   .  Feld - pointer to DOUBLE, first element of array
   .  NOfN - information about array size
   .  value - constant to add

   DESCRIPTION:
   This procedure adds the value value to each entry of the array Feld.

   This procedure is invoked by Field_genStochField to add the desired
   mean value to the stochastic field with mean value zero if
   Field_genStochField was called without the correct option. Otherwise
   this procedure will be called by correct so, that the field afterwards
   has the desired mean value.

   SEE ALSO:
   mult, correct, Field_genStochField

   RETURN VALUE:
   void
 */
/************************************************************************/

static void add(DOUBLE *Feld, INT *NOfN, DOUBLE value)
{
#ifdef __THREEDIM__
  int i,j,k;
#else
  int i,j;
#endif

  for (i = 0; i < NOfN[0]; i++)
    for (j = 0; j < NOfN[1]; j++)
#ifdef __THREEDIM__
      for (k = 0; k < NOfN[2]; k++)
        REALPART(NOfN, Feld ,i, j, k) += value;
#else
      REALPART(NOfN, Feld, i, j) += value;
#endif
}

/************************************************************************/
/*
   mult - multiply each entry of an array with a number

   SYNOPSIS:
   void mult(DOUBLE *Feld, INT *N, DOUBLE factor);

   PARAMETERS:
   .  Feld - pointer to DOUBLE, first element of array
   .  NOfN - information about array size
   .  factor - number for multiplication

   DESCRIPTION:
   This procedure multiplies each entry of the array Feld with
   the number factor.

   It is invoked by the procedure correct to archive the right
   variance of the array.

   SEE ALSO:
   add, correct

   RETURN VALUE:
   void
 */
/************************************************************************/

static void mult (DOUBLE *Feld, INT *NOfN, DOUBLE factor)
{
#ifdef __THREEDIM__
  int i,j,k;
#else
  int i,j;
#endif

  for (i = 0; i < NOfN[0]; i++)
    for (j = 0; j < NOfN[1]; j++)
#ifdef __THREEDIM__
      for (k = 0; k < NOfN[2]; k++)
        REALPART(NOfN, Feld ,i, j, k) *= factor;
#else
      REALPART(NOfN, Feld, i, j) *= factor;
#endif
}

/************************************************************************/
/*
   correct - corrects the mean and the variance of an array

   SYNOPSIS:
   void correct(DOUBLE *Feld, INT *NOfN, DOUBLE wantedMean, DOUBLE wantedVariance);

   PARAMETERS:
   .  Feld - pointer to DOUBLE, first element of array
   .  NOfN - information about array size
   .  wantedMean - the wanted mean value of the array
   .  wantedVariance - the wanted variance of the array

   DESCRIPTION:
   This procedure multiplies an array with a constant and adds another
   value to each entry of the array to get an array with predefined
   mean value and variance. Therefore this procedure first computes
   the mean and the variance of the array and than shifts and multiplies
   it such that it gets the right mean and variance.

   It is invoked by Field_genStochField if the correct option is
   set after the computation of the stochastic field.

   SEE ALSO:
   add, mult, Field_genStochField

   RETURN VALUE:
   void
 */
/************************************************************************/

static void correct (DOUBLE *Feld, INT *NOfN, DOUBLE wantedMean, DOUBLE wantedVariance)
{
  INT i,j, Anzahl;
#ifdef __THREEDIM__
  INT k;
#endif
  DOUBLE realMean, realVar;
  DOUBLE addition, factor;

  realMean = 0.0; realVar = 0.0;

#ifdef __THREEDIM__
  Anzahl = NOfN[0] * NOfN[1] * NOfN[2];
  for (i = 0; i < NOfN[0]; i++)
    for (j = 0; j < NOfN[1]; j++)
      for (k = 0; k < NOfN[2]; k++)
      {
        realMean += REALPART(NOfN, Feld, i, j, k);
        realVar  += REALPART(NOfN, Feld, i, j, k) * REALPART(NOfN, Feld, i, j, k);
      }
#else
  Anzahl = NOfN[0] * NOfN[1];
  for (i = 0; i < NOfN[0]; i++)
    for (j = 0; j < NOfN[1]; j++)
    {
      realMean += REALPART(NOfN, Feld, i, j);
      realVar  += REALPART(NOfN, Feld, i, j) * REALPART(NOfN, Feld, i, j);
    }
#endif

  realMean /= Anzahl;
  realVar  /= Anzahl;
  realVar  -= realMean * realMean;

  factor = sqrt(wantedVariance/realVar);
  addition = wantedMean - factor*realMean;

  mult(Feld, NOfN, factor);
  add(Feld, NOfN, addition);
}

/************************************************************************/
/*
   Field_genStochField - calculate an array with prescribed stochastic
   properties

   SYNOPSIS:
   INT Field_genSochField(NP_STOCH_FIELD *np);

   PARAMETERS:
   .  np - numproc with properties of the field and pointer to array for the
        field

   DESCRIPTION:
   This procedure is called by the NPStochFieldInit to get
   an array with prescribed stochastical properties.

   At first this procedure allocates the memory for the working array.
   Than it calls stochFourierFld to compute the fourier coefficients
   of the desired array. This values are fourier-transform, shifted to the
   mean value and depending on the options corrected.

   The units of space are measured in numbers of node points. This
   includes, that we have an isotropic , regular grid for the generation
   of the stochastic array.

   SEE ALSO:
   add, correct, stochFourierFld, fFt, NPStochFieldInit

   RETURN VALUE:
   INT error code
   .n zero if none
 */
/************************************************************************/

INT Field_genStochField(NP_STOCH_FIELD *np)
{
  DOUBLE *FieldH;
  DOUBLE F[DIM], Cor[DIM+1];       /* the last one is the product of the others */
  INT i, NOfN[DIM+1];
  MULTIGRID *theMg;
  HEAP *theHeap;
  INT MarkKey;

  theMg = NP_MG(np);
  theHeap = MGHEAP(theMg);
  MarkTmpMem(theHeap,&MarkKey);

  NOfN[DIM] = 1;
  Cor[DIM] = 1.0;
  for (i=0; i<DIM; i++)
  {
    NOfN[DIM] *= (NOfN[i] = np->size[i]);
    F[i] = 1.0/NOfN[i];
    Cor[DIM] *= (Cor[i] = np->cor[i]/np->cs[i]);
  }
  theSeed = np->initial;
  initialize = 0;

  if ((FieldH = (DOUBLE *) GetTmpMem(theHeap, 2*NOfN[DIM]*sizeof(DOUBLE),MarkKey)) == NULL)
    return(1);
  H = FieldH;

  stochFourierFld(NOfN, np->var, np->nugget,
                  Cor, np->actype, F);

  fFT(H, &(np->size[0]));

  mult(H, &(np->size[0]), 1.0/sqrt((DOUBLE) NOfN[DIM]));
  add(H, &(np->size[0]),np->mean);

  CopyTo(H, &(np->size[0]),np->Fld);

  ReleaseTmpMem(theHeap,MarkKey);
  return(0);
}

static INT ReadArgvINTVec (const char *name,  INT *ivals, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  int ix, iy, iz;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %d %d %d",option,&ix,&iy,&iz)!=DIM+1)
        continue;
      if (strcmp(option,name) == 0)
      {
        ivals[0] = ix;
        ivals[1] = iy;
#ifdef __THREEDIM__
        ivals[2] = iz;
#endif
        return(0);
      }
    }

  REP_ERR_RETURN(1);
}

static INT testPow2(INT number)
{
  INT compare;

  for(compare=1; compare<number; compare <<= 1)
  {}

  return ((compare==number));
}

/****************************************************************************/
/*D
        NP_STOCH_FIELD - type definition for stochastic fields

   DESCRIPTION:
   This numproc type is used for the generation and storage of
   stochastic fields. It must be initialized by

   'INT NPStochFieldInit (NP_STOCH_FIELD *theNP, INT argc , char **argv);'

   During initialization the stochastic field is allocated and generated.
   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPStochFieldDisplay (NP_STOCH_FIELD *theNP);'
   'INT NPStochDataExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   npinit $s <size> [<size2> <size3>] [$m <mean>] [$v <var>] [$n <nugget>]
       [$c <cor> [<cor2> <cor3>]] [$d <dist> [<dist2> <dist3>]
       [$e|$g] [$i <initial>] [$lin|$const]

   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

static INT NPStochFieldInit(NP_BASE *theNP, INT argc , char **argv)
{
  NP_STOCH_FIELD *np;
  MULTIGRID *theMG;
  HEAP *theHeap;
  INT i, NumOfDouble, sopt, ival[DIM], ret;
  DOUBLE dval[DIM];
  DOUBLE *theField;

  np = (NP_STOCH_FIELD *) theNP;
  theMG = NP_MG(theNP);
  if(theMG == NULL)
    return(NP_NOT_ACTIVE);
  theHeap = theMG->theHeap;

  ret = NP_ACTIVE;
  sopt = 0;
  if (ReadArgvINTVec("s",ival,argc,argv))
    if (ReadArgvINT("s",ival,argc,argv))
    {
      for (i=0; i<DIM; i++)
        if (np->size[i] <= 0)
          ret = NP_NOT_ACTIVE;
    }
    else
    if ((ival[0] <= 0) || (! testPow2(ival[0])))
    {
      PrintErrorMessage('E',"NPStochFieldInit","size must be a power of 2");
      ret = NP_NOT_ACTIVE;
    }
    else
    {
      for (i=0; i<DIM; i++)
        if (np->size[i] != ival[0])
        {
          np->size[i] = ival[0];
          sopt = 1;
        }
    }
  else
    for (i=0; i<DIM; i++)
      if ((ival[i] > 0) && (testPow2(ival[1])))
        if (np->size[i] != ival[i])
        {
          np->size[i] = ival[i];
          sopt = 1;
        }
        else
        {
          PrintErrorMessage('E',"NPStochFieldInit","size must be a power of 2");
          ret = NP_NOT_ACTIVE;
        }

  if (ReadArgvDOUBLE("m",dval,argc,argv))
  {
    if (np->mean == 0.0)
      ret = NP_NOT_ACTIVE;
  }
  else
  if (dval[0] == 0.0)
  {
    PrintErrorMessage('E',"NPStochFieldInit","vanishing mean");
    ret = NP_NOT_ACTIVE;
  }
  else
    np->mean = dval[0];

  if (ReadArgvDOUBLE("v",dval,argc,argv))
  {
    if (np->var < 0.0)
      ret = NP_NOT_ACTIVE;
  }
  else
  if (dval[0] < 0.0)
  {
    PrintErrorMessage('E',"NPStochFieldInit","negative variance");
    ret = NP_NOT_ACTIVE;
  }
  else
    np->var = dval[0];

  if (!(ReadArgvDOUBLE("n",dval,argc,argv)))
  {
    if (dval[0] < 0.0)
    {
      PrintErrorMessage('E',"NPStochFieldInit","negative nugget");
      ret = NP_NOT_ACTIVE;
    }
    else
      np->var = dval[0];
  }

  if (ReadArgvPosition("c",argc,argv,dval))
    if (ReadArgvDOUBLE("c",dval,argc,argv))
    {
      for (i=0; i<DIM; i++)
        if (np->cor[i] <= 0.0)
          ret = NP_NOT_ACTIVE;
    }
    else
    if (dval[0] <= 0.0)
    {
      PrintErrorMessage('E',"NPStochFieldInit","correlation must be positiv");
      ret = NP_NOT_ACTIVE;
    }
    else
      for (i=0; i<DIM; i++)
        np->cor[i] = dval[0];
  else
    for (i=0; i<DIM; i++)
      if (dval[i] > 0.0)
        np->cor[i] = dval[i];
      else
      {
        PrintErrorMessage('E',"NPStochFieldInit","correlation must be positiv");
        ret = NP_NOT_ACTIVE;
      }

  if (ReadArgvPosition("d",argc,argv,dval))
    if (ReadArgvDOUBLE("d",dval,argc,argv))
    {
      for (i=0; i<DIM; i++)
        if (np->cs[i] <= 0.0)
          ret = NP_NOT_ACTIVE;
    }
    else
    if (dval[0] <= 0.0)
    {
      PrintErrorMessage('E',"NPStochFieldInit","cell size must be positiv");
      ret = NP_NOT_ACTIVE;
    }
    else
      for (i=0; i<DIM; i++)
        np->cs[i] = dval[0];
  else
    for (i=0; i<DIM; i++)
      if (dval[i] > 0.0)
        np->cs[i] = dval[i];
      else
      {
        PrintErrorMessage('E',"NPStochFieldInit","cell size must be positiv");
        ret = NP_NOT_ACTIVE;
      }

  if (ReadArgvOption("e",argc,argv))
    if(ReadArgvOption("b",argc,argv))
    {
      PrintErrorMessage('E',"NPStochFieldInit","bell-shaped exclusive or exponential autocor.");
      ret = NP_NOT_ACTIVE;
    }
    else
      np->actype = FIELD_EXPONENTIAL;
  else if (ReadArgvOption("b",argc,argv))
    np->actype = FIELD_GAUSSIAN;
  else if ((np->actype != FIELD_EXPONENTIAL) && (np->actype != FIELD_GAUSSIAN))
    ret = NP_NOT_ACTIVE;

  if (ReadArgvINT("i",ival,argc,argv))
  {
    if (np->initial <= 0)
      ret = NP_NOT_ACTIVE;
  }
  else
  {
    if (ival[0] < 0)
    {
      PrintErrorMessage('E',"NPStochFieldInit","positive initial value");
      ret = NP_NOT_ACTIVE;
    }
    if (ival[0] == 0)
      np->initial = (INT) time (NULL);
    else
      np->initial = ival[0];
  }

  if (ReadArgvOption("lin",argc,argv))
    if(ReadArgvOption("const",argc,argv))
    {
      PrintErrorMessage('E',"NPStochFieldInit","linear interpolation exclusive or constant value");
      ret = NP_NOT_ACTIVE;
    }
    else
      np->inttype = TRUE;
  else if (ReadArgvOption("const",argc,argv))
    np->inttype = FALSE;
  else if ((np->inttype != TRUE) && (np->inttype != FALSE))
    ret = NP_NOT_ACTIVE;

  if (sopt == 1)
  {
    if (np->Fld != NULL)
      PutFreelistMemory(theHeap, np->Fld, np->FldSize);

    NumOfDouble = 1;
    for (i=0; i<DIM; i++)
      NumOfDouble *= np->size[i];
    np->FldSize = NumOfDouble*sizeof(DOUBLE);
    if ((theField=GetFreelistMemory(theHeap,np->FldSize))==NULL)
    {
      PrintErrorMessage('E',"NPStochFieldInit","not enough memory");
      ret = NP_NOT_ACTIVE;
    }
    else
      np->Fld=theField;
  }

  if (ret == NP_ACTIVE)
    if (Field_genStochField(np))
    {
      PrintErrorMessage('E',"NPStochFieldInit","Cannot initialize the stoch. field");
      ret = NP_NOT_ACTIVE;
    }

  return(ret);
}

static INT NPStochFieldDisplay(NP_BASE *theNP)
{
  NP_STOCH_FIELD *np;

  np = (NP_STOCH_FIELD *) theNP;

#ifdef __THREEDIM__
  UserWriteF(DISPLAY_NP_FORMAT_SIII,"Size",np->size[0],np->size[1]
             ,np->size[2]);
#else
  UserWriteF(DISPLAY_NP_FORMAT_SII,"Size",np->size[0],np->size[1]);
#endif

  UserWriteF(DISPLAY_NP_FORMAT_SF,"Mean value",np->mean);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Variance",np->var);
#ifdef __THREEDIM__
  UserWriteF(DISPLAY_NP_FORMAT_SFFF,"Cor. lengths",np->cor[0],np->cor[1]
             ,np->cor[2]);
  UserWriteF(DISPLAY_NP_FORMAT_SFFF,"Cell size",np->cs[0],np->cs[1]
             ,np->cs[2]);
#else
  UserWriteF(DISPLAY_NP_FORMAT_SFF,"Cor. lengths",np->cor[0],np->cor[1]);
  UserWriteF(DISPLAY_NP_FORMAT_SFF,"Cell size",np->cs[0],np->cs[1]);
#endif
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Nugget",np->nugget);
  if (np->actype==FIELD_EXPONENTIAL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Autocorrelation","exponential");
  else if (np->actype==FIELD_GAUSSIAN)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Autocorrelation","gaussian");

  if (np->initial > 0)
    UserWriteF(DISPLAY_NP_FORMAT_SI,"Random initial",np->initial);
  else
    UserWriteF(DISPLAY_NP_FORMAT_S,"Random initial");

  if (np->inttype==TRUE)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Interpolation","linear in each dir");
  else if (np->inttype==FALSE)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Interpolation","constant on cells");

  return(0);
}

INT Field_RandomValues (NP_FIELD *theField, DOUBLE *Pos, DOUBLE *out)
{
  NP_STOCH_FIELD *np;
  INT i, node[DIM];
  DOUBLE cornerValues[8], alpha[DIM], beta;

  np = (NP_STOCH_FIELD*)theField;
  if (np->Fld == NULL) return(1);

  for (i=0; i<DIM; i++)
  {
    alpha[i] = Pos[i] * np->cor[i] / np->cs[i];
    node[i] = ((INT) alpha[i]) % np->size[i];
    if (node[i] < 0)
    {
      node[i] += np->size[i];
      alpha[i] = ((INT) alpha[i]) - alpha[i];
    }
    else
      alpha[i] -= (INT) alpha[i];
  }
#ifdef __THREEDIM__
  switch(np->inttype)
  {
  case FALSE :
    out[0] = (REALPART(np->size,np->Fld,node[0],node[1],node[2]) - np->mean) / sqrt(np->var);
    break;
  case TRUE :
  {
    INT extraNode[DIM];

    for (i=0; i<DIM; i++)
      extraNode[i] = (node[i]+1) % np->size[i];
    cornerValues[0] =  REALPART(np->size,np->Fld,node[0],node[1],node[2]);
    cornerValues[1] =  REALPART(np->size,np->Fld,extraNode[0],node[1],node[2]);
    cornerValues[2] =  REALPART(np->size,np->Fld,node[0],extraNode[1],node[2]);
    cornerValues[3] =  REALPART(np->size,np->Fld,extraNode[0],extraNode[1],node[2]);
    cornerValues[4] =  REALPART(np->size,np->Fld,node[0],node[1],extraNode[2]);
    cornerValues[5] =  REALPART(np->size,np->Fld,extraNode[0],node[1],extraNode[2]);
    cornerValues[6] =  REALPART(np->size,np->Fld,node[0],extraNode[1],extraNode[2]);
    cornerValues[7] =  REALPART(np->size,np->Fld,extraNode[0],extraNode[1],extraNode[2]);
    beta = 1.0 - alpha[2];
    for (i=0; i<4; i++)
      cornerValues[i] = beta*cornerValues[i] + alpha[2]*cornerValues[i+4];
    beta = 1.0 - alpha[1];
    for (i=0; i<2; i++)
      cornerValues[i] = beta*cornerValues[i] + alpha[1]*cornerValues[i+2];
    beta = 1.0 - alpha[0];
    out[0] = beta*cornerValues[0] + alpha[0]*cornerValues[1];
  }
  break;
  default :
    return(1);
  }
#else
  switch(np->inttype)
  {
  case FALSE :
    out[0] = (REALPART(np->size,np->Fld,node[0],node[1]) - np->mean) / sqrt(np->var);
    break;
  case TRUE :
  {
    INT extraNode[DIM];

    for (i=0; i<DIM; i++)
      extraNode[i] = (node[i]+1) % np->size[i];
    cornerValues[0] =  REALPART(np->size,np->Fld,node[0],node[1]);
    cornerValues[1] =  REALPART(np->size,np->Fld,extraNode[0],node[1]);
    cornerValues[2] =  REALPART(np->size,np->Fld,node[0],extraNode[1]);
    cornerValues[3] =  REALPART(np->size,np->Fld,extraNode[0],extraNode[1]);
    beta = 1.0 - alpha[1];
    for (i=0; i<2; i++)
      cornerValues[i] = beta*cornerValues[i] + alpha[1]*cornerValues[i+2];
    beta = 1.0 - alpha[0];
    out[0] = beta*cornerValues[0] + alpha[0]*cornerValues[1];
  }
  break;
  default :
    return(1);
  }
#endif

  return(0);
}

static INT StochFieldConstruct (NP_BASE *theNP)
{
  NP_FIELD *theField;
  NP_STOCH_FIELD *np;
  INT i;

  theNP->Init     = NPStochFieldInit;
  theNP->Display  = NPStochFieldDisplay;
  theNP->Execute  = NULL;

  /* field entry */
  theField = (NP_FIELD *) theNP;
  theField->Evaluate = Field_RandomValues;

  /* stochastic parameters */
  np = (NP_STOCH_FIELD *) theNP;
  for (i=0; i<DIM; i++)
  {
    np->size[i] =  0;
    np->cor[i]  = -1.0;
    np->cs[i]   =  1.0;
  }
  np->mean = 0.0;
  np->var = -1.0;
  np->nugget = -0.0;
  np->actype = 0;
  np->initial = -1;
  np->inttype = 0;
  np->Fld = NULL;

  return(0);
}

/****************************************************************************/
/*D
   NP_GET_FIELD - type definition for stochastic fields evaluation

   DESCRIPTION:
   This numproc type is used for the evaluation of
   stochastic fields. It must be initialized by

   'INT NPGetFieldInit (NP_BASE *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPGetFieldDisplay (NP_BASE *theNP);'
   'INT NPGetFldExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   npinit $F <npfield> [$M <mean>] [$V <var>] [$C <cor> [<cor2> <cor3>]]
       [$NOR|$lLOGNOR]

   .ve

   SEE ALSO:
   num_proc
   D*/
/****************************************************************************/

static INT NPGetFieldInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_GET_FIELD *np;
  char field[VALUELEN];
  DOUBLE dval[3];
  INT i, ret;

  ret = NP_ACTIVE;

  np = (NP_GET_FIELD *) theNP;

  for (i=1; i<argc; i++)
    if (argv[i][0]=='F')
    {
      if (sscanf(argv[i],"F %s",field)!=1)
      {
        PrintErrorMessage('E',"NPGetFieldInit","stoch data np missing");
        ret = NP_NOT_ACTIVE;
      }
      else
        np->FldNp = (NP_FIELD *)GetNumProcByName(NP_MG(theNP),field,FIELD_CLASS_NAME);
    }

  if (ReadArgvDOUBLE("M",dval,argc,argv))
  {
    if (np->mean == 0.0)
      ret = NP_NOT_ACTIVE;
  }
  else
  if (dval[0] == 0.0)
  {
    PrintErrorMessage('E',"NPGetFieldInit","vanishing mean");
    ret= NP_NOT_ACTIVE;
  }
  else
    np->mean = dval[0];

  if (ReadArgvDOUBLE("V",dval,argc,argv))
  {
    if (np->var < 0.0)
      ret = NP_NOT_ACTIVE;
  }
  else
  if (dval[0] < 0.0)
  {
    PrintErrorMessage('E',"NPGetFieldInit","negative variance");
    ret = NP_NOT_ACTIVE;
  }
  else
    np->var = dval[0];

  if (ReadArgvPosition("C",argc,argv,dval))
    if (ReadArgvDOUBLE("C",dval,argc,argv))
    {
      for (i=0; i<DIM; i++)
        if (np->cor[i] <= 0.0)
          ret = NP_NOT_ACTIVE;
    }
    else
    if (dval[0] <= 0.0)
    {
      PrintErrorMessage('E',"NPGetFieldInit","correlation must be positiv");
      ret = NP_NOT_ACTIVE;
    }
    else
      for (i=0; i<DIM; i++)
        np->cor[i] = dval[0];
  else
    for (i=0; i<DIM; i++)
      if (dval[i] > 0.0)
        np->cor[i] = dval[i];
      else
      {
        PrintErrorMessage('E',"NPGetFieldInit","correlation must be positiv");
        ret = NP_NOT_ACTIVE;
      }

  if (ReadArgvOption("NOR",argc,argv))
    if (ReadArgvOption("LOGNOR",argc,argv))
    {
      PrintErrorMessage('E',"NPGetFieldInit","normal- and lognormaldistributed are exclusive");
      ret= NP_NOT_ACTIVE;
    }
    else
      np->dtype = FIELD_NORMDIST;
  else if (ReadArgvOption("LOGNOR",argc,argv))
    np->dtype = FIELD_LOGNORM;
  else if ((np->dtype != FIELD_NORMDIST) && (np->dtype != FIELD_LOGNORM))
    ret = NP_NOT_ACTIVE;

  if (np->FldNp == NULL)
    ret = NP_NOT_ACTIVE;

  return(ret);
}

static INT NPGetFieldDisplay(NP_BASE *theNP)
{
  NP_GET_FIELD *np;

  np = (NP_GET_FIELD *) theNP;
  if (np->FldNp != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"field",np->FldNp->base.v.name);
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"field","---");
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Mean value",np->mean);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Variance",np->var);
#ifdef __THREEDIM__
  UserWriteF(DISPLAY_NP_FORMAT_SFFF,"Cor. lengths",np->cor[0],np->cor[1],np->cor[2]);
#else
  UserWriteF(DISPLAY_NP_FORMAT_SFF,"Cor. lengths",np->cor[0],np->cor[1]);
#endif
  if (np->dtype==FIELD_NORMDIST) UserWriteF(DISPLAY_NP_FORMAT_SS,"Distribution","normal distributed");
  else if (np->dtype==FIELD_LOGNORM) UserWriteF(DISPLAY_NP_FORMAT_SS,"Distribution","lognormal");

  return(0);
}

INT Field_GetFieldAtPoint (NP_FIELD *theField, DOUBLE *Pos, DOUBLE *out)
{
  NP_GET_FIELD *np;
  NP_FIELD *npsd;
  INT i;
  DOUBLE zeta, lambda, transX[DIM], value[1];

  np = (NP_GET_FIELD *) theField;
  if ((npsd = np->FldNp) == NULL) return(1);

  for (i=0; i<DIM; i++)
    transX[i] = Pos[i]/np->cor[i];

  if ((*npsd->Evaluate)(npsd, transX, value)) return(1);

  switch (np->dtype)
  {
  case FIELD_NORMDIST :
    out[0] = np->mean + sqrt(np->var) * value[0];
    break;
  case FIELD_LOGNORM :
    zeta =  sqrt(log(np->var/(np->mean * np->mean) + 1.0));
    lambda = log(np->mean) - (zeta * zeta) / 2.0;
    out[0] = exp(lambda + zeta * value[0]);
    break;
  default :
    return(1);
  }

  return(0);
}

static INT GetFieldConstruct (NP_BASE *theNP)
{
  NP_FIELD *theField;
  NP_GET_FIELD *np;
  INT i;

  theNP->Init             = NPGetFieldInit;
  theNP->Display          = NPGetFieldDisplay;
  theNP->Execute          = NULL;

  /* field evaluation */
  theField = (NP_FIELD *) theNP;
  theField->Evaluate = Field_GetFieldAtPoint;

  np = (NP_GET_FIELD *) theNP;
  for (i=0; i<DIM; i++)
    np->cor[i] = -1.0;
  np->mean = 0.0;
  np->var = -1.0;
  np->dtype = 0;
  np->FldNp = NULL;

  return(0);
}

/****************************************************************************/
/*D
   NP_ANISO_FLD - type definition for anisotropic stochastic fields

   DESCRIPTION:
   This numproc type is used for the evaluation of
   stochastic fields. It must be initialized by

   'INT NPGetFieldInit (NP_BASE *theNP, INT argc , char **argv);'

   This routine returns 'EXECUTABLE' if the initizialization is complete
   and  'ACTIVE' else.
   The data they can be displayed and the num proc can be executed by

   'INT NPGetFieldDisplay (NP_BASE *theNP);'
   'INT NPGetFldExecute (NP_BASE *theNP, INT argc , char **argv);'

   .vb
   npinit $F <npfield> [$M <mean>] [$V <var>] [$E <eang1> <eang2> <eang3>]
       [$A <angle>] [$C <cor> [<cor2> <cor3>]] [$NOR|$lLOGNOR] [$I|$NI]

   .ve

   SEE ALSO:
   num_proc, NP_GET_FIELD
   D*/
/****************************************************************************/

static INT NPanisoFldInit(NP_BASE *theNP, INT argc , char **argv)
{
  NP_ANISO_FIELD *np;
  DOUBLE dval[DIM];
  INT ret;
#ifdef __THREEDIM__
  INT i;
#endif

  ret = NPGetFieldInit (theNP, argc, argv);

  np = (NP_ANISO_FIELD *) theNP;

#ifdef __THREEDIM__
  if (ReadArgvPosition("E",argc,argv,dval))
  {
    for (i=0; i<DIM; i++)
      if ((np->euler[i] < -180.0) || (np->euler[i] > 360.0))
        ret = NP_NOT_ACTIVE;
  }
  else
    for (i=0; i<DIM; i++)
      if ((dval[i] < -180.0) || (dval[i] > 360.0))
      {
        PrintErrorMessage('E',"NPGetFieldInit","Euler angle in -180..360");
        ret = NP_NOT_ACTIVE;
      }
      else
        np->euler[i] = dval[i];
#else
  if (ReadArgvDOUBLE("A",dval, argc, argv))
  {
    if ((np->angle < -180.0) || (np->angle > 360.0))
      ret = NP_NOT_ACTIVE;
  }
  else
  if ((dval[0] < -180.0) || (dval[0] > 360.0))
  {
    PrintErrorMessage('E',"NPGetFieldInit","Angle should be in -180..360");
    ret = NP_NOT_ACTIVE;
  }
  else
    np->angle = dval[0];
#endif

  return(ret);
}


static INT NPanisoFldDisplay(NP_BASE *theNP)
{
  NP_ANISO_FIELD *np;

  NPGetFieldDisplay(theNP);

  np = (NP_ANISO_FIELD *) theNP;

#ifdef __THREEDIM__
  UserWriteF(DISPLAY_NP_FORMAT_SFFF,"Euler angle",
             np->euler[0],np->euler[1],np->euler[2]);
#else
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Rotation",np->angle);
#endif

  return(0);
}

INT Field_RotateAndGetField (NP_FIELD *theField, DOUBLE *Pos, DOUBLE *out)
{
  NP_ANISO_FIELD  *np;
  DOUBLE transX[DIM];
#ifdef __THREEDIM__
  DOUBLE sinus[DIM], cosinus[DIM];
  INT i;
#else
  DOUBLE sinus, cosinus;
#endif

  np      = (NP_ANISO_FIELD *) theField;

#ifdef __TWODIM__
  sinus = sin( - np->angle*PI/180.0);
  cosinus = cos(np->angle*PI/180.0);
  transX[0] = Pos[0] * cosinus - Pos[1] * sinus;
  transX[1] = Pos[0] * sinus + Pos[1] * cosinus;
#else
  for (i=0; i<DIM; i++)
  {
    sinus[i] = sin(np->euler[i]*PI/180.0);
    cosinus[i] = cos(np->euler[i]*PI/180.0);
  }
  transX[0] = (cosinus[2]*cosinus[0] - cosinus[1]*sinus[0]*sinus[2]) * Pos[0]
              - (sinus[2]*cosinus[0] - cosinus[1]*sinus[0]*cosinus[2]) * Pos[1]
              + sinus[1]*sinus[0] * Pos[2];

  transX[1] = (cosinus[2]*sinus[0] + cosinus[1]*cosinus[0]*sinus[2]) * Pos[0]
              - (sinus[2]*sinus[0] + cosinus[1]*cosinus[0]*cosinus[2]) * Pos[1]
              - sinus[1]*cosinus[0] * Pos[2];

  transX[2] = sinus[1]*sinus[2] * Pos[0]
              + sinus[1]*cosinus[2] * Pos[1]
              + cosinus[1] * Pos[2];
#endif

  return(Field_GetFieldAtPoint(theField, transX, out));
}


static INT GetAnisoFieldConstruct       (NP_BASE *theNP)
{
  NP_FIELD *theField;
  NP_GET_FIELD *fldnp;
  NP_ANISO_FIELD *anisonp;
  INT i;

  theNP->Init = NPanisoFldInit;
  theNP->Display = NPanisoFldDisplay;
  theNP->Execute = NULL;

  /* field entry */
  theField = (NP_FIELD *) theNP;
  theField->Evaluate = Field_RotateAndGetField;

  fldnp = (NP_GET_FIELD *) theNP;
  for (i=0; i<DIM; i++)
    fldnp->cor[i] = -1.0;
  fldnp->mean = 0.0;
  fldnp->var = -1.0;
  fldnp->dtype = 0;
  fldnp->FldNp = NULL;

  anisonp =(NP_ANISO_FIELD *) theNP;
#ifdef __THREEDIM__
  for (i=0; i<DIM; i++)
    anisonp->euler[i] = 0.0;
#else
  anisonp->angle = 0.0;
#endif

  return(0);
}

INT InitStochField (void)
{
  if (CreateClass(FIELD_CLASS_NAME ".stoch",sizeof(NP_STOCH_FIELD),StochFieldConstruct)) return(__LINE__);
  if (CreateClass(FIELD_CLASS_NAME ".scale",sizeof(NP_GET_FIELD),GetFieldConstruct)) return(__LINE__);
  if (CreateClass(FIELD_CLASS_NAME ".rot",sizeof(NP_ANISO_FIELD),GetAnisoFieldConstruct)) return(__LINE__);

  return(0);
}
