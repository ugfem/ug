// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*									                                                                        */
/* File:	field.h						                                                        */
/*									                                                                        */
/* Purpose:	definition of the field num proc type		                                */
/*									                                                                        */
/* Author:	Carsten Schwarz						                                                */
/*		Institut fuer Hydromechanik und Wasserwirtschaft	                        */
/*		ETH Hoenggerberg					                                                        */
/*		8093 Zuerich						                                                        */
/*		email: schwarz@ihw.baum.ethz.ch			                                        */
/*									                                                                        */
/* History:	February 1997						                                                */
/*									                                                                        */
/* Remarks:                                                                                                                     */
/*									                                                                        */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*									                                                                        */
/* auto include mechanism and other include files			                        */
/*									                                                                        */
/****************************************************************************/

#ifndef __NPFIELD__
#define __NPFIELD__

#include "np.h"

/****************************************************************************/
/*									                                                                        */
/* defines in the following order					                                        */
/*									                                                                        */
/*	compile time constants defining static data size (i.e. arrays)	        */
/*	other constants							                                                                */
/*	macros								                                                                */
/*									                                                                        */
/****************************************************************************/

#define FIELD_CLASS_NAME        "field"

#define FIELD_EXPONENTIAL       1
#define FIELD_GAUSSIAN          2

#define FIELD_LOGNORM           1
#define FIELD_NORMDIST          2

/****************************************************************************/
/*									                                                                        */
/* definition of exported data structures				                                */
/*									                                                                        */
/****************************************************************************/

struct np_field
{
  NP_BASE base;
  INT (*Evaluate)(struct np_field *np, DOUBLE *Pos, DOUBLE *out);
};

typedef struct np_field NP_FIELD;

struct np_stoch_field
{
  NP_FIELD field;

  /* configuration */
  INT size[DIM];
  DOUBLE mean;
  DOUBLE var;
  DOUBLE cor[DIM];
  DOUBLE cs[DIM];
  DOUBLE nugget;
  INT actype;
  INT inttype;
  INT initial;
  DOUBLE *Fld;
  MEM FldSize;                          /* current allocation size of Fld in byte		*/
};


struct np_get_fld
{
  NP_FIELD field;

  /* configuration */
  DOUBLE mean;
  DOUBLE var;
  DOUBLE cor[DIM];
  INT dtype;
  NP_FIELD *FldNp;
};

struct np_aniso_fld
{
  struct np_get_fld field;

  /* configuration */
#ifdef __THREEDIM__
  DOUBLE euler[3];
#else
  DOUBLE angle;
#endif
};

typedef struct np_stoch_field NP_STOCH_FIELD;
typedef struct np_get_fld NP_GET_FIELD;
typedef struct np_aniso_fld NP_ANISO_FIELD;

/****************************************************************************/
/*									                                                                        */
/* definition of exported functions					                                        */
/*									                                                                        */
/****************************************************************************/

INT Field_genStochField(NP_STOCH_FIELD *np);
INT Field_RandomValues (NP_FIELD *theField, DOUBLE *Pos, DOUBLE *out);
INT Field_GetFieldAtPoint (NP_FIELD *theField, DOUBLE *Pos, DOUBLE *out);
INT Field_RotateAndGetField (NP_FIELD *theField, DOUBLE *Pos, DOUBLE *out);

INT InitStochField (void);

#endif
