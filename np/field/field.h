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

#define FIELD_CLASS_NAME                        "field"

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

/****************************************************************************/
/*									                                                                        */
/* definition of exported functions					                                        */
/*									                                                                        */
/****************************************************************************/

INT InitStochField ();

#endif
