// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  connectuggrapeOFF.c		                                                                        */
/*																			*/
/* Purpose:   Provides procedure dummies if GRAPE is not included           */
/*																			*/
/* Author:	  Klaus Johannsen				                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/*			  Monika Wierse/Martin Metscher									*/
/*            Universitaet Freiburg											*/
/*            Institut fuer Angewandte Mathematik							*/
/*            Hermann--Herder--Str. 10										*/
/*            D-79104 Freiburg												*/
/*																			*/
/* History:   27.04.96 begin, ug version 3.1								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* ug includes */
#include "defs.h"
#include "gm.h"
#include "ugdevices.h"
#include "evm.h"
#include "general.h"
#include "connectuggrape.h"


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);



int CallGrape (MULTIGRID *theMG)
{
  UserWrite("Grape library not included!\nIf Grape is available set GRAPE=ON in ug.conf and recompile.\n");
  return(0);
}

INT InitGrape (void)
{
  return (0);
}

#ifndef __AIX__
void usleep (unsigned long time)
{
  return;
}
#endif
