// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      f2c.c                                                         */
/*                                                                          */
/* Purpose:   init & exit functions for the ddd fortran interface           */
/*                                                                          */
/* Author:    Jens Boenisch                                                 */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: boenisch@rus.uni-stuttgart.de                       */
/*                                                                          */
/* History:   95/12/07 jb  begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifdef F_FRONTEND



/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*                                                                          */
/****************************************************************************/

#include <strings.h>
#include <sys/types.h>
#include <stdio.h>

#include "dddi.h"


#define FORTRAN_MAIN   F77SYM(fortran_main,FORTRAN_MAIN)


/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/****************************************************************************/

/*
        with F_FRONTEND, ddd has the main program, and calls the F77
        main program as a Fortran subroutine.
 */

int main (int argc, char **argv)
{
  DDD_Init (&argc, &argv);

  /* call the fortran main routine */
  FORTRAN_MAIN ();

  DDD_Exit ();
}


/****************************************************************************/

void DDD_SetConfig (int *me1, int *master1, int *procs1)
{
  *me1 = me;
  *master1 = master;
  *procs1 = procs;
}



/****************************************************************************/

#endif
