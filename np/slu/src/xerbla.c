// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*******************************************************/
/* This file was modified according to be integrated   */
/* in the ug-software-package. The changes are due to  */
/*                                                     */
/* Klaus Johannsen,                                    */
/* University of Heidelberg,                           */
/* Im Neuenheimer Feld 368,                            */
/* Germany.                                            */
/*******************************************************/

/* Subroutine */ int xerbla_(char *srname, int *info)
{
  /*  -- LAPACK auxiliary routine (version 2.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         September 30, 1994


      Purpose
      =======

      XERBLA  is an error handler for the LAPACK routines.
      It is called by an LAPACK routine if an input parameter has an
      invalid value.  A message is printed and execution stops.

      Installers may consider modifying the STOP statement in order to
      call system-specific exception-handling facilities.

      Arguments
      =========

      SRNAME  (input) CHARACTER*6
              The name of the routine which called XERBLA.

      INFO    (input) INT
              The position of the invalid parameter in the parameter list

              of the calling routine.

     =====================================================================
   */

  printf("** On entry to %6s, parameter number %2d had an illegal value\n",
         srname, *info);

  /*     End of XERBLA */

  return 0;
} /* xerbla_ */
