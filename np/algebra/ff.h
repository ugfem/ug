// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  tff.h															*/
/*																			*/
/* Purpose:   (tangential) frequency filtering decompostion routines        */
/*																			*/
/* Author:	  Christian Wrobel                                                                              */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
/*																			*/
/* History:   14.11.95 begin, ug version 3.1								*/
/*																			*/
/* Remarks:   FF is used as the abbreviation for "frequency filtering"		*/
/*                        TFF is used as the abbreviation for the "tangential frequency */
/*			       filtering" method due to Christian Wagner, 1995			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __TFF__
#define __TFF__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT TFFDecomp( DOUBLE wavenr,
               DOUBLE wavenr3D,
               const BLOCKVECTOR *bv,
               const BV_DESC *bvd,
               const BV_DESC_FORMAT *bvdf,
               INT tv_comp,
               GRID *grid );

INT FFDecomp( DOUBLE wavenr,
              DOUBLE wavenr3D,
              const BLOCKVECTOR *bv,
              const BV_DESC *bvd,
              const BV_DESC_FORMAT *bvdf,
              INT tv1_comp,
              INT tv2_comp,
              GRID *grid );


INT TFFPrepareSolver( GRID *grid,
                      INT K_comp,
                      INT u_comp,
                      INT f_comp,
                      const BV_DESC_FORMAT *bvdf );

INT TFFSolve( const BLOCKVECTOR *bv,
              const BV_DESC *bvd,
              const BV_DESC_FORMAT *bvdf,
              INT K_comp,
              INT u_comp,
              INT f_comp,
              INT cor_comp,
              INT FF_comp,
              INT LU_comp,
              INT tv_comp,
              INT aux_comp,
              INT aux3D_comp,
              INT FF3D_comp,
              DOUBLE meshwidth,
              DOUBLE eps,
              GRID *grid );
#endif
