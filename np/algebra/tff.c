// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  tff.c                                                                                                                 */
/*																			*/
/* Purpose:   tangential frequency filtering decompostion routines          */
/*                                                                                                                                              */
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
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <assert.h>
#include <math.h>
#include <time.h>

#include "switch.h"             /* for  __TWODIM__ and  __THREEDIM__ */
#include "compiler.h"
#include "gm.h"          /* for data structure               */
#include "ugstruct.h"    /* for GetStringValue               */
#include "misc.h"        /* for MIN, MAX, PI, ...            */
#include "devices.h"     /* for UserWrite, PrintErrorMessage */
#include "commands.h"    /* for GetCurrentMultigrid          */
#include "debug.h"

#include "num.h"
#include "ugblas.h"

#include "ff_gen.h"
#include "tff.h"

#ifdef T
#undef T
#endif

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
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);

/* value below them a division is refused calculating a testvector */
DOUBLE TFFsmallTV = 1e-3;

/* ratio for a jump to be detected */
DOUBLE TFFmuchBigger = 100.0;

/* value below them a number is considered as 0.0 */
DOUBLE TFFEPS = 1e-16;

/* value below them an approximation error is considered as ok */
DOUBLE TFFaccuracy = 1e-10;

static DOUBLE glob_h;   /* used if THETA_ANA is defined */

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

#ifdef __BLOCK_VECTOR_DESC__

void TFFConstructTestvector( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr, DOUBLE wavenr_3D );

INT TFFCalculateTheta( const BLOCKVECTOR *bv_dest, const BLOCKVECTOR *bv_source, const BV_DESC *bvd_dest, const BV_DESC *bvd_source, const BV_DESC_FORMAT *bvdf, INT Theta, INT Tinv, INT L, INT tv_comp, INT aux_comp, INT auxsub_comp, INT Lsub_comp );

INT TFFUpdateDiagBlock( const BLOCKVECTOR *bv_dest, const BV_DESC *bvd_dest, const BV_DESC *bvd_source, const BV_DESC_FORMAT *bvdf, INT T, INT DL, INT Theta, GRID *grid );

INT TFFMultWithM( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT y_comp, INT T_comp, INT L_comp, INT Tinv_comp, INT x_comp, INT aux_comp, INT auxsub_comp, INT Lsub_comp );

INT InitTFF( void );

/****************************************************************************/
/*																			*/
/* auxiliary routines for frequency filtering								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   TFFConstructTestvector - construct a sine-shaped testvector (2D/3D)

   SYNOPSIS:
   INT TFFConstructTestvector( const BLOCKVECTOR *bv, INT tv_comp,
                                                           DOUBLE wavenr, DOUBLE wavenr_3D );

   PARAMETERS:
   .  bv - root of the blockvector tree in 2D/3D (plane/cube)
   .  tv_comp - component for the testvector in the vector-data
   .  wavenr - number of sine half-oszillations along one gridline (2D (sub)problems)
   .  wavenr3D - number of sine half-oszillations along the planes (only 3D problems)

   DESCRIPTION:
   Calculate the entries of a testvector to a given frequency with respect
   to the meshwidth (calculated as the number of vectors within 1 gridline).
   2D: let h be the meshwidth, (ih,jh) a grid point,
   then tv(i,j) = sin( i*h*wavenr*pi )
   3D: let h be the meshwidth, (ih,jh,kh) a grid point,
   then tv(i,j,k) = sin( i*h*wavenr*pi ) * sin( j*h*wavenr*pi )

   'wavenr3D' is used only in the 3D case.

   RESTRICTIONS:
   The grid must be ordered by a blockvector structure as created by
   'CreateBVStripe2D'/'3D'. Especially the leaf-blockvectors must have
   the right 'BVNUMBEROFVECTORS' and in 3D the sons of 'bv' must be numbered
   consecutively to determine the right meshwidth.

   The spcial dimension is determined staticly by '__TWODIM'/'__THREEDIM'.

   The grid should be regular (i.e. each line/plane contains the same number of
   vectors).

   WARNING:
   In 3D some lines may consist entirely of testvector entries 0.0. To
   circumvent this problem use 'TFFConstructTestvector_loc'.

   SEE ALSO:
   TFFConstructTestvector_loc, CreateBVStripe2D, CreateBVStripe3D

   RETURN VALUE:
   .n   void
   D*/
/*************************************************************************/

void TFFConstructTestvector( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr, DOUBLE wavenr_3D )
{
  register DOUBLE hkpi, pos, plane_pos, tensor;
  register VECTOR *v, *end_v;
  register BLOCKVECTOR *bv_i, *bv_end, *bv_glob_end;
  INT length, plane_length;
  DOUBLE plane_hkpi;

  bv_glob_end = BVDOWNBVEND(bv);
  bv = BVDOWNBV(bv);
  for ( ; bv != bv_glob_end; bv = BVSUCC(bv) )       /* over all lines resp. planes */
  {
#ifdef __TWODIM__
    assert( BVDOWNTYPE(bv) == BVDOWNTYPEVECTOR );
    length = BVNUMBEROFVECTORS(bv) + 1;
    hkpi = pos = ( PI * wavenr ) / (DOUBLE)length;
    end_v = BVENDVECTOR( bv );
    BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )                /* over all points in the line */
    {
      VVALUE( v, tv_comp ) = sin( pos );
      pos += hkpi;
    }
#else
    assert( BVDOWNTYPE(bv) == BVDOWNTYPEBV );
    bv_i = BVDOWNBV(bv);
    bv_end = BVDOWNBVEND(bv);

    plane_length = BVNUMBER(BVDOWNBVLAST(bv)) - BVNUMBER(bv_i) + 2;
    plane_hkpi = plane_pos = (PI * wavenr_3D) / (double)plane_length;
    tensor = sin ( plane_pos );

    for ( ; bv_i != bv_end; bv_i = BVSUCC(bv_i) )               /* over all lines */
    {
      length = BVNUMBEROFVECTORS(bv_i) + 1;
      hkpi = pos = ( PI * wavenr ) / (DOUBLE)length;
      end_v = BVENDVECTOR( bv_i );
      BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv_i), end_v )                    /* over all points in the line */
      {
        VVALUE( v, tv_comp ) = tensor * sin( pos );
        pos += hkpi;
      }
      plane_pos += plane_hkpi;
      tensor = sin( plane_pos );
    }
#endif
  }
}

/****************************************************************************/
/*D
   TFFConstructTestvector_loc - construct a sine-shaped testvector (2D/3D)

   SYNOPSIS:
   INT TFFConstructTestvector_loc( const BLOCKVECTOR *bv, INT tv_comp,
                                                                   DOUBLE wavenr, DOUBLE wavenr_3D );

   PARAMETERS:
   .  bv - root of the blockvector tree in 2D/3D (plane/cube)
   .  tv_comp - component for the testvector in the vector-data
   .  wavenr - number of sine half-oszillations along one gridline (2D (sub)problems)
   .  wavenr3D - number of sine half-oszillations along the planes (only 3D problems)

   DESCRIPTION:
   Calculate the entries of a testvector to a given frequency with respect
   to the meshwidth (calculated as the number of vectors within 1 gridline).
   The distinction whether 'bv' represents a plane (2D) or a cube (3D) is made
   dynamically by the blockvector structure itself: if 'bv' has only one
   further level of blockvectors (representing grid lines!) it is
   considered as a plane, otherwise as a cube.
   plane: let h be the meshwidth, (ih,jh) a grid point,
   then tv(i,j) = sin( i*h*wavenr*pi )
   cube: let h be the meshwidth, (ih,jh,kh) a grid point,
   then tv(i,j,k) = sin( i*h*wavenr*pi ) * sin( j*h*wavenr*pi )

   In 3D some lines may consist entirely of testvector entries 0.0 only
   in the case of the global decomposition; in the case of decomposing
   local problems no 0-lines can occur!
   'wavenr3D' is used only in the 3D case.

   WARNING:
   The local subproblems do not use the global testvector, but a local one.

   RESTRICTIONS:
   The grid must be ordered by a blockvector structure as created by
   'CreateBVStripe2D'/'3D'.

   The grid should be regular (i.e. each line/plane contains the same number of
   vectors).

   SEE ALSO:
   TFFConstructTestvector, CreateBVStripe2D, CreateBVStripe3D

   RETURN VALUE:
   .n   void
   D*/
/*************************************************************************/

void TFFConstructTestvector_loc( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr, DOUBLE wavenr_3D )
{
  register DOUBLE hkpi, pos, plane_pos, tensor;
  register VECTOR *v, *end_v;
  register BLOCKVECTOR *bv_i, *bv_end;
  INT length, plane_length;
  DOUBLE plane_hkpi;

  if ( BV_IS_LEAF_BV( bv ) )
  /* 2D block */
  {
    length = BVNUMBEROFVECTORS(bv) + 1;
    hkpi = pos = ( PI * wavenr ) / (DOUBLE)length;
    end_v = BVENDVECTOR( bv );
    BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )                /* over all points in the line */
    {
      /*printf("pos %g  sin %g\n", pos, sin(pos) );*/
      VVALUE( v, tv_comp ) = sin( pos );
      pos += hkpi;
    }
  }
  else
  /* 3D block */
  {
    assert( BVDOWNTYPE(bv) == BVDOWNTYPEBV );
    bv_i = BVDOWNBV(bv);
    bv_end = BVDOWNBVEND(bv);

    plane_length = BVNUMBER(BVDOWNBVLAST(bv)) - BVNUMBER(bv_i) + 2;
    plane_hkpi = plane_pos = (PI * wavenr_3D) / (double)plane_length;
    tensor = sin ( plane_pos );

    for ( ; bv_i != bv_end; bv_i = BVSUCC(bv_i) )               /* over all lines */
    {
      length = BVNUMBEROFVECTORS(bv_i) + 1;
      hkpi = pos = ( PI * wavenr ) / (DOUBLE)length;
      end_v = BVENDVECTOR( bv_i );
      /*printf("tensor %g\n", tensor );*/
      BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv_i), end_v )                    /* over all points in the line */
      {
        VVALUE( v, tv_comp ) = tensor * sin( pos );
        pos += hkpi;
      }
      plane_pos += plane_hkpi;
      tensor = sin( plane_pos );
    }

  }
}


/****************************************************************************/
/*D
   TFFCalculateTheta - calculate the coupling matrix Theta

   SYNOPSIS:
   INT TFFCalculateTheta( const BLOCKVECTOR *bv_dest,
                                                  const BLOCKVECTOR *bv_source,
                                                  const BV_DESC *bvd_dest,
                                                  const BV_DESC *bvd_source,
                                                  const BV_DESC_FORMAT *bvdf,
                                                  INT Theta,
                                                  INT Tinv,
                                                  INT L,
                                                  INT tv_comp,
                                                  INT aux_comp,
                                                  INT auxsub_comp,
                                                  INT Lsub_comp );

   PARAMETERS:
   .  bv_dest - row-blockvector of Theta
   .  bv_source - column-blockvector of Theta (where the testvector lives)
   .  bvd_dest - description of the row-blockvector
   .  bvd_source - description of the column-blockvector
   .  bvdf - format to interpret the 'bvd's
   .  Theta - position of the Theta-component in the MATRIX-data
   .  Tinv - position of the component for the LU-decomposed diagonal blocks
   .  L - position of the off-diagonal component of the matrix to be decomposed
   .  tv_comp - position of the testvector-component in the VECTOR-data
   .  aux_comp - position of the auxiliary-component in the VECTOR-data
   .  auxsub_comp - position of the auxiliary-component for subproblems (only 3D)
   .  Lsub - position of the off-diagonal component of the subproblem-matrix to be decomposed (only 3D)

   DESCRIPTION:
   The diagonal entries of matrix Theta and Theta^T are calculated such that
   they fullfil the filtering condition (d=dest, s=source)
           Theta_(d,s) * tv_s == T_(d,d)^-1 * L_(d,s) * tv_s
   If a entry of the testvector is 0 the associated Theta-entry is set as the
   average of the neighbouring entries (for details see the code).

   To do certain experiments there can be activated variants of the usual
   algorithm by defining macro-names. See the code.

   'auxsub_comp' and 'Lsub_comp' are used only in the 3D case; in 2D you can
   use an arbitrary number.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if the testvector is constant 0.0
   D*/
/****************************************************************************/

INT TFFCalculateTheta( const BLOCKVECTOR *bv_dest, const BLOCKVECTOR *bv_source, const BV_DESC *bvd_dest, const BV_DESC *bvd_source, const BV_DESC_FORMAT *bvdf, INT Theta, INT Tinv, INT L, INT tv_comp, INT aux_comp, INT auxsub_comp, INT Lsub_comp )
/* calculates Theta_(i,i-1) and Theta_(i-1,i) */
/* only for regular grids */
/* auxsub and Lsub only in the 3D case */
{
  register VECTOR *v_dest, *v_source, *end_dest, *pred_dest, *pred_source, *succ_dest, *succ_source, *start_dest;
  register DOUBLE val, pred_val, succ_val;
  register MATRIX *m;
  INT missed, pred_found, succ_found, pred_in_block, succ_in_block;

  /* aux_source := L_(source,dest) * tv_dest */
  dsetBS( bv_source, aux_comp, 0.0 );
  dmatmulBS( bv_source, bvd_dest, bvdf, aux_comp, L, tv_comp );

  /* aux_source = (T_source)^-1 * aux_source */
#ifdef THETA_EXACT
  dcopyBS( bv_source, aux4_COMP, aux_comp );
  dsetBS( bv_source, aux_comp, 0.0 );
  gs_solveBS ( bv_source, bvd_source, bvdf, 1e-16, 100, Theta, aux_comp, aux4_COMP, aux5_COMP, TRUE );
#else
  TFFMultWithMInv( bv_source, bvd_source, bvdf, aux_comp, Lsub_comp, Tinv, aux_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );
#endif

  /* calculate Theta */
  /* Theta must fulfill the equation Theta * tv = aux = T^-1 * L * tv */

  missed = 0;
  end_dest = BVENDVECTOR( bv_dest );
  for ( v_dest = start_dest = BVFIRSTVECTOR( bv_dest ), v_source = BVFIRSTVECTOR( bv_source );
        v_dest != end_dest;
        v_dest = SUCCVC( v_dest ), v_source = SUCCVC( v_source ) )
  {
    val = VVALUE( v_dest, tv_comp );
    if ( fabs( val ) < TFFsmallTV )
    {
      SETVCUSED( v_dest, TRUE );
      missed++;
    }
    else
    {
      SETVCUSED( v_dest, FALSE );
      m = GetMatrix( v_source, v_dest );
      assert( m != NULL );
      MVALUE( m, Theta ) = MVALUE( MADJ(m), Theta ) = VVALUE( v_source, aux_comp ) / val;
    }
  }

  /* treat the missed members */

  v_dest = start_dest;
  v_source = BVFIRSTVECTOR( bv_source );
  start_dest = PREDVC( start_dest );

  /*printf( "%4d. block row, theta = %12g lambda = %12g\n", BVNUMBER( bv_dest ), MVALUE(GetMatrix(v_source,v_dest),Theta),1/MVALUE(GetMatrix(v_source,v_dest),Theta));*/

  while( missed > 0 )
  {
    if ( VCUSED( v_dest ) )
    {
#ifndef NDEBUG
      if ( mute_level >= 50 )
        UserWrite( "Missed vector in TFFCalculateTheta.\n" );
#endif
      /* the vector was missed */
      /* search a not missed neighbor to calculate an approximate value */

      missed--;
      pred_dest = succ_dest = v_dest;
      pred_source = succ_source = v_source;
      pred_found = succ_found = FALSE;
      pred_in_block = succ_in_block = TRUE;

      while ( !pred_found && !succ_found && (pred_in_block || succ_in_block) )
      {
        if ( !VCUSED( pred_dest ) && pred_in_block )
        {
          pred_val = MVALUE( GetMatrix( pred_source, pred_dest ), Theta );
          pred_found = TRUE;
        }

        if ( !VCUSED( succ_dest ) && succ_in_block )
        {
          succ_val = MVALUE( GetMatrix( succ_source, succ_dest ), Theta );
          succ_found = TRUE;
        }

        if ( pred_in_block )
        {
          pred_dest = PREDVC( pred_dest );
          pred_source = PREDVC( pred_source );
          pred_in_block = (pred_dest != start_dest);
        }

        if ( succ_in_block )
        {
          succ_dest = SUCCVC( succ_dest );
          succ_source = SUCCVC( succ_source );
          succ_in_block = (succ_dest != end_dest);
        }
      }                   /* while !pred_found */

      if ( pred_found )
        if ( succ_found )
        {
          if ( fabs(pred_val) > TFFmuchBigger * fabs(succ_val ) )
            val = succ_val;
          else if ( fabs(succ_val) > TFFmuchBigger * fabs(pred_val ) )
            val = pred_val;
          else
            val = (pred_val + succ_val) * 0.5;
        }
        else
          val = pred_val;
      else
      if ( succ_found )
        val = succ_val;
      else
      {
        UserWrite( "Testvector was zero in TFFCalculateTheta.\n" );
        m = GetMatrix( v_source, v_dest );
        assert( m != NULL );
        MVALUE( m, Theta ) = MVALUE( MADJ(m), Theta ) = 1e11;
        return NUM_ERROR;
      }

      m = GetMatrix( v_source, v_dest );
      assert( m != NULL );
      MVALUE( m, Theta ) = MVALUE( MADJ(m), Theta ) = val;

    }             /* if !VCUSED */

    v_dest = SUCCVC( v_dest );
    v_source = SUCCVC( v_source );

  }       /* while missed */

  return NUM_OK;
}


/****************************************************************************/
/*D
   TFFUpdateDiagBlock - calculate the filtered diagonal block

   SYNOPSIS:
   INT TFFUpdateDiagBlock( const BLOCKVECTOR *bv_dest,
                                                   const BV_DESC *bvd_dest,
                                                   const BV_DESC *bvd_source,
                                                   const BV_DESC_FORMAT *bvdf,
                                                   INT T,
                                                   INT DL,
                                                   INT Theta,
                                                   GRID *grid )

   PARAMETERS:
   .  bv_dest - new diagonal block and row-blockvector of Theta
   .  bvd_dest - description of the row-blockvector
   .  bvd_source - description of the last diagonal block and column-blockvector of Theta
   .  bvdf - format to interpret the 'bvd's
   .  T - position of the component of the filtered matrix in the MATRIX-data
   .  DL - position of the component of the given (stiffness) matrix
   .  Theta - position of the Theta-component in the MATRIX-data
   .  grid - grid where the matrixes lives on (for allocating extra connections)

   DESCRIPTION:
   Calculates for given Theta and T (d=dest, s=source)
                T_dd += Theta_ds*T_ss*Theta_sd - Theta_ds*L_sd - L_ds*Theta_sd
   If neccessary new connections are allocated as 'extra connections'.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if the testvector is constant 0.0
   D*/
/****************************************************************************/

INT TFFUpdateDiagBlock( const BLOCKVECTOR *bv_dest, const BV_DESC *bvd_dest, const BV_DESC *bvd_source, const BV_DESC_FORMAT *bvdf, INT T, INT DL, INT Theta, GRID *grid )
{
  /* T_dest += Theta_(dest,source)*T_source*Theta_(source,dest) */
  d3matmulBS( bv_dest, bvd_source, bvd_source, bvd_dest, bvdf, T, Theta, T, Theta, grid );

  /* T_dest -= Theta_(dest,source)*L_(source,dest) */
  d2matmul_minusBS( bv_dest, bvd_source, bvd_dest, bvdf, T, Theta, DL, grid );

  /* T_dest -= L_(dest,source)*Theta_(source,dest) */
  d2matmul_minusBS( bv_dest, bvd_source, bvd_dest, bvdf, T, DL, Theta, grid );

  return NUM_OK;
}


/****************************************************************************/
/*D
   TFFMultWithM - for a frequency filtered matrix M calculate y := M * x

   SYNOPSIS:
   INT TFFMultWithM( const BLOCKVECTOR *bv,
                                         const BV_DESC *bvd,
                                         const BV_DESC_FORMAT *bvdf,
                                         INT y_comp,
                                         INT T_comp,
                                         INT L_comp,
                                         INT Tinv_comp,
                                         INT x_comp,
                                         INT aux_comp,
                                         INT auxsub_comp,
                                         INT Lsub_comp );

   PARAMETERS:
   .  bv - blockvector covering M
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  y_comp - position of the result in the VECTOR-data
   .  T_comp - position of the frequency filtered diagonal blocks of the stiffness matrix
   .  L_comp - position of the off-diagonal blocks of the stiffness matrix
   .  Tinv_comp - position of the LU-decomposed diagonal blocks
   .  x_comp - position of the vector to be multiplied in the VECTOR-data
   .  aux_comp - position of the auxiliary-component in the VECTOR-data
   .  auxsub_comp - position of the auxiliary-component for subproblems (only 3D)
   .  Lsub_comp - position of the off-diagonal blocks of the subproblem-matrix (only 3D)

   DESCRIPTION:
   'M' is an frequency filtered decomposed matrix as produced by 'TFFDecomp',
   i.e. M = ( L + T ) * T^-1 * ( T + U ).
   This function calculates the product 'y := M * x'. Do not mix up this
   function with 'TFFMultWithMInv'!

   Both the lower ('L') and upper ('U') off-diagonal part is stored in the
   same component 'L_comp' in the MATRIX-data.

   'auxsub_comp', 'Tsub_comp' and 'Lsub_comp' are used only in the 3D case;
   in 2D you can use an arbitrary number.

   WARNING:
   'x_comp', 'y_comp' and 'aux_comp' must be pairwise different!

   REMARK:
   This function is needed only for experimenting, not for the pure solving
   of a PDE.

   SEE ALSO:
   TFFDecomp, TFFMultWithMInv

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   D*/
/****************************************************************************/

INT TFFMultWithM( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT y_comp, INT T_comp, INT L_comp, INT Tinv_comp, INT x_comp, INT aux_comp, INT auxsub_comp, INT Lsub_comp )
{
  register BLOCKVECTOR *bv_i, *bv_ip1, *bv_stop;
  register BV_DESC *bvd_i, *bvd_ip1, *bvd_temp;
  BV_DESC bvd1, bvd2;

  /* To minimize the incrementation of BVDs there are used two (one for
     index i and the other for index i+1) which are swapped in the loop;
     thus only one incrementation by 2 is necessary */
  bvd1 = bvd2 = *bvd;
  bvd_i = &bvd1;
  bvd_ip1 = &bvd2;
  BVD_PUSH_ENTRY( bvd_i, 0, bvdf );
  BVD_PUSH_ENTRY( bvd_ip1, 1, bvdf );

  /* calculate aux := T^-1*(U + T) * x = (T^-1 * U) * x + x */
  bv_stop = BVDOWNBVLAST( bv );
  for ( bv_i = BVDOWNBV(bv), bv_ip1 = BVSUCC( bv_i );
        bv_i != bv_stop;                        /* except the last blockvector */
        bv_i = bv_ip1, bv_ip1 = BVSUCC( bv_ip1 ) )
  {
    /* aux_i := 0 */
    dsetBS( bv_i, aux_comp, 0.0 );

    /* aux_i += U_(i,i+1) * x_i+1 */
    dmatmulBS( bv_i, bvd_ip1, bvdf, aux_comp, L_comp, x_comp );

    /* aux_i := (T_i)^-1 * aux_i */
    TFFMultWithMInv( bv_i, bvd_i, bvdf, aux_comp, Lsub_comp, Tinv_comp, aux_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );

    /* aux_i += x_i */
    daddBS( bv_i, aux_comp, x_comp );

    /* prepare BVDs for next loop */
    SWAP( bvd_i, bvd_ip1, bvd_temp );
    BVD_INC_LAST_ENTRY( bvd_ip1, 2, bvdf );
  }
  /* aux_last :=  T^-1*(U + T) * x_last = x_last */
  dcopyBS( bv_i, aux_comp, x_comp );

  bv_ip1 = bv_i;
  bv_i = BVPRED(bv_ip1);
  SWAP( bvd_i, bvd_ip1, bvd_temp );
  BVD_DEC_LAST_ENTRY( bvd_i, 2, bvdf );
  /* now i+1 points to the last blockvector, i to its predecessor */

  /* calculate y := (L + T) * aux */
  /* note the reverse direction of the calculation; this prevents y_i-1
     from beeing overriden in the i.th step, since it is needed in the
     step i-1 again. */
  bv_stop = BVDOWNBV(bv);
  for ( ;
        bv_ip1 != bv_stop;              /* except blockvector 0 */
        bv_ip1 = bv_i, bv_i = BVPRED( bv_ip1 ) )
  {
    /* y_i+1 := 0 */
    dsetBS( bv_ip1, y_comp, 0.0 );

    /* y_i+1 += T_i+1 * aux_i+1 */
    dmatmulBS( bv_ip1, bvd_ip1, bvdf, y_comp, T_comp, aux_comp );

    /* y_i+1 += L_(i+1,i) * aux_i */
    dmatmulBS( bv_ip1, bvd_i, bvdf, y_comp, L_comp, aux_comp );

    /* prepare BVDs for next loop */
    SWAP( bvd_i, bvd_ip1, bvd_temp );
    BVD_DEC_LAST_ENTRY( bvd_i, 2, bvdf );
  }
  /* y_0 := 0 */
  dsetBS( bv_ip1, y_comp, 0.0 );

  /* y_0 += T_0 * aux_0 */
  dmatmulBS( bv_ip1, bvd_ip1, bvdf, y_comp, T_comp, aux_comp );

  return NUM_OK;
}


/****************************************************************************/
/*D
   TFFMultWithMInv - for a frequency filtered matrix M calculate v := M^-1 * b

   SYNOPSIS:
   INT TFFMultWithMInv( const BLOCKVECTOR *bv,
                                                const BV_DESC *bvd,
                                                const BV_DESC_FORMAT *bvdf,
                                                INT v_comp,
                                                INT L_comp,
                                                INT Tinv_comp,
                                                INT b_comp,
                                                INT aux_comp,
                                                INT auxsub_comp,
                                                INT Lsub_comp );

   PARAMETERS:
   .  bv - blockvector covering M
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  v_comp - position of the result in the VECTOR-data
   .  L_comp - position of the off-diagonal blocks of the stiffness matrix
   .  Tinv_comp - position of the LU-decomposed diagonal blocks
   .  b_comp - position of the right hand side vector in the VECTOR-data
   .  aux_comp - position of the auxiliary-component in the VECTOR-data
   .  auxsub_comp - position of the auxiliary-component for subproblems (only 3D)
   .  Lsub_comp - position of the off-diagonal blocks of the subproblem-matrix (only 3D)

   DESCRIPTION:
   'M' is an frequency filtered decomposed matrix as produced by 'TFFDecomp',
   i.e. M = ( L + T ) * T^-1 * ( T + U ).
   This function calculates 'v := M^-1 * b', i.e. solves 'M * v = b' for v.
   Do not mix up this function with 'TFFMultWithM'!

   Both the lower ('L') and upper ('U') off-diagonal part is stored in the
   same component 'L_comp' in the MATRIX-data.

   'v_comp' and 'b_comp' may be equal.

   Destroys 'b_comp'.

   'auxsub_comp' and 'Lsub_comp' are used only in the 3D case;
   in 2D you can use an arbitrary number.
   The distinction whether 'bv' represents a plane (2D) or a cube (3D) is made
   dynamically by the blockvector structure itself: if 'bv' has only one
   further level of blockvectors (representing grid lines!) it is
   considered as a plane, otherwise as a cube.

   To do certain experiments there can be activated variants of the usual
   algorithm by defining macro-names. See the code.

   SEE ALSO:
   TFFDecomp, TFFMultWithM

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    error code from 'solveLUMatBS'
   D*/
/****************************************************************************/

INT TFFMultWithMInv( const BLOCKVECTOR *bv,
                     const BV_DESC *bvd,
                     const BV_DESC_FORMAT *bvdf,
                     INT v_comp,
                     INT L_comp,
                     INT Tinv_comp,
                     INT b_comp,
                     INT aux_comp,
                     INT auxsub_comp,
                     INT Lsub_comp )
{
  register BLOCKVECTOR *bv_i, *bv_ip1, *bv_stop;
  register BV_DESC *bvd_i, *bvd_ip1, *bvd_temp;
  BV_DESC bvd1, bvd2;

  if ( BV_IS_LEAF_BV(bv) )
    return solveLUMatBS( bv, bvd, bvdf, v_comp, Tinv_comp, b_comp );

  ASSERT( v_comp != aux_comp );
  ASSERT( b_comp != aux_comp );
  ASSERT( L_comp != Tinv_comp );

  IFDEBUG(np,0)
  if ( auxsub_comp != DUMMY_COMP )
  {
    ASSERT( auxsub_comp != v_comp );
    ASSERT( auxsub_comp != b_comp );
    ASSERT( auxsub_comp != aux_comp );
    ASSERT( Lsub_comp != DUMMY_COMP );
    ASSERT( Lsub_comp != L_comp );
    ASSERT( Lsub_comp != Tinv_comp );
  }
  ENDDEBUG

  /* To minimize the incrementation of BVDs there are used two (one for
     index i and the other for index i+1) which are swapped in the loop;
     thus only one incrementation by 2 is necessary */
    bvd1 = bvd2 = *bvd;
  bvd_i = &bvd1;
  bvd_ip1 = &bvd2;
  BVD_PUSH_ENTRY( &bvd1, 0, bvdf );
  BVD_PUSH_ENTRY( &bvd2, 1, bvdf );

  /* solve lower triangular matrix; except the last block in this loop */
  /* aux := ( L + T )^-1 * b */
  bv_stop = BVDOWNBVLAST(bv);
  for ( bv_i = BVDOWNBV(bv), bv_ip1 = BVSUCC( bv_i );
        bv_i != bv_stop;
        bv_i = bv_ip1, bv_ip1 = BVSUCC( bv_ip1 ) )
  {
    /* aux_i := (T_i)^-1 * b_i */
#ifdef MINV_2D_EXACT
    if( L_comp == 0 )
      gs_solveBS ( bv_i, bvd_i, bvdf, 1e-16, 100, FF_COMP, aux_comp, b_comp, aux5_COMP, TRUE );
    else
      /*gs_solveBS ( bv_i, bvd_i, bvdf, 1e-16, 100, FF3D_COMP, aux_comp, b_comp, aux5_COMP, TRUE );*/
      TFFMultWithMInv( bv_i, bvd_i, bvdf, aux_comp, Lsub_comp, Tinv_comp, b_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );
#else
    TFFMultWithMInv( bv_i, bvd_i, bvdf, aux_comp, Lsub_comp, Tinv_comp, b_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );
#endif

    /* b_i+1 -= L_(i+1,i) * aux_i */
    dmatmul_minusBS( bv_ip1, bvd_i, bvdf, b_comp, L_comp, aux_comp );

    /* prepare BVDs for next loop */
    SWAP( bvd_i, bvd_ip1, bvd_temp );
    BVD_INC_LAST_ENTRY( bvd_ip1, 2, bvdf );
  }
  /* special treatment: v_last = (T_last)^-1 * b_last */
#ifdef MINV_2D_EXACT
  if( L_comp == 0 )
    gs_solveBS ( bv_i, bvd_i, bvdf, 1e-16, 100, FF_COMP, aux_comp, b_comp, aux5_COMP, TRUE );
  else
    /*gssolveBS ( bv_i, bvd_i, bvdf, 1e-16, 100, FF3D_COMP, aux_comp, b_comp, aux5_COMP, TRUE );*/
    TFFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, Lsub_comp, Tinv_comp, b_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );
#else
  TFFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, Lsub_comp, Tinv_comp, b_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );
#endif

  /* solve upper triangular matrix; the last block is already calculated */
  /* v := (T^-1*U + I )^-1 * aux */
  assert( bv_i == BVDOWNBVLAST(bv) );
  SWAP( bvd_i, bvd_ip1, bvd_temp );
  BVD_DEC_LAST_ENTRY( bvd_i, 2, bvdf );
  bv_stop = BVPRED( BVDOWNBV(bv) );
  for ( bv_i = BVPRED( bv_i ); bv_i != bv_stop; bv_i = BVPRED( bv_i ) )
  {
    /* v_i := L_(i,i+1) * v_i+1 */
    dsetBS( bv_i, v_comp, 0.0 );
    dmatmulBS( bv_i, bvd_ip1, bvdf, v_comp, L_comp, v_comp );

    /* v_i := (T_i)^-1 * v_i */
#ifdef MINV_2D_EXACT
    dcopyBS( bv_i, aux4_COMP, v_comp );
    if( L_comp == 0 )
      gs_solveBS ( bv_i, bvd_i, bvdf, 1e-16, 100, FF_COMP, v_comp, aux4_COMP, aux5_COMP, TRUE );
    else
      /*gs_solveBS ( bv_i, bvd_i, bvdf, 1e-16, 100, FF3D_COMP, v_comp, aux4_COMP, aux5_COMP, TRUE );*/
      TFFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, Lsub_comp, Tinv_comp, v_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );
#else
    TFFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, Lsub_comp, Tinv_comp, v_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP );
#endif

    /* v_i := aux_i - v_i */
    dminusaddBS( bv_i, v_comp, aux_comp );

    /* prepare BVDs for next loop */
    SWAP( bvd_i, bvd_ip1, bvd_temp );
    BVD_DEC_LAST_ENTRY( bvd_i, 2, bvdf );
  }

  return NUM_OK;
}

/****************************************************************************/
/*																			*/
/* global exported routines													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   TFFDecomp - calculate the tangential frequency filtering decomposition of a matrix

   SYNOPSIS:
   INT TFFDecomp( DOUBLE wavenr,
                                  DOUBLE wavenr3D,
                                  const BLOCKVECTOR *bv,
                                  const BV_DESC *bvd,
                                  const BV_DESC_FORMAT *bvdf,
                                  INT LU_comp,
                                  INT FF_comp,
                                  INT K_comp,
                                  INT tv_comp,
                                  INT aux_comp,
                                  INT auxsub_comp,
                                  INT FFsub_comp,
                                  GRID *grid );

   PARAMETERS:
   .  wavenr - number of sine half-oszillations along one gridline (2D (sub)problems)
   .  wavenr3D - number of sine half-oszillations along the planes (only 3D problems)
   .  bv - blockvector covering the matrix to be decomposed
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  LU_comp - position of the component for the LU-decomposed (diagonal) blocks
   .  FF_comp - position of the Theta and filtered diagonal blocks in the MATRIX-data
   .  K_comp - position of the stiffness matrix to be decomposed
   .  tv_comp - position of the testvector-component in the VECTOR-data
   .  aux_comp - position of the auxiliary-component in the VECTOR-data
   .  auxsub_comp - position of the auxiliary-component for subproblems (only 3D)
   .  FFsub_comp - position of the Theta and filtered diagonal blocks of the subproblem-matrix (only 3D)
   .  grid - grid where the matrixes lives on (for allocating extra connections)

   DESCRIPTION:
   This function calculates the tangential frequency filtered decomposition
   of the given matrix 'K'. It expects a blockvector structure according to a
   linewise (2D) resp. plane-/linewise (3D) decomposition of the domain as
   it is constructed by 'CreateBVStripe2D'/'CreateBVStripe3D'. The result is
   a matrix 'M' in the form 'M = ( L + T ) * T^-1 * ( T + U )' where the
   lower ('L') and upper ('U') off-diagonal blocks are the same as in
   the given matrix ('K_comp'), 'T' is the frequency filtered diagonal
   blockmatrix ('FF_comp'). The auxiliary matrix Theta is stored in the off-diagonal
   blocks of 'FF_comp'.

   In the 3D case the frequency filtered diagonal blocks are again decomposed
   by TFF. The matrix to be decomposed is in the 'FF_comp' and the frequency
   filtered diagonal blocks are stored in 'FFsub_comp' (their LU decomposition
   finally in 'LU_comp'). 'auxsub_comp' and 'FFsub_comp' are used only in
   the 3D case; in 2D you can use an arbitrary number.

   The distinction whether 'bv' represents a plane (2D) or a cube (3D) is made
   dynamically by the blockvector structure itself: if 'bv' has only one
   further level of blockvectors (representing grid lines!) it is
   considered as a plane, otherwise as a cube.

   To do certain experiments there can be activated variants of the usual
   algorithm by defining macro-names. See the code.

   SEE ALSO:
   TFFCalculateTheta, TFFConstructTestvector_loc, TFFConstructTestvector, TFFMultWithM, TFFMultWithMInv

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    error code from 'LUDecomposeDiagBS'
   D*/
/****************************************************************************/

INT TFFDecomp( DOUBLE wavenr,
               DOUBLE wavenr3D,
               const BLOCKVECTOR *bv,
               const BV_DESC *bvd,
               const BV_DESC_FORMAT *bvdf,
               INT LU_comp,
               INT FF_comp,
               INT K_comp,
               INT tv_comp,
               INT aux_comp,
               INT auxsub_comp,
               INT FFsub_comp,
               GRID *grid )
{
  register BLOCKVECTOR *bv_i, *bv_im1, *bv_end;
  register BV_DESC *bvd_i, *bvd_im1, *bvd_temp;
  BV_DESC bvd1, bvd2;
#ifdef THETA_ANA
  DOUBLE lambda, li, a, bx, by, bz;
  MATRIX *m;
#endif

  if ( BV_IS_LEAF_BV(bv) )
  {
    dmatcopyBS( bv, bvd, bvdf, LU_comp, K_comp );
    return LUDecomposeDiagBS( bv, bvd, bvdf, LU_comp, grid );
  }

  /* initialize the BVDs */
  bvd1 = bvd2 = *bvd;                   /* copy of the BVD */
  bvd_i = &bvd1;
  bvd_im1 = &bvd2;
  BVD_PUSH_ENTRY( bvd_im1, 0, bvdf );
  BVD_PUSH_ENTRY( bvd_i, 1, bvdf );

#ifdef THETA_ANA
  bv_im1 = BVDOWNBV(bv);
  bv_i = BVSUCC( bv_im1 );
  a = MVALUE( VSTART(BVFIRSTVECTOR( bv_im1 )), K_comp );
  m = GetMatrix( BVFIRSTVECTOR(bv_im1), SUCCVC(BVFIRSTVECTOR(bv_im1)));
  assert( m != NULL );
  bx = MVALUE( m, K_comp );
  if ( BV_IS_LEAF_BV( bv_im1 ) )
  {             /* 2D block */
    m = GetMatrix( BVFIRSTVECTOR(bv_im1), BVFIRSTVECTOR(bv_i));
    assert( m != NULL );
    by = MVALUE( m, K_comp );
    li = lambda = (a + 2.0 * bx * cos(wavenr*PI*glob_h)) / by;
    /*printf("a = %g bx = %g by = %g  lambda = %g  h= %g i = %g\n", a,bx, by,lambda, glob_h, wavenr );*/
  }
  else
  {             /* 3D block */
    m = GetMatrix( BVFIRSTVECTOR(BVDOWNBV(bv_im1)), BVFIRSTVECTOR(BVSUCC(BVDOWNBV(bv_im1))));
    assert( m != NULL );
    by = MVALUE( m, K_comp );
    m = GetMatrix( BVFIRSTVECTOR(bv_im1), BVFIRSTVECTOR(bv_i));
    assert( m != NULL );
    bz = MVALUE( m, K_comp );
    li = lambda = (a + 2.0 * bx * cos(wavenr*PI*glob_h) + 2.0 * by * cos(wavenr3D*PI*glob_h)) / bz;
    /*printf("a = %g  bx = %g  by = %g  bz = %g  lambda = %g  h= %g i = %g  j = %g\n", a,bx,by,bz,lambda, glob_h, wavenr, wavenr3D );*/
  }
#endif

  /* the loop calculates (T_i-1)^-1 and Theta(i,j) for i=1..n-1 */
  /* T_0 := D_0 */
  bv_im1 = BVDOWNBV(bv);
  dmatcopyBS( bv_im1, bvd_im1, bvdf, FF_comp, K_comp );
  bv_end = BVDOWNBVEND(bv);
  for ( bv_i = BVSUCC( bv_im1 );
        bv_i != bv_end;
        bv_im1 = bv_i, bv_i = BVSUCC( bv_i ) )
  {
    /* T_i-1 decompose */
    TFFDecomp( wavenr, wavenr3D, bv_im1, bvd_im1, bvdf, LU_comp, FFsub_comp, FF_comp, tv_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP, grid );

    /* calculate Theta_(i-1,i) and Theta_(i,i-1);
       result on off-digonal blocks of FF_comp */
    TFFConstructTestvector_loc( bv_i, tv_comp, wavenr, wavenr3D );

#ifdef THETA_ANA
    dmatsetBS( bv_i, bvd_im1, bvdf, FF_comp, 1.0/li );
    dmatsetBS( bv_im1, bvd_i, bvdf, FF_comp, 1.0/li );
    /*printf("theta = %g\n", 1.0 / li );*/
    li = lambda - 1.0 / li;
#else
    TFFCalculateTheta( bv_i, bv_im1, bvd_i, bvd_im1, bvdf, FF_comp, LU_comp, K_comp, tv_comp, aux_comp, auxsub_comp, FF_comp );
#endif

    /* T_i := D_i + Theta_(i,i-1)*T_i-1*Theta_(i-1,i) -
                      Theta_(i,i-1)*L_(i-1,i) - L_(i,i-1)*Theta_(i-1,i) */
    dmatcopyBS( bv_i, bvd_i, bvdf, FF_comp, K_comp );
    TFFUpdateDiagBlock( bv_i, bvd_i, bvd_im1, bvdf, FF_comp, K_comp, FF_comp, grid );

    /* update BVDs for the next loop */
    SWAP( bvd_i, bvd_im1, bvd_temp );
    BVD_INC_LAST_ENTRY( bvd_i, 2, bvdf );
  }
  /* now bv_im1 and bvd_im1 points to the last block */

  /* calculate the last (T_n)^-1 */
  TFFDecomp( wavenr, wavenr3D, bv_im1, bvd_im1, bvdf, LU_comp, FFsub_comp, FF_comp, tv_comp, auxsub_comp, DUMMY_COMP, DUMMY_COMP, grid );


#ifdef CHECK_CALCULATION
  {
    /* notice: for the first block exists no testvector. But it is still 0
                   thus it do not disturb the check */
    DOUBLE norm;

    /* to get the right global testvector for the tests */
    if( K_comp == 0 )
      TFFConstructTestvector( bv, tv_comp, wavenr, wavenr3D );

    /* aux2 := M * tv */
    TFFMultWithM( bv, bvd, bvdf, aux2_COMP, FF_comp, K_comp, LU_comp, tv_comp, aux_comp, auxsub_comp, FF_comp );

    /* aux := K * tv */
    dsetBS( bv, aux_comp, 0.0 );
    dmatmulBS( bv, bvd, bvdf, aux_comp, K_comp, tv_comp );

    /* aux2 -= aux */
    dsubBS( bv, aux2_COMP, aux_comp );

    /* norm of aux2 */
    eunormBS( bv, aux2_COMP, &norm );

    if ( norm > TFFaccuracy )
    {
      printf( "approximation M *tv == K * tv (%g) is not good enough!!!!!!!!!!!!!!!!!\n", norm );
      UserWriteF( "approximation M *tv == K * tv (%g) is not good enough!!!!!!!!!!!!!!!!!\n", norm );
    }
    else
    if ( mute_level >= 50 )

    {
      printf( "approximation M *tv == K * tv (%g) is good\n", norm );
      UserWriteF( "approximation M *tv == K * tv (%g) is good\n", norm );
    }
  }
#endif

  /*printvBS(bv, tv_comp); printmBS(bv,bv,K_comp); printmBS(bv,bv,FF_comp); printmBS(bv,bv,LU_comp); printf("\n");*/
  return NUM_OK;
}


/****************************************************************************/
/*D
   TFFSolve - solves an linear equation with the tangential frequency filtering method

   SYNOPSIS:
   INT  TFFSolve( const BLOCKVECTOR *bv,
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
                                  INT auxsub_comp,
                                  INT FFsub_comp,
                                  DOUBLE meshwidth,
                                  DOUBLE eps,
                                  GRID *grid );

   PARAMETERS:
   .  bv - blockvector covering the matrix to be decomposed
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  K_comp - position of the stiffness matrix in the MATRIX-data to be decomposed
   .  u_comp - position of the solution in the VECTOR-data
   .  f_comp - position of the right hand side in the VECTOR-data
   .  cor_comp - position of the correction of the solution in the VECTOR-data
   .  FF_comp - position of the Theta and filtered diagonal blocks in the MATRIX-data
   .  LU_comp - position of the component for the LU-decomposed (diagonal) blocks
   .  tv_comp - position of the testvector-component in the VECTOR-data
   .  aux_comp - position of the auxiliary-component in the VECTOR-data
   .  auxsub_comp - position of the auxiliary-component for subproblems (only 3D)
   .  FFsub_comp - position of the Theta and filtered diagonal blocks of the subproblem-matrix (only 3D)
   .  meshwidth - meshwidth of the regular grid
   .  eps - accuracy of the final solution (stopping criterion)
   .  grid - grid where the matrixes lives on (for allocating extra connections)

   DESCRIPTION:
   This function solves the linear equation 'K*u = f' with the tangential
   frequency filtering method (see the literature). 'u_comp' contains the
   first approximation of the solution, the starting point of the iteration.
   'f_comp' imports the right hand side vector and will be destroyed during
   the iteration because it is used as the defect-vector.

   'auxsub_comp' and 'FFsub_comp' are used only in the 3D case; in 2D you
   can use an arbitrary number.

   To do certain experiments there can be activated variants of the usual
   algorithm by defining macro-names. See the code. If CHECK_CALC is active
   the constants aux2_COMP and aux3_COMP are used as further postions in the
   VECTOR-data for checking calculations. This is dirty!

   WARNING:
   Destroys the right hand side 'f_comp'.

   Works only for regular meshes. Blockvector structure according to
   'CreateBVStripe2D'/'CreateBVStripe3D' must be prepared.

   SEE ALSO:
   TFFCalculateTheta, TFFConstructTestvector_loc, TFFConstructTestvector, TFFMultWithM, TFFMultWithMInv, TFFDecomp, CreateBVStripe2D, CreateBVStripe3D

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   D*/
/****************************************************************************/

INT  TFFSolve( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT K_comp, INT u_comp, INT f_comp, INT cor_comp, INT FF_comp, INT LU_comp, INT tv_comp, INT aux_comp, INT auxsub_comp, INT FFsub_comp, DOUBLE meshwidth, DOUBLE eps, GRID *grid )
{
  DOUBLE old_norm, new_norm, start_norm, step_norm, final_acc;
  INT i, j, it, nr_TFFs;

  /* calculate i with 1/pow(2,i) == meshwidth */
  nr_TFFs = (INT)( log(1.0/meshwidth)/M_LN2 + 0.5 );
  UserWriteF( "meshwidth %g = 1/%g  nr_TFFs %d\n", meshwidth, 1.0/meshwidth, nr_TFFs );

  /* f := f - K * u */
  start_norm = new_norm = CalculateDefectAndNormBS( bv, bvd, bvdf, f_comp, f_comp, K_comp, u_comp );
  UserWriteF( "start defect %g\n", start_norm );
  /* now the rhs is destroyed and is used as the defect */

  /* different possibilities for accuracy */
  /*final_acc = eps * meshwidth * meshwidth;*/
  /*final_acc = eps * meshwidth * meshwidth * start_norm;*/
  final_acc = eps;
  it = 0;
  while ( new_norm > final_acc )
  {
    it++;

    old_norm = new_norm;

#ifdef __THREEDIM__
    for ( j = 0; j < nr_TFFs; j++ )
    {
      i = j;
#endif

#ifdef __TWODIM__
    for ( i = 0; i < nr_TFFs; i++ )
    {
#endif

#ifdef THETA_ANA
    glob_h = meshwidth;
#endif
    step_norm = new_norm;

    /* construct testvector */
    /*TFFConstructTestvector( bv, tv_comp, (double)(1<<i), (double)(1<<j) );*/

    /* Calculates the TFF decomposition of M = (L + T) * T^-1 * (L^T + T) */
    TFFDecomp( (double)(1<<i), (double)(1<<j), bv, bvd, bvdf, LU_comp, FF_comp, K_comp, tv_comp, aux_comp, auxsub_comp, FFsub_comp, grid );

    /* cor := M^-1 * f */
    dcopyBS( bv, cor_comp, f_comp );
    TFFMultWithMInv( bv, bvd, bvdf, cor_comp, K_comp, LU_comp, cor_comp, aux_comp, auxsub_comp, FF_comp );

#ifdef CHECK_CALCULATION
    {
      DOUBLE norm;

      /* check invertation of testvector tv == M^-1 * K * tv;
              must always be true since it is equivalent to
              S * tv == 0, i.e. solve exact on the tv-subspace*/

      /* aux2 := K * tv */
      dsetBS( bv, aux2_COMP, 0.0 );
      dmatmulBS( bv, bvd, bvdf, aux2_COMP, K_comp, tv_comp );

      /* aux2 := M^-1 * aux2 */
      TFFMultWithMInv( bv, bvd, bvdf, aux2_COMP, K_comp, LU_comp, aux2_COMP, aux3_COMP, auxsub_comp, FF_comp );

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      eunormBS( bv, aux2_COMP, &norm );

      if ( norm > TFFaccuracy )
      {
        printf( "M^-1*K*tv==tv is not good (%g)!!!!!!!!!!!!!!!!!\n", norm );
        UserWriteF( "M^-1*K*tv==tv is not good (%g)!!!!!!!!!!!!!!!!!\n", norm );
      }
      else
      if ( mute_level >= 50 )
      {
        printf( "M^-1*K*tv==tv is good  (%g)\n", norm );
        UserWriteF( "M^-1*K*tv==tv is good (%g)\n", norm );
      }

      /* tests how good tv == M^-1 * M * tv is fullfiled; must always */

      /* aux2 := M * tv */
      TFFMultWithM( bv, bvd, bvdf, aux2_COMP, FF_comp, K_comp, LU_comp, tv_comp, aux3_COMP, auxsub_comp, FF_comp );

      /* aux2 := M^-1 * aux2 */
      dcopyBS( bv, aux_comp, tv_comp );
      TFFMultWithMInv( bv, bvd, bvdf, aux2_COMP, K_comp, LU_comp, aux2_COMP, aux3_COMP, auxsub_comp, FF_comp );

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      eunormBS( bv, aux2_COMP, &norm );

      if ( norm > TFFaccuracy )
      {
        printf( "M^-1*M*tv==tv is not good enough (%g)!!!!!!!!!!!!!!!!!\n", norm );
        UserWriteF( "M^-1*M*tv==tv is not good enough (%g)!!!!!!!!!!!!!!!!!\n", norm );
      }
      else
      if ( mute_level >= 50 )
      {
        printf( "M^-1*M*tv==tv is good  (%g)\n", norm );
        UserWriteF( "M^-1*M*tv==tv is good (%g)\n", norm );
      }

      /* tests how good tv == M * M^-1 * tv is fullfiled */

      /* aux := M^-1 * tv */
      dcopyBS( bv, aux_comp, tv_comp );
      TFFMultWithMInv( bv, bvd, bvdf, aux_comp, K_comp, LU_comp, aux_comp, aux3_COMP, auxsub_comp, FF_comp );

      /* aux2 := M * aux */
      TFFMultWithM( bv, bvd, bvdf, aux2_COMP, FF_comp, K_comp, LU_comp, aux_comp, aux3_COMP, auxsub_comp, FF_comp );

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      eunormBS( bv, aux2_COMP, &norm );

      if ( norm > TFFaccuracy )
      {
        printf( "M*M^-1*tv==tv is not good enough (%g)..................\n", norm );
        UserWriteF( "M*M^-1*tv==tv is not good enough (%g)..................\n", norm );
      }
      else
      if ( mute_level >= 50 )
      {
        printf( "M*M^-1*tv==tv is good  (%g)\n", norm );
        UserWriteF( "M*M^-1*tv==tv is good (%g)\n", norm );
      }

      /* check invertation of testvector tv == K * M^-1 * tv */

      /* aux3 := M^-1 * tv */
      dcopyBS( bv, aux3_COMP, tv_comp );
      TFFMultWithMInv( bv, bvd, bvdf, aux3_COMP, K_comp, LU_comp, aux3_COMP, aux2_COMP, auxsub_comp, FF_comp );

      /* aux2 := K * aux3 */
      dsetBS( bv, aux2_COMP, 0.0 );
      dmatmulBS( bv, bvd, bvdf, aux2_COMP, K_comp, aux3_COMP );

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      eunormBS( bv, aux2_COMP, &norm );

      if ( norm > TFFaccuracy )
      {
        printf( "K*M^-1*tv==tv is not good (%g)..................\n", norm );
        UserWriteF( "K*M^-1*tv==tv is not good (%g)..................\n", norm );
      }
      else
      if ( mute_level >= 50 )
      {
        printf( "K*M^-1*tv==tv is good  (%g)\n", norm );
        UserWriteF( "K*M^-1*tv==tv is good (%g)\n", norm );
      }


      if ( mute_level >= 50 )
      {
        /* another test: compare with given stiffness matrix
           K * cor - M * cor */

        /* aux2 := K * cor */
        dsetBS( bv, aux2_COMP, 0.0 );
        dmatmulBS( bv, bvd, bvdf, aux2_COMP, K_comp, cor_comp );

        /* aux2 -= aux */
        dsubBS( bv, aux2_COMP, aux_comp );

        /* norm of aux2 */
        eunormBS( bv, aux2_COMP, &norm );

        if ( norm > TFFaccuracy )
        {
          printf( "Invertation is not good enough (%g) resp. K\n", norm );
          UserWriteF( "Invertation is not good enough (%g) resp. K\n", norm );
        }
        else
        {
          printf( "Invertation is good  (%g) resp. K\n", norm );
          UserWriteF( "Invertation is good (%g) resp. K\n", norm );
        }
      }
    }
#endif

    /* u := u + cor */
    daddBS( bv, u_comp, cor_comp );

    /* f := f - K * cor */
    new_norm = CalculateDefectAndNormBS( bv, bvd, bvdf, f_comp, f_comp, K_comp, cor_comp );

#ifdef __TWODIM__
    UserWriteF( "Wavenumber = %2d new defect = %12lg conv. rate = %12lg\n", 1<<i, new_norm, new_norm/step_norm );
#endif

#ifdef __THREEDIM__
    UserWriteF( "Wnr plane = %2d Wnr line = %2d new defect = %12lg conv. rate = %12lg\n", 1<<j, 1<<i, new_norm, new_norm/step_norm );
#endif
  }               /* for */

  UserWriteF( "new defect = %4lg conv. rate = %12lg\n", new_norm, new_norm/old_norm );
}         /* while */

UserWriteF( "avarage of convergency rate ( %d iterations) = %12lg\n", it, pow( new_norm / start_norm, 1.0 / (DOUBLE)it ) );

return NUM_OK;
}

#endif /* __BLOCK_VECTOR_DESC__ */
