// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ff.c                                                                                                                  */
/*																			*/
/* Purpose:   (tangential) frequency filtering decompostion routine         */
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
/* Remarks:   FF is used as the abbreviation for "frequency filtering" in	*/
/*			      genral and for Wittums method in particular special		*/
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
#include "ff.h"

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

REP_ERR_FILE;

static DOUBLE glob_h;   /* used if THETA_ANA is defined */

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

#ifdef __BLOCK_VECTOR_DESC__

INT TFFCalculateTheta( const BLOCKVECTOR *bv_dest, const BLOCKVECTOR *bv_source, const BV_DESC *bvd_dest, const BV_DESC *bvd_source, const BV_DESC_FORMAT *bvdf, INT tv_comp );

INT TFFUpdateDiagBlock( const BLOCKVECTOR *bv_dest, const BV_DESC *bvd_dest, const BV_DESC *bvd_source, const BV_DESC_FORMAT *bvdf, INT T, INT DL, INT Theta, GRID *grid );

INT InitFF( void );

/****************************************************************************/
/*D
   TFFCalculateTheta - calculate the coupling matrix Theta for tangential frequency filtering

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
   .  bv_source - column-blockvector of Theta
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
           Theta_(s,d) * tv_d == T_(s,s)^-1 * L_(s,d) * tv_d
   If a entry of the testvector is 0 the associated Theta-entry is set as the
   average of the neighbouring entries (for details see the code).

   To do certain experiments there can be activated variants of the usual
   algorithm by defining macro-names. See the code.

   'auxsub_comp' and 'Lsub_comp' are used only in the 3D case; in 2D you can
   use an arbitrary number.

   SEE ALSO:
   FFCalculateTheta

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if the testvector is constant 0.0
   D*/
/****************************************************************************/

INT TFFCalculateTheta( const BLOCKVECTOR *bv_dest,
                       const BLOCKVECTOR *bv_source,
                       const BV_DESC *bvd_dest,
                       const BV_DESC *bvd_source,
                       const BV_DESC_FORMAT *bvdf,
                       INT tv_comp )
{
  register VECTOR *v_dest, *v_source, *end_dest, *pred_dest, *pred_source, *succ_dest, *succ_source, *start_dest;
  register DOUBLE val, pred_val, succ_val;
  register MATRIX *m;
  INT missed, pred_found, succ_found, pred_in_block, succ_in_block;
  INT aux_comp, L_comp, Theta_comp;

  aux_comp = GET_AUX_VEC;
  L_comp = STIFFMAT_ON_LEVEL_BLOCKWISE(bv_dest);
  Theta_comp = DECOMPMAT_ON_LEVEL_BLOCKWISE(bv_dest);

  ASSERT( aux_comp != DUMMY_COMP );
  ASSERT( tv_comp != DUMMY_COMP );
  ASSERT( L_comp != DUMMY_COMP );
  ASSERT( Theta_comp != DUMMY_COMP );
  /*ASSERT( L_comp == STIFFMAT_ON_LEVEL_BLOCKWISE(bv_source) );*/
  /*ASSERT( Theta_comp == DECOMPMAT_ON_LEVEL_BLOCKWISE(bv_source) );*/

  /* aux_source := L_(source,dest) * tv_dest */
  dsetBS( bv_source, aux_comp, 0.0 );
  dmatmul_addBS( bv_source, bvd_dest, bvdf, aux_comp, L_comp, tv_comp );

  /* aux_source = (T_source)^-1 * aux_source */
#ifdef THETA_EXACT
  dcopyBS( bv_source, aux4_COMP, aux_comp );
  dsetBS( bv_source, aux_comp, 0.0 );
  gs_solveBS ( bv_source, bvd_source, bvdf, 1e-16, 100, Theta_comp, aux_comp, aux4_COMP, aux5_COMP, TRUE );
#else
        #ifdef ModelP
  FFMultWithMInv( bv_source, bvd_source, bvdf, aux_comp, aux_comp, NULL, NULL );
        #else
  FFMultWithMInv( bv_source, bvd_source, bvdf, aux_comp, aux_comp );
        #endif
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
    if ( fabs( val ) < FFsmallTV )
    {
      SETVCUSED( v_dest, TRUE );
      missed++;
    }
    else
    {
      SETVCUSED( v_dest, FALSE );
      m = GetMatrix( v_source, v_dest );
      ASSERT( m != NULL );
      MVALUE( m, Theta_comp ) = MVALUE( MADJ(m), Theta_comp ) = VVALUE( v_source, aux_comp ) / val;
    }
  }

  /* treat the missed members */

  v_dest = start_dest;
  v_source = BVFIRSTVECTOR( bv_source );
  start_dest = PREDVC( start_dest );

  /*printf( "%4d. block row, theta = %12g lambda = %12g\n", BVNUMBER( bv_dest ), MVALUE(GetMatrix(v_source,v_dest),Theta_comp),1/MVALUE(GetMatrix(v_source,v_dest),Theta_comp));*/

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
          pred_val = MVALUE( GetMatrix( pred_source, pred_dest ), Theta_comp );
          pred_found = TRUE;
        }

        if ( !VCUSED( succ_dest ) && succ_in_block )
        {
          succ_val = MVALUE( GetMatrix( succ_source, succ_dest ), Theta_comp );
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
          if ( fabs(pred_val) > FFmuchBigger * fabs(succ_val ) )
            val = succ_val;
          else if ( fabs(succ_val) > FFmuchBigger * fabs(pred_val ) )
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
        MVALUE( m, Theta_comp ) = MVALUE( MADJ(m), Theta_comp ) = 1e11;
        FREE_AUX_VEC(aux_comp);
        REP_ERR_RETURN(NUM_ERROR);
      }

      m = GetMatrix( v_source, v_dest );
      assert( m != NULL );
      MVALUE( m, Theta_comp ) = MVALUE( MADJ(m), Theta_comp ) = val;

    }             /* if !VCUSED */

    v_dest = SUCCVC( v_dest );
    v_source = SUCCVC( v_source );

  }       /* while missed */

  FREE_AUX_VEC(aux_comp);

  return(NUM_OK);
}


/****************************************************************************/
/*D
   FFCalculateThetaAndUpdate - calculate the FF approximation of a block

   SYNOPSIS:
   INT FFCalculateThetaAndUpdate( const BLOCKVECTOR *bv_dest,
                         const BLOCKVECTOR *bv_source,
                                                 const BV_DESC *bvd_dest,
                                                 const BV_DESC *bvd_source,
                                                 const BV_DESC_FORMAT *bvdf,
                                                 INT Theta_comp,
                                                 INT Tinv_comp,
                                                 INT K_comp,
                                                 INT tv1_comp,
                                                 INT tv2_comp,
                                                 INT aux1_comp,
                                                 INT aux2_comp,
                                                 INT auxsub_comp,
                                                 INT Ksub_comp,
                                                 GRID *grid);

   PARAMETERS:
   .  bv_dest - blockvector for testvector
   .  bv_source - previous blockvector
   .  bvd_dest - description of the destination-blockvector
   .  bvd_source - description of the source-blockvector
   .  bvdf - format to interpret the 'bvd's
   .  Theta_comp - position of the Theta-component in the MATRIX-data
   .  Tinv_comp - position of the component for the LU-decomposed diagonal blocks
   .  K_comp - position of the diagonal and off-diagonal component of the matrix to be decomposed
   .  tv1_comp - position of the 1. testvector-component in the VECTOR-data
   .  tv2_comp - position of the 2. testvector-component in the VECTOR-data
   .  aux1_comp - position of the 1. auxiliary-component in the VECTOR-data
   .  aux2_comp - position of the 2. auxiliary-component in the VECTOR-data
   .  auxsub_comp - position of the auxiliary-component for subproblems (only 3D)
   .  Ksub_comp - position of the diagonal and off-diagonal component of the subproblem-matrix to be decomposed (only 3D)
   .  grid - to allocate new entries if necessary

   DESCRIPTION:
   The tridiagonal matrix Theta is calculated such that
   they fullfil the filtering condition (d=dest, s=source)
           Theta * tv_d == L(d,s) * T_(s,s)^-1 * L_(s,d) * tv_d

   And then set
       T_(d,d) := D_(d,d) - Theta

   'auxsub_comp' and 'Lsub_comp' are used only in the 3D case; in 2D you can
   use an arbitrary number.

   SEE ALSO:
   TFFCalculateTheta

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    NUM_ERROR if the testvectors are incompatible
   D*/
/****************************************************************************/

INT FFCalculateThetaAndUpdate( const BLOCKVECTOR *bv_dest,
                               const BLOCKVECTOR *bv_source,
                               const BV_DESC *bvd_dest,
                               const BV_DESC *bvd_source,
                               const BV_DESC_FORMAT *bvdf,
                               INT tv1_comp,
                               INT tv2_comp,
                               GRID *grid)
{
  register VECTOR *vi, *vip1, *end_v;
  register DOUBLE e1i, e2i, e1ip1, e2ip1, a1, a2, det, off_val;
  register MATRIX *m_offdiag;
  INT aux1_comp, aux2_comp, K_comp, T_comp;
  CONNECTION *con;

  aux1_comp = GET_AUX_VEC;
  ASSERT( aux1_comp>-1 );
  aux2_comp = GET_AUX_VEC;
  ASSERT( aux2_comp>-1 );
  K_comp = STIFFMAT_ON_LEVEL_BLOCKWISE(bv_dest);
  T_comp = DECOMPMAT_ON_LEVEL_BLOCKWISE(bv_dest);

  ASSERT( aux1_comp != DUMMY_COMP );
  ASSERT( aux2_comp != DUMMY_COMP );
  ASSERT( tv1_comp != DUMMY_COMP );
  ASSERT( tv2_comp != DUMMY_COMP );
  ASSERT( K_comp != DUMMY_COMP );
  ASSERT( T_comp != DUMMY_COMP );
  ASSERT( K_comp != T_comp );
  ASSERT( aux1_comp != aux2_comp );
  ASSERT( aux1_comp != tv1_comp );
  ASSERT( aux1_comp != tv2_comp );
  ASSERT( aux2_comp != tv1_comp );
  ASSERT( aux2_comp != tv2_comp );
  ASSERT( tv1_comp != tv2_comp );
  /*ASSERT( K_comp == STIFFMAT_ON_LEVEL_BLOCKWISE(bv_source) );*/
  /*ASSERT( T_comp == DECOMPMAT_ON_LEVEL_BLOCKWISE(bv_source) );*/

  /* aux_i := L_(i,i+1) * tv_i+1 */
  dsetBS( bv_source, aux1_comp, 0.0 );
  dsetBS( bv_source, aux2_comp, 0.0 );
  dmatmul_addBS( bv_source, bvd_dest, bvdf, aux1_comp, K_comp, tv1_comp );
  dmatmul_addBS( bv_source, bvd_dest, bvdf, aux2_comp, K_comp, tv2_comp );

  /* aux_i = (T_i)^-1 * aux_i */
#ifdef ModelP
  FFMultWithMInv( bv_source, bvd_source, bvdf, aux1_comp, aux1_comp, NULL, NULL );
  FFMultWithMInv( bv_source, bvd_source, bvdf, aux2_comp, aux2_comp, NULL, NULL );
  if ( BVNUMBER(bv_source) == -100 )       /* bv_source == Lines */
  {
#ifdef FFCOMM
    FFVectorConsistent( (BLOCKVECTOR*)bv_source, aux1_comp );
    FFVectorConsistent( (BLOCKVECTOR*)bv_source, aux2_comp );
#else
    if( l_vector_consistentBS( grid, bvd_source, bvdf, aux1_comp )!=NUM_OK ) REP_ERR_RETURN (1);
    if( l_vector_consistentBS( grid, bvd_source, bvdf, aux2_comp )!=NUM_OK ) REP_ERR_RETURN (1);
#endif
  }
#else
  FFMultWithMInv( bv_source, bvd_source, bvdf, aux1_comp, aux1_comp );
  FFMultWithMInv( bv_source, bvd_source, bvdf, aux2_comp, aux2_comp );
#endif

  /* aux_i+1 := L_(i+1,i) * aux_i */
  dsetBS( bv_dest, aux1_comp, 0.0 );
  dsetBS( bv_dest, aux2_comp, 0.0 );
  dmatmul_addBS( bv_dest, bvd_source, bvdf, aux1_comp, K_comp, aux1_comp );
  dmatmul_addBS( bv_dest, bvd_source, bvdf, aux2_comp, K_comp, aux2_comp );


  /* calculate Theta */
  /* Theta must fulfill the equation Theta * tv = aux = L * T^-1 * L * tv */

  dmatsetBS( bv_dest, bvd_dest, bvdf, T_comp, 0.0 );

  vi = BVFIRSTVECTOR( bv_dest );
  end_v = BVLASTVECTOR( bv_dest );              /* stop before last vector */

  /* prepare for loop */
  e1ip1 = VVALUE( vi, tv1_comp ); e2ip1 = VVALUE( vi, tv2_comp );
  a1 = VVALUE( vi, aux1_comp );   a2 = VVALUE( vi, aux2_comp );

  BLOCK_L_VLOOP( vi, BVFIRSTVECTOR(bv_dest), end_v )
  {
    vip1 = SUCCVC( vi );

    e1i = e1ip1;                                            e2i = e2ip1;
    e1ip1 = VVALUE( vip1, tv1_comp );       e2ip1 = VVALUE( vip1, tv2_comp );

    det = e1i*e2ip1 - e1ip1*e2i;
    if ( fabs(det) < SMALL_D )
    {
      PRINTDEBUG( np, 0, ("FFCalculateThetaAndUpdate: testvectors incompatible\n") );
      printf("tv1\n"); printvBS(bv_dest, tv1_comp);printf("tv2\n"); printvBS(bv_dest, tv2_comp);
      REP_ERR_RETURN(NUM_ERROR);
    }

    MVALUE( VSTART(vi), T_comp ) = MVALUE( VSTART(vi), K_comp ) -
                                   (a1*e2ip1 - a2*e1ip1) / det;

    /*m_offdiag = GetMatrix( vi, vip1 );*/
    /*ASSERT( m_offdiag != NULL );*/
    if ( (m_offdiag = GetMatrix( vi, vip1 )) == NULL )
    {
      if ( (con = CreateExtraConnection( grid, vi, vip1 )) == NULL )
      {
        PrintErrorMessage( 'E', "FFCalculateThetaAndUpdate", "Not enough memory" );
        REP_ERR_RETURN (NUM_OUT_OF_MEM);
      }
      m_offdiag = CMATRIX0( con );
    }
    off_val = (e1i*a2 - e2i*a1) / det;
    if( m_offdiag != NULL )
    {
      MVALUE( m_offdiag, T_comp ) = MVALUE( m_offdiag, K_comp ) - off_val;
      MVALUE( MADJ(m_offdiag), T_comp ) = MVALUE(MADJ(m_offdiag), K_comp ) - off_val;
    }
    else
      ASSERT(FALSE);

    a1 = VVALUE( vip1, aux1_comp ) - off_val*e1i;
    a2 = VVALUE( vip1, aux2_comp ) - off_val*e2i;
  }

  /* end_v points to last vector and its value is e?ip1 */
  MVALUE( VSTART(end_v), T_comp ) = MVALUE( VSTART(end_v), K_comp ) -
                                    (e1ip1*a1 + e2ip1*a2) / (e1ip1*e1ip1 + e2ip1 * e2ip1);

  FREE_AUX_VEC(aux2_comp);
  FREE_AUX_VEC(aux1_comp);

  return(NUM_OK);
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
   If necessary new connections are allocated as 'extra connections'.

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

  return(NUM_OK);
}


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
   FFDecomp, TFFCalculateTheta, FFConstructTestvector_loc, FFConstructTestvector, FFMultWithM, FFMultWithMInv

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
               INT tv_comp,
               GRID *grid )
{
  register BLOCKVECTOR *bv_i, *bv_im1, *bv_end;
  register BV_DESC *bvd_i, *bvd_im1, *bvd_temp;
  BV_DESC bvd1, bvd2;
  INT K_comp, FF_comp;
#ifdef THETA_ANA
  DOUBLE lambda, li, a, bx, by, bz;
  MATRIX *m;
#endif

  ASSERT( !BV_IS_EMPTY(bv) );

  K_comp = STIFFMAT_ON_LEVEL(bv);
  FF_comp = DECOMPMAT_ON_LEVEL(bv);

  if ( BV_IS_LEAF_BV(bv) )
  {
    dmatcopyBS( bv, bvd, bvdf, FF_comp, K_comp );
    return LUDecomposeDiagBS( bv, bvd, bvdf, FF_comp, grid );
  }

  if ( BV_IS_DIAG_BV(bv) )
  {
    bvd1 = *bvd;
    bv_end = BVDOWNBVEND(bv);
    for ( bv_i = BVDOWNBV( bv ); bv_i != bv_end; bv_i = BVSUCC( bv_i ) )
      if( !BV_IS_EMPTY(bv_i) )
      {
        BVD_PUSH_ENTRY( &bvd1, BVNUMBER(bv_i), bvdf );
        TFFDecomp( wavenr, wavenr3D, bv_i, &bvd1, bvdf, tv_comp, grid );
        BVD_DISCARD_LAST_ENTRY(&bvd1);
      }

    return(NUM_OK);
  }

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

  /* initialize the BVDs */
  bvd1 = bvd2 = *bvd;                   /* copy of the BVD */
  bvd_i = &bvd1;
  bvd_im1 = &bvd2;

  /* the loop calculates (T_i-1)^-1 and Theta(i,j) for i=1..n-1 */
  /* T_0 := D_0 */

  bv_end = BVDOWNBVEND(bv);

  /* set up first block */
  bv_im1 = BVDOWNBV(bv);
  while( BV_IS_EMPTY( bv_im1 ) && (bv_im1 != bv_end) )          /* search first nonempty bv */
    bv_im1 = BVSUCC( bv_im1 );
  ASSERT( bv_im1 != bv_end );
  BVD_PUSH_ENTRY( bvd_im1, BVNUMBER(bv_im1), bvdf );

  /* set up second block */
  bv_i = BVSUCC( bv_im1 );
  while( (bv_i != bv_end) && BV_IS_EMPTY( bv_i ) )              /* search first nonempty bv */
    bv_i = BVSUCC( bv_i );
  if ( bv_i != bv_end )
    BVD_PUSH_ENTRY( bvd_i, BVNUMBER(bv_i), bvdf );
  /* else: bv_i and bvd_i are never used;
     thus the content of them have not to be updated */

  dmatcopyBS( bv_im1, bvd_im1, bvdf, FF_comp, K_comp );
  for ( ; bv_i != bv_end; )
  {
    /* T_i-1 decompose */
    TFFDecomp( wavenr, wavenr3D, bv_im1, bvd_im1, bvdf, tv_comp, grid );

    /* calculate Theta_(i-1,i) and Theta_(i,i-1);
       result on off-digonal blocks of FF_comp */
    FFConstructTestvector_loc( bv_i, tv_comp, wavenr, wavenr3D );

#ifdef THETA_ANA
    dmatsetBS( bv_i, bvd_im1, bvdf, FF_comp, 1.0/li );
    dmatsetBS( bv_im1, bvd_i, bvdf, FF_comp, 1.0/li );
    /*printf("theta = %g\n", 1.0 / li );*/
    li = lambda - 1.0 / li;
    /*
       #elif defined FF_PARALLEL_SIMULATION
                    BoxTFFCalculateTheta( bv_i, bv_im1, bvd_i, bvd_im1, bvdf, tv_comp, grid );
     */
#else
    TFFCalculateTheta( bv_i, bv_im1, bvd_i, bvd_im1, bvdf, tv_comp );
#endif

    /* T_i := D_i + Theta_(i,i-1)*T_i-1*Theta_(i-1,i) -
                      Theta_(i,i-1)*L_(i-1,i) - L_(i,i-1)*Theta_(i-1,i) */
    dmatcopyBS( bv_i, bvd_i, bvdf, FF_comp, K_comp );
    TFFUpdateDiagBlock( bv_i, bvd_i, bvd_im1, bvdf, FF_comp, K_comp, FF_comp, grid );

    /* update BVDs for the next loop */
    bv_im1 = bv_i;
    SWAP( bvd_i, bvd_im1, bvd_temp );
    bv_i = BVSUCC( bv_i );
    while( (bv_i != bv_end) && BV_IS_EMPTY( bv_i ) )                    /* search first nonempty bv */
      bv_i = BVSUCC( bv_i );
    if( bv_i != bv_end )
    {
      BVD_DISCARD_LAST_ENTRY( bvd_i );
      BVD_PUSH_ENTRY( bvd_i, BVNUMBER(bv_i), bvdf );
    }
    /* else: bv_i and bvd_i are never used;
       thus the content of them have not to be updated */
  }
  /* now bv_im1 and bvd_im1 points to the last block */

  /* calculate the last (T_n)^-1 */
  TFFDecomp( wavenr, wavenr3D, bv_im1, bvd_im1, bvdf, tv_comp, grid );


#ifdef CHECK_CALCULATION
  {
    /* notice: for the first block exists no testvector. But it is still 0
                   thus it do not disturb the check */
    DOUBLE norm;

    /* to get the right global testvector for the tests */
    if( K_comp == 0 )
      FFConstructTestvector( bv, tv_comp, wavenr, wavenr3D );

    /* aux2 := M * tv */
    FFMultWithM( bv, bvd, bvdf, aux2_COMP, FF_comp, K_comp, LU_comp, tv_comp, aux_comp, auxsub_comp, FF_comp );

    /* aux := K * tv */
    dsetBS( bv, aux_comp, 0.0 );
    dmatmul_addBS( bv, bvd, bvdf, aux_comp, K_comp, tv_comp );

    /* aux2 -= aux */
    dsubBS( bv, aux2_COMP, aux_comp );

    /* norm of aux2 */
    dnrm2BS( bv, aux2_COMP, &norm );

    if ( norm > FFaccuracy )
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
  return(NUM_OK);
}


/****************************************************************************/
/*D
   FFDecomp - calculate the tangential frequency filtering decomposition of a matrix

   SYNOPSIS:
   INT FFDecomp( DOUBLE wavenr,
                                  DOUBLE wavenr3D,
                                  const BLOCKVECTOR *bv,
                                  const BV_DESC *bvd,
                                  const BV_DESC_FORMAT *bvdf,
                                  INT LU_comp,
                                  INT FF_comp,
                                  INT K_comp,
                                  INT tv1_comp,
                                  INT tv2_comp,
                                  INT aux1_comp,
                                  INT aux2_comp,
                                  INT auxsub1_comp,
                                  INT auxsub2_comp,
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
   .  tv1_comp - position of the 1. testvector-component in the VECTOR-data
   .  tv2_comp - position of the 2. testvector-component in the VECTOR-data
   .  aux1_comp - position of the 1. auxiliary-component in the VECTOR-data
   .  aux2_comp - position of the 2. auxiliary-component in the VECTOR-data
   .  auxsub1_comp - position of the auxiliary-component for subproblems (only 3D)
   .  auxsub2_comp - position of the auxiliary-component for subproblems (only 3D)
   .  FFsub_comp - position of the Theta and filtered diagonal blocks of the subproblem-matrix (only 3D)
   .  grid - grid where the matrixes lives on (for allocating extra connections)

   DESCRIPTION:
   This function calculates the (classical) frequency filtered decomposition
   of the given matrix 'K'. It expects a blockvector structure according to a
   linewise (2D) resp. plane-/linewise (3D) decomposition of the domain as
   it is constructed by 'CreateBVStripe2D'/'CreateBVStripe3D'. The result is
   a matrix 'M' in the form 'M = ( L + T ) * T^-1 * ( T + U )' where the
   lower ('L') and upper ('U') off-diagonal blocks are the same as in
   the given matrix ('K_comp'), 'T' is the frequency filtered diagonal
   blockmatrix ('FF_comp'). For the auxiliary matrix Theta the diagonal
   blocks of 'FF_comp' are used temporarely.

   In the 3D case the frequency filtered diagonal blocks are again decomposed
   by FF. The matrix to be decomposed is in the 'FF_comp' and the frequency
   filtered diagonal blocks are stored in 'FFsub_comp' (their LU decomposition
   finally in 'LU_comp'). 'auxsub_comp' and 'FFsub_comp' are used only in
   the 3D case; in 2D you can use an arbitrary number.

   The distinction whether 'bv' represents a plane (2D) or a cube (3D) is made
   dynamically by the blockvector structure itself: if 'bv' has only one
   further level of blockvectors (representing grid lines!) it is
   considered as a plane, otherwise as a cube.

   The sonlist of 'bv' must contain at least one non-empty blockvector.

   SEE ALSO:
   TFFDecomp, FFCalculateTheta, FFConstructTestvector_loc, FFConstructTestvector, FFMultWithM, FFMultWithMInv

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    error code from 'LUDecomposeDiagBS'
   D*/
/****************************************************************************/

INT FFDecomp( DOUBLE wavenr,
              DOUBLE wavenr3D,
              const BLOCKVECTOR *bv,
              const BV_DESC *bvd,
              const BV_DESC_FORMAT *bvdf,
              INT tv1_comp,
              INT tv2_comp,
              GRID *grid )
{
  register BLOCKVECTOR *bv_i, *bv_im1, *bv_end;
  register BV_DESC *bvd_i, *bvd_im1, *bvd_temp;
  BV_DESC bvd1, bvd2;
  INT K_comp, FF_comp;

  ASSERT( !BV_IS_EMPTY(bv) );

  K_comp = STIFFMAT_ON_LEVEL(bv);
  FF_comp = DECOMPMAT_ON_LEVEL(bv);

  ASSERT( FF_comp != DUMMY_COMP );
  ASSERT( K_comp != DUMMY_COMP );
  ASSERT( tv1_comp != DUMMY_COMP );
  ASSERT( tv2_comp != DUMMY_COMP );
  ASSERT( tv1_comp != tv2_comp );

  if ( BV_IS_LEAF_BV(bv) )
  {
    dmatcopyBS( bv, bvd, bvdf, FF_comp, K_comp );
    return LUDecomposeDiagBS( bv, bvd, bvdf, FF_comp, grid );
  }

  if ( BV_IS_DIAG_BV(bv) )
  {
    bvd1 = *bvd;
    bv_end = BVDOWNBVEND(bv);
    for ( bv_i = BVDOWNBV( bv ); bv_i != bv_end; bv_i = BVSUCC( bv_i ) )
      if( !BV_IS_EMPTY(bv_i) )
      {
        BVD_PUSH_ENTRY( &bvd1, BVNUMBER(bv_i), bvdf );
        FFDecomp( wavenr, wavenr3D, bv_i, &bvd1, bvdf, tv1_comp, tv2_comp, grid );
        BVD_DISCARD_LAST_ENTRY(&bvd1);
      }

    return(NUM_OK);
  }

  /* initialize the BVDs */
  bvd1 = bvd2 = *bvd;                   /* copy of the BVD */
  bvd_i = &bvd1;
  bvd_im1 = &bvd2;

  /* the loop calculates (T_i-1)^-1 and Theta(i,j) for i=1..n-1 */
  /* T_0 := D_0 */

  bv_end = BVDOWNBVEND(bv);

  /* set up first block */
  bv_im1 = BVDOWNBV(bv);
  while( BV_IS_EMPTY( bv_im1 ) && (bv_im1 != bv_end) )          /* search first nonempty bv */
    bv_im1 = BVSUCC( bv_im1 );
  ASSERT( bv_im1 != bv_end );
  BVD_PUSH_ENTRY( bvd_im1, BVNUMBER(bv_im1), bvdf );

  /* set up second block */
  bv_i = BVSUCC( bv_im1 );
  while( (bv_i != bv_end) && BV_IS_EMPTY( bv_i ) )              /* search first nonempty bv */
    bv_i = BVSUCC( bv_i );
  if ( bv_i != bv_end )
    BVD_PUSH_ENTRY( bvd_i, BVNUMBER(bv_i), bvdf );
  /* else: bv_i and bvd_i are never used;
     thus the content of them have not to be updated */

  dmatcopyBS( bv_im1, bvd_im1, bvdf, FF_comp, K_comp );
  for ( ; bv_i != bv_end; )
  {
    /* T_i-1 decompose */
    FFDecomp( wavenr, wavenr3D, bv_im1, bvd_im1, bvdf, tv1_comp, tv2_comp, grid );

    /* calculate Theta_(i,i);
       result on digonal blocks of FF_comp */
    if ( BVNUMBER(bv_i)== -101 )
    {                   /* cross point sysytem aus ptff.c */
      FFConstructTestvector_loc( bv_i, tv1_comp, 1.0, 1.0 );
      FFConstructTestvector_loc( bv_i, tv2_comp, 2.0, 2.0 );
      printf( "special crosspoint tv\n");
      /*printvBS( bv_i, tv1_comp ); printvBS( bv_i, tv2_comp );*/
    }
    else
    {
      FFConstructTestvector_loc( bv_i, tv1_comp, wavenr, wavenr3D );
      FFConstructTestvector_loc( bv_i, tv2_comp, wavenr+1.0, wavenr3D );
    }

#if (defined FF_ModelP) || (defined FF_PARALLEL_SIMULATION)
    if ( BVNUMBER(bv_i) == -100 )               /* lines */
    {
      /* construct FF filtered approximation of the leaf blocks of the schur complement */
      /* the determination of the "lines" block is only a quick hack!!! */

      ConstructSchurFFApprox( bv_im1, bv_i, bvd_im1, bvd_i, bvdf, tv1_comp, tv2_comp, grid );

                        #ifdef ModelP
      /* make Schur complement matrix for lines consistent */
      /*printf(PFMT" vor\n",me);printmBS(bv_i,bv_i,FF_comp);*/
#ifdef FFCOMM
      FFTridiagMatConsistent( bv_i, FF_comp );
#else
      ASSERT(grid!=NULL);
      if( l_matrix_consistent( grid, DECOMP_MATDATA_DESC_ON_LEVEL(bv), MAT_CONS )!=NUM_OK ) REP_ERR_RETURN(NUM_ERROR);
#endif
      /*printf(PFMT" nach\n",me);printmBS(bv_i,bv_i,FF_comp);*/
                        #endif
    }
    else
    {
      FFCalculateThetaAndUpdate( bv_i, bv_im1, bvd_i, bvd_im1, bvdf, tv1_comp, tv2_comp, grid );

                        #ifdef ModelP
      if ( BVNUMBER(bv_i) == -101 )                     /* crosspoints */
      {
        /* make Schur complement matrix for cross points consistent */
#ifdef FFCOMM
        FFTridiagMatConsistent( bv_i, FF_comp );
#else
        ASSERT(grid!=NULL);
        if( l_matrix_consistent( grid, DECOMP_MATDATA_DESC_ON_LEVEL(bv), MAT_CONS )!=NUM_OK ) REP_ERR_RETURN(NUM_ERROR);
#endif
      }
                        #endif
    }
#else
    FFCalculateThetaAndUpdate( bv_i, bv_im1, bvd_i, bvd_im1, bvdf, tv1_comp, tv2_comp, grid );
#endif
    /* update BVDs for the next loop */
    bv_im1 = bv_i;
    SWAP( bvd_i, bvd_im1, bvd_temp );
    bv_i = BVSUCC( bv_i );
    while( (bv_i != bv_end) && BV_IS_EMPTY( bv_i ) )                    /* search first nonempty bv */
      bv_i = BVSUCC( bv_i );
    if( bv_i != bv_end )
    {
      BVD_DISCARD_LAST_ENTRY( bvd_i );
      BVD_PUSH_ENTRY( bvd_i, BVNUMBER(bv_i), bvdf );
    }
    /* else: bv_i and bvd_i are never used;
       thus the content of them have not to be updated */
  }
  /* now bv_im1 and bvd_im1 points to the last block */

  /* calculate the last (T_n)^-1 */
  FFDecomp( wavenr, wavenr3D, bv_im1, bvd_im1, bvdf, tv1_comp, tv2_comp, grid );

  /*	printvBS(bv, tv1_comp);printvBS(bv, tv2_comp);*/
  /*printf("K_comp\n"); printmBS(bv,bv,K_comp);printf("FF_comp\n"); printmBS(bv,bv,FF_comp); printf("\n");*/
  return(NUM_OK);
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
   TFFCalculateTheta, FFConstructTestvector_loc, FFConstructTestvector, FFMultWithM, FFMultWithMInv, TFFDecomp, CreateBVStripe2D, CreateBVStripe3D

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   D*/
/****************************************************************************/

INT  TFFSolve( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT K_comp, INT u_comp, INT f_comp, INT cor_comp, INT FF_comp, INT LU_comp, INT tv_comp, INT aux_comp, INT auxsub_comp, INT FFsub_comp, DOUBLE meshwidth, DOUBLE eps, GRID *grid )
/* not necessary for the new np */
{
  DOUBLE old_norm, new_norm, start_norm, step_norm, final_acc;
  INT i, j=0, it, nr_TFFs;

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
    /*FFConstructTestvector( bv, tv_comp, (double)(1<<i), (double)(1<<j) );*/

    /* Calculates the TFF decomposition of M = (L + T) * T^-1 * (L^T + T) */
    TFFDecomp( (double)(1<<i), (double)(1<<j), bv, bvd, bvdf, tv_comp, grid );

    /* cor := M^-1 * f */
    dcopyBS( bv, cor_comp, f_comp );
#ifdef ModelP
    FFMultWithMInv( bv, bvd, bvdf, cor_comp, cor_comp, NULL, NULL );
#else
    FFMultWithMInv( bv, bvd, bvdf, cor_comp, cor_comp );
#endif

#ifdef CHECK_CALCULATION
    {
      DOUBLE norm;

      /* check invertation of testvector tv == M^-1 * K * tv;
              must always be true since it is equivalent to
              S * tv == 0, i.e. solve exact on the tv-subspace*/

      /* aux2 := K * tv */
      dsetBS( bv, aux2_COMP, 0.0 );
      dmatmul_addBS( bv, bvd, bvdf, aux2_COMP, K_comp, tv_comp );

      /* aux2 := M^-1 * aux2 */
#ifdef ModelP
      FFMultWithMInv( bv, bvd, bvdf, aux2_COMP, K_comp, LU_comp, aux2_COMP, aux3_COMP, auxsub_comp, FF_comp, NULL, NULL );
#else
      FFMultWithMInv( bv, bvd, bvdf, aux2_COMP, K_comp, LU_comp, aux2_COMP, aux3_COMP, auxsub_comp, FF_comp );
#endif

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      dnrm2BS( bv, aux2_COMP, &norm );

      if ( norm > FFaccuracy )
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
      FFMultWithM( bv, bvd, bvdf, aux2_COMP, FF_comp, K_comp, LU_comp, tv_comp, aux3_COMP, auxsub_comp, FF_comp );

      /* aux2 := M^-1 * aux2 */
      dcopyBS( bv, aux_comp, tv_comp );
#ifdef ModelP
      FFMultWithMInv( bv, bvd, bvdf, aux2_COMP, K_comp, LU_comp, aux2_COMP, aux3_COMP, auxsub_comp, FF_comp, NULL, NULL );
#else
      FFMultWithMInv( bv, bvd, bvdf, aux2_COMP, K_comp, LU_comp, aux2_COMP, aux3_COMP, auxsub_comp, FF_comp );
#endif

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      dnrm2BS( bv, aux2_COMP, &norm );

      if ( norm > FFaccuracy )
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
#ifdef ModelP
      FFMultWithMInv( bv, bvd, bvdf, aux_comp, K_comp, LU_comp, aux_comp, aux3_COMP, auxsub_comp, FF_comp, NULL, NULL );
#else
      FFMultWithMInv( bv, bvd, bvdf, aux_comp, K_comp, LU_comp, aux_comp, aux3_COMP, auxsub_comp, FF_comp );
#endif

      /* aux2 := M * aux */
      FFMultWithM( bv, bvd, bvdf, aux2_COMP, FF_comp, K_comp, LU_comp, aux_comp, aux3_COMP, auxsub_comp, FF_comp );

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      dnrm2BS( bv, aux2_COMP, &norm );

      if ( norm > FFaccuracy )
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
#ifdef ModelP
      FFMultWithMInv( bv, bvd, bvdf, aux3_COMP, K_comp, LU_comp, aux3_COMP, aux2_COMP, auxsub_comp, FF_comp, NULL, NULL );
#else
      FFMultWithMInv( bv, bvd, bvdf, aux3_COMP, K_comp, LU_comp, aux3_COMP, aux2_COMP, auxsub_comp, FF_comp );
#endif

      /* aux2 := K * aux3 */
      dsetBS( bv, aux2_COMP, 0.0 );
      dmatmul_addBS( bv, bvd, bvdf, aux2_COMP, K_comp, aux3_COMP );

      /* aux2 -= tv */
      dsubBS( bv, aux2_COMP, tv_comp );

      /* norm of aux2 */
      dnrm2BS( bv, aux2_COMP, &norm );

      if ( norm > FFaccuracy )
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
        dmatmul_addBS( bv, bvd, bvdf, aux2_COMP, K_comp, cor_comp );

        /* aux2 -= aux */
        dsubBS( bv, aux2_COMP, aux_comp );

        /* norm of aux2 */
        dnrm2BS( bv, aux2_COMP, &norm );

        if ( norm > FFaccuracy )
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

return(NUM_OK);
}

#endif /* __BLOCK_VECTOR_DESC__ */
