// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ff_gen.c                                                                                                      */
/*																			*/
/* Purpose:   general frequency filtering decompostion routines             */
/*                                                                                                                                              */
/* Author:	  Christian Wrobel                                                                              */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
/*																			*/
/* History:   14.11.95 begin, ug version 3.1								*/
/*			  11.03.96 restructuring										*/
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
#include <string.h>

#include "compiler.h"
#include "gm.h"          /* for data structure               */
#include "devices.h"     /* for UserWrite, PrintErrorMessage */
#include "commands.h"    /* for GetCurrentMultigrid              */
#include "cmdint.h"      /* for CreateCommand                */
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
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
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

/* mute level for FFs */
INT mute_level  = 0;

/* value below them a division is refused calculating a testvector */
DOUBLE FFsmallTV = 1e-3;

/* ratio for a jump to be detected */
DOUBLE FFmuchBigger = 100.0;

/* value below them a number is considered as 0.0 */
DOUBLE FFEPS = 1e-16;

/* value below them an approximation error is considered as ok */
DOUBLE FFaccuracy = 1e-10;

/* global array to hold the matrix hierarchy */
INT FF_Mats[FF_MAX_MATS];
MATDATA_DESC *FF_MATDATA_DESC_ARRAY[FF_MAX_MATS];

/* global array to hold the auxiliary vectors */
INT FF_Vecs[FF_MAX_VECS];
INT TOS_FF_Vecs = 0;


/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
INT aux2_COMP = DUMMY_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
INT aux3_COMP = DUMMY_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
INT aux4_COMP = DUMMY_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
INT aux5_COMP = DUMMY_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
INT aux6_COMP = DUMMY_COMP;


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* special blas routines for BLOCKVECTOR									*/
/*																			*/
/****************************************************************************/

INT     storeVectorBS( BLOCKVECTOR *bv, INT x_comp, GRID *grid )
/* store the given vector in an allocated buffer in the blockvector */
{
  register DOUBLE *mem;
  register VECTOR *v, *end_v;

  if ( (mem = BACKUP_MEM(bv)) == NULL )
  {
    if ( grid == NULL )
    {
      PrintErrorMessage( 'E', "storeVectorBS", "No memory allocated in blockvector" );
      REP_ERR_RETURN(GM_OUT_OF_MEM);

    }
    if ( (mem = GetMem( MGHEAP(MYMG(grid)), BVNUMBEROFVECTORS(bv)*sizeof(DOUBLE), FROM_BOTTOM )) == NULL )
    {
      PrintErrorMessage( 'E', "storeVectorBS", "Not enough memory to store the vector" );
      REP_ERR_RETURN(GM_OUT_OF_MEM);
    }
    BVUSERDATA( bv ) = mem;
  }

  end_v = BVENDVECTOR( bv );
  BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )
  *mem++ = VVALUE( v, x_comp );

  /* the number of copied elements must be equal the number of vectors
          stored in the blockvector for them the memory is allocated */
  ASSERT( ( mem - BACKUP_MEM(bv) ) == BVNUMBEROFVECTORS( bv ) );

  return NUM_OK;
}


INT     restoreVectorBS( BLOCKVECTOR *bv, INT x_comp )
/* restore the given vector from an allocated buffer in the blockvector */
{
  register DOUBLE *mem;
  register VECTOR *v, *end_v;

  mem = BACKUP_MEM( bv );

  ASSERT( mem != NULL );

  end_v = BVENDVECTOR( bv );
  BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )
  VVALUE( v, x_comp ) = *mem++;

  return NUM_OK;
}


#ifdef NIE_DAA    /* nur fuer veraltete Datenstruktur geeignet */
static INT HasDirichletNeighbour( const VECTOR *v )
{
  register NODE *my_node, *nb_node;
  register LINK *link;
  register VERTEX *nb_vertex;

  ASSERT(FALSE);       /* neues boundary konzept einbauen */
  REP_ERR_RETURN(1);

  register VSEGMENT *vseg;
  register BndCondProcPtr bnd_cond;
  register BNDSEGDESC *bnd_desc;

  DOUBLE value;
  INT type;

  /* only node vectors are considered */
  ASSERT( VOTYPE(v) == NODEVEC );

  my_node = (NODE*)VOBJECT( v );

  for ( link = START( my_node ); link != NULL; link = NEXT( link ) )
  {
    nb_node = NBNODE( link );
    if ( nb_node == my_node )
      continue;

    nb_vertex = MYVERTEX( nb_node );

    if ( OBJT(nb_vertex) != BVOBJ )
      continue;                         /* not a boundary vertex, since no dirichlet neighbour */

    for ( vseg = VSEG(nb_vertex); vseg != NULL; vseg = NEXTSEG(vseg) )
    {
      bnd_desc = BSEG(vseg);
      bnd_cond = (BndCondProcPtr)BNDCONDFUNC(bnd_desc);

      (*bnd_cond)( SEGDATA(bnd_desc), CONDDATA(bnd_desc), PVECT(vseg),
                   &value, &type );

      if ( type == DIRICHLET )
        return TRUE;
    }
  }

  return FALSE;
}
#endif /* NIE_DAA */

void printv( INT x_nr )
/* for calling from a debugger */
{
  register VECTOR *v;
  DOUBLE pos[DIM];
  GRID *g;

  g = GRID_ON_LEVEL(GetCurrentMultigrid(), CURRENTLEVEL(GetCurrentMultigrid()));

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    VectorPosition(v,pos);
#ifdef ModelP
    printf(PFMT,me);
#endif
    printf("x=%5.2f y=%5.2f ",pos[0],pos[1]);
#ifdef __THREEDIM__
    printf("z=%5.2f ",pos[2]);
#endif
    printf("  index = %d  ", VINDEX( v ) );
    printf("u[%d]=%15.8lf ",x_nr,VVALUE(v,x_nr));
    /*printf("   cl %d %d sk ",VCLASS(v),VNCLASS(v));*/
    /*for (j=0; j<ncomp; j++)
            printf("%d ",((VECSKIP(v) & (1<<j))!=0));*/
    printf("\n");
  }
  return;
}

void printvgrid( GRID *g, INT x_nr )
/* for calling from a debugger */
{
  register VECTOR *v;
  DOUBLE pos[DIM];

  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
    VectorPosition(v,pos);
#ifdef ModelP
    printf(PFMT,me);
#endif
    printf("x=%5.2f y=%5.2f ",pos[0],pos[1]);
#ifdef __THREEDIM__
    printf("z=%5.2f ",pos[2]);
#endif
    printf("  index = %d  ", VINDEX( v ) );
    printf("u[%d]=%15.8lf ",x_nr,VVALUE(v,x_nr));
    /*printf("   cl %d %d sk ",VCLASS(v),VNCLASS(v));*/
    /*for (j=0; j<ncomp; j++)
            printf("%d ",((VECSKIP(v) & (1<<j))!=0));*/
    printf("\n");
  }
  return;
}

void printvBS( const BLOCKVECTOR *bv, INT x_nr )
/* for calling from a debugger */
{
  register VECTOR *v;
  DOUBLE pos[DIM];

  for (v=BVFIRSTVECTOR(bv); v!= BVENDVECTOR(bv); v=SUCCVC(v))
  {
    VectorPosition(v,pos);
#ifdef ModelP
    printf(PFMT,me);
#endif
    printf("x=%5.2f y=%5.2f ",pos[0],pos[1]);
#ifdef __THREEDIM__
    printf("z=%5.2f ",pos[2]);
#endif
    printf("  index = %d  ", VINDEX( v ) );
    printf("u[%d]=%15g ",x_nr,VVALUE(v,x_nr));
    /*printf("   cl %d %d sk ",VCLASS(v),VNCLASS(v));*/
    /*for (j=0; j<ncomp; j++)
            printf("%d ",((VECSKIP(v) & (1<<j))!=0));*/
    printf("\n");
  }
  return;
}

void printm( INT m_nr )
/* for calling from a debugger */
{
  register VECTOR *v, *w;
  register MATRIX *m;
  GRID *g;

  g = GRID_ON_LEVEL(GetCurrentMultigrid(), CURRENTLEVEL(GetCurrentMultigrid()));

#ifdef ModelP
  printf(PFMT,me);
#endif
  printf("comp (%d)\n",m_nr);
  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    for ( w = FIRSTVECTOR( g ); w != NULL; w = SUCCVC( w ) )
    {
      for (m=VSTART(v); m!=NULL; m = MNEXT(m))
        if ( MDEST( m ) == w )
        {
          printf("%5.2lf",(double)MVALUE(m,m_nr));
          break;
        }
      if ( m == NULL )
        printf("     ");
    }
    printf("\n");
  }

  return;
}

void printmgrid( GRID *g, INT m_nr )
/* for calling from a debugger */
{
  register VECTOR *v, *w;
  register MATRIX *m;

#ifdef ModelP
  printf(PFMT,me);
#endif
  printf("comp (%d)\n",m_nr);
  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    for ( w = FIRSTVECTOR( g ); w != NULL; w = SUCCVC( w ) )
    {
      for (m=VSTART(v); m!=NULL; m = MNEXT(m))
        if ( MDEST( m ) == w )
        {
          printf("%5.2lf",(double)MVALUE(m,m_nr));
          break;
        }
      if ( m == NULL )
        printf("     ");
    }
    printf("\n");
  }

  return;
}

void printmMG( MULTIGRID *theMG, INT m_nr )
/* for calling from a debugger */
{
  register VECTOR *v, *w;
  register MATRIX *m;
  GRID *g; int i;

  for ( i=0; i <= TOPLEVEL(theMG); i++ )
  {
    g = GRID_ON_LEVEL( theMG, i);
    printf("comp (%d)\n",m_nr);
    for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
    {
#ifdef ModelP
      printf(PFMT,me);
#endif
      for ( w = FIRSTVECTOR( g ); w != NULL; w = SUCCVC( w ) )
      {
        for (m=VSTART(v); m!=NULL; m = MNEXT(m))
          if ( MDEST( m ) == w )
          {
            printf("%5.2lf",(double)MVALUE(m,m_nr));
            break;
          }
        if ( m == NULL )
          printf("     ");
      }
      printf("\n");
    }
  }

  return;
}

void printmBS( const BLOCKVECTOR *bv_row, const BLOCKVECTOR *bv_col, INT m_nr )
/* for calling from a debugger */
{
  register VECTOR *v, *w;
  register MATRIX *m;

#ifdef ModelP
  printf(PFMT,me);
#endif
  printf("comp (%d)\n",m_nr);

  if ( BVNUMBEROFVECTORS( bv_row ) == 0 || BVNUMBEROFVECTORS( bv_col ) == 0 )
  {
    printf( "empty\n" );
    return;
  }

  for (v=BVFIRSTVECTOR(bv_row); v!= BVENDVECTOR(bv_row); v=SUCCVC(v))
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    for ( w = BVFIRSTVECTOR(bv_col); w !=BVENDVECTOR(bv_col); w = SUCCVC( w ) )
    {
      for (m=VSTART(v); m!=NULL; m = MNEXT(m))
        if ( MDEST( m ) == w )
        {
          printf("%7.4lf",(double)MVALUE(m,m_nr));
          break;
        }
      if ( m == NULL )
        printf("       ");
    }
    printf("\n");
  }

  return;
}


void printPatternBS( const BLOCKVECTOR *bv_row, const BLOCKVECTOR *bv_col, INT m_nr )
/* for calling from a debugger */
{
  register VECTOR *v, *w;
  register MATRIX *m;

#ifdef ModelP
  printf(PFMT,me);
#endif
  printf("comp (%d)\n",m_nr);

  if ( BVNUMBEROFVECTORS( bv_row ) == 0 )
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    printf( "empty\n" );
    return;
  }

  for (v=BVFIRSTVECTOR(bv_row); v!= BVENDVECTOR(bv_row); v=SUCCVC(v))
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    for ( w = BVFIRSTVECTOR(bv_col); w !=BVENDVECTOR(bv_col); w = SUCCVC( w ) )
    {
      for (m=VSTART(v); m!=NULL; m = MNEXT(m))
        if ( MDEST( m ) == w )
        {
          if ( MVALUE(m,m_nr) == 0.0 )
            printf(".");
          else
            printf("*");
          break;
        }
      if ( m == NULL )
        printf(" ");
    }
    printf("\n");
  }
  printf("\f");

  return;
}

#ifdef __BLOCK_VECTOR_DESC__

static void printBVrec( BLOCKVECTOR *bv_first, char *indent, const BV_DESC *bvd_parent, const BV_DESC_FORMAT *bvdf )
{
  register VECTOR *v;
  BLOCKVECTOR *bv;
  char indent_rec[200];
  BV_DESC bvd_bv;

  if ( bvdf != NULL )
  {
    bvd_bv = *bvd_parent;
    BVD_PUSH_ENTRY( &bvd_bv, 0, bvdf );
  }

  strcpy( indent_rec, indent );
  strcat( indent_rec, "    " );

  for ( bv = bv_first; bv != NULL; bv = BVSUCC( bv ) )
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    printf( "%s Nr. %d ", indent, BVNUMBER(bv) );
    if ( BVNUMBEROFVECTORS(bv) != 0 )
    {
      printf( "number of vectors %2d ", BVNUMBEROFVECTORS(bv) );
      printf( "first vector %3d ", VINDEX(BVFIRSTVECTOR(bv)) );
      printf( "last vector %3d ", VINDEX(BVLASTVECTOR(bv)) );
      printf( "level %2d", BVLEVEL(bv) );
      printf( " %s", BVORIENTATION(bv)==BVNOORIENTATION ? "(N)" : (BVORIENTATION(bv)==BVVERTICAL ? "(V)" : (BVORIENTATION(bv)==BVHORIZONTAL ? "(H)" : "")));
      printf( "%s\n", BV_IS_DIAG_BV(bv) ? "(D)" : "");

      if ( bvdf != NULL )
      {
        BVD_DISCARD_LAST_ENTRY( &bvd_bv );
        BVD_PUSH_ENTRY( &bvd_bv, BVNUMBER(bv), bvdf );
        for ( v = BVFIRSTVECTOR( bv ); v != BVENDVECTOR( bv ); v = SUCCVC( v ) )
        {
          /*printf("%d %x %d\n", VINDEX(v), VBVD(v).entry, VBVD(v).current );*/
          if ( !VMATCH( v, &bvd_bv, bvdf ) )
            printf( "%s     vector %d doesn't match the blockvector\n", indent, VINDEX(v) );
        }
      }
    }
    else
    {
      printf( "No vectors" );
      if ( BVFIRSTVECTOR( bv ) != NULL || BVLASTVECTOR( bv ) != NULL )
        printf( " but the vector pointers are set ????????" );
      printf( "\n" );
    }

    if ( !BV_IS_LEAF_BV( bv ) )
      printBVrec( BVDOWNBV(bv), indent_rec, &bvd_bv, bvdf );
  }

  return;
}

void printBV( const BV_DESC_FORMAT *bvdf )
{
  BLOCKVECTOR *bv;
  BV_DESC bvd;

  bv = GFIRSTBV( GRID_ON_LEVEL(GetCurrentMultigrid(), CURRENTLEVEL(GetCurrentMultigrid())) );

  if ( bv == NULL )
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    printf( "No blockvectors\n" );
    return;
  }

  BVD_INIT( &bvd );
  printBVrec( bv, "", &bvd, bvdf );
  return;
}


void printBVgrid( GRID *grid, const BV_DESC_FORMAT *bvdf )
{
  BLOCKVECTOR *bv;
  BV_DESC bvd;

  bv = GFIRSTBV( grid );

  if ( bv == NULL )
  {
#ifdef ModelP
    printf(PFMT,me);
#endif
    printf( "No blockvectors\n" );
    return;
  }

  BVD_INIT( &bvd );
  printBVrec( bv, "", &bvd, bvdf );
  return;
}


DOUBLE FFMeshwidthOfGrid( GRID *grid )
/* determine the meshwidth for a regular mesh */
{
  VERTEX *vertex1, *vertex2;
  LINK *link;
  NODE *nb;
  DOUBLE meshwidth;

  /* determine the meshwidth for a regular mesh */
  vertex1 = MYVERTEX( FIRSTNODE( grid ) );

  link = START( FIRSTNODE( grid ) );
  do
  {
    ASSERT( link != NULL );
    nb = NBNODE( link );
    vertex2 = MYVERTEX( nb );
    link = NEXT( link );
  }       /* while exactly in one of the x and y coordinate is a difference */
  while ( ((fabs(XC(vertex1)-XC(vertex2)) > 1e-6) || (fabs(YC(vertex1)-YC(vertex2)) <= 1e-6))
          && ((fabs(XC(vertex1)-XC(vertex2)) <= 1e-6) || (fabs(YC(vertex1)-YC(vertex2)) > 1e-6)) );
  /* now vertex1 and vertex2 are along a horizontal or vertical grid line and
     not on a digonal line */

  meshwidth = fabs( XC(vertex1) - XC(vertex2) );
  if ( meshwidth < 1e-6 )
    meshwidth = fabs( YC(vertex1) - YC(vertex2) );

  return meshwidth;
}

INT FF_PrepareGrid( GRID *grid, DOUBLE *meshwidth, INT init, INT K_comp, INT x_comp, INT b_comp, const BV_DESC_FORMAT *bvdf )
/* traditional FF solver for stripewise decomposition */
/* points must be ordered lexicographic, boundary nodes at the end of the list */
/* must be a square grid */
{
  INT n, nn;
  register VECTOR *v;
  register MATRIX *m, *mnext;
  BLOCKVECTOR *bv_inner;
  BV_DESC bvd_dirichlet;

  /* determine the meshwidth for a regular mesh */
  *meshwidth = FFMeshwidthOfGrid( grid );

  nn = NVEC( grid );

  /* if there exists already a blockvector list in the grid, dispose it first */
  FreeAllBV( grid );

#ifdef __TWODIM__
  n = (INT)(sqrt((double)nn) + 0.00001 );       /* rounding */

  if ( n * n != nn )
  {
    PrintErrorMessage( 'E', "FF_PrepareGrid", "grid is not a square!" );
    REP_ERR_RETURN(1);
  }

  /* remove the boundary nodes from the calculations */
  n = n - 2;

  if ( CreateBVStripe2D( grid, n*n, n ) != GM_OK )
  {
    PrintErrorMessage( 'F', "FF_PrepareGrid", "can not build blockvector structure" );
    REP_ERR_RETURN(1);
  }

#endif

#ifdef __THREEDIM__
  n = (INT)(pow( (double)nn, 1.0/3.0 ) + 0.00001 );       /* rounding */

  if ( n * n * n != nn )
  {
    PrintErrorMessage( 'E', "FF_PrepareGrid", "grid is not a cube!" );
    REP_ERR_RETURN(1);
  }

  /* remove the boundary nodes from the calculations */
  n = n - 2;

  if ( CreateBVStripe3D( grid, n*n*n, n, n ) != GM_OK )
  {
    PrintErrorMessage( 'F', "FF_PrepareGrid", "can not build blockvector structure" );
    REP_ERR_RETURN(1);
  }

#endif

  /* check consistency of the blockvector structure */
  bv_inner = GFIRSTBV( grid );
  ASSERT( !BV_IS_LEAF_BV( bv_inner ) );
  ASSERT( BV_IS_LEAF_BV( BVSUCC(bv_inner) ) );      /* blockvector for the boundary vectors */
  ASSERT( BVSUCC(BVSUCC(bv_inner)) == NULL );       /* the boundary vectors are the last */
  BVD_INIT( &bvd_dirichlet );
  BVD_PUSH_ENTRY( &bvd_dirichlet, BVNUMBER( BVSUCC(bv_inner) ), bvdf );

  if ( init )
  {
    /* put the dirichlet values (stored in u_dirichlet) into the rhs for
       the inner points: f_inner -= K * u_dirichlet
       and correct the stiffness matrix */
    dmatmul_minusBS( bv_inner, &bvd_dirichlet, bvdf, b_comp, K_comp, x_comp );
    dmatsetBS( bv_inner, &bvd_dirichlet, bvdf, K_comp, 0.0 );

    /* remove all 0-connection caused by the triangles */
    for ( v = FIRSTVECTOR( grid ); v != NULL; v = SUCCVC( v ) )
      for ( m = VSTART( v ); m != NULL; m = mnext )
      {
        /*printf( "K[%d,%d] = %g\n", VINDEX(v), VINDEX(MDEST(m)), MVALUE(m,K_comp));*/
        mnext = MNEXT( m );
        if ( (fabs(MVALUE(m,K_comp)) < 1e-15) && (fabs(MVALUE(MADJ(m),K_comp)) < 1e-15) )
          if ( DisposeConnection( grid, MMYCON( m ) ) != 0 )
            PrintErrorMessage( 'E', "FF_PrepareGrid", "error in disposing connection ############\n" );
      }
  }

  FFmuchBigger = 100.0;
  FFsmallTV = 1e-3;
  FFaccuracy = 1e-10;

  mute_level = GetMuteLevel();

#ifdef weglassen

#ifdef __TWODIM__
  UserWrite( "FF2D prepared.\n" );
#endif

#ifdef __THREEDIM__
  UserWrite( "FF3D prepared.\n" );
#endif

#ifdef THETA_EXACT
  UserWrite( "THETA_EXACT is active\n" );
  printf( "THETA_EXACT is active\n" )
#else
  UserWrite( "THETA_EXACT is not active\n" );
  /*printf( "THETA_EXACT is not active\n" );*/
#endif

#ifdef MINV_2D_EXACT
  UserWrite( "MINV_2D_EXACT is active\n" );
  printf( "MINV_2D_EXACT is active\n" );
#else
  UserWrite( "MINV_2D_EXACT is not active\n" );
  /*printf( "MINV_2D_EXACT is not active\n" );*/
#endif

#ifdef THETA_ANA
  UserWrite( "THETA_ANA is active\n" );
  printf( "THETA_ANA is active\n" );
#else
  UserWrite( "THETA_ANA is not active\n" );
  /*printf( "THETA_ANA is not active\n" );*/
#endif

#endif /* weglassen */

  /*printmBS( bv_inner, bv_inner, 0 );*/
  printBVgrid( grid, bvdf );

  return 0;
}


/****************************************************************************/
/*D
   FFConstructTestvector - construct a sine-shaped testvector (2D/3D)

   SYNOPSIS:
   void FFConstructTestvector( const BLOCKVECTOR *bv, INT tv_comp,
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
   circumvent this problem use 'FFConstructTestvector_loc'.

   SEE ALSO:
   FFConstructTestvector_loc, CreateBVStripe2D, CreateBVStripe3D

   RETURN VALUE:
   .n   void
   D*/
/*************************************************************************/

void FFConstructTestvector( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr, DOUBLE wavenr_3D )
{
  register DOUBLE hkpi, pos;
  register VECTOR *v, *end_v;
  register BLOCKVECTOR *bv_glob_end;
  INT length;
#       ifdef __THREEDIM__
  register BLOCKVECTOR *bv_end,*bv_i;
  register DOUBLE plane_pos, tensor;
  DOUBLE plane_hkpi;
  INT plane_length;
#       endif

  bv_glob_end = BVDOWNBVEND(bv);
  bv = BVDOWNBV(bv);
  for ( ; bv != bv_glob_end; bv = BVSUCC(bv) )       /* over all lines resp. planes */
  {
#ifdef __TWODIM__
    ASSERT( BVDOWNTYPE(bv) == BVDOWNTYPEVECTOR );
    length = BVNUMBEROFVECTORS(bv) + 1;
    hkpi = pos = ( PI * wavenr ) / (DOUBLE)length;
    end_v = BVENDVECTOR( bv );
    BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )                /* over all points in the line */
    {
      VVALUE( v, tv_comp ) = sin( pos );
      pos += hkpi;
    }
#else
    ASSERT( !BV_IS_LEAF_BV(bv) );
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
   MeshwidthForFFConstructTestvector_loc - calculate meshwodth an other variables

   SYNOPSIS:
   void MeshwidthForFFConstructTestvector_loc( const VECTOR *v,
             const VECTOR *vn, DOUBLE *meshwidth, DOUBLE *pos );

   PARAMETERS:
   .  v - first vector
   .  vn - second vector
   .  meshwidth - meshwidth between 1. and 2. vector
   .  coord - coordinate in the line of the 1. vector relative to the begining

   DESCRIPTION:
   Determines weather the 2 vectors are horizontal or vertical neighbours
   (the coordinates are considered). The distance along the (horizontal
   or vertical) line is returned as 'meshwidth'.

   Also the 'coord' is returned as the distance of the first vector from the
   beginning of its line (the line is assumed to be in a unit square/cube,
   i.e. the line is [0,1]).

   RESTRICTIONS:
   The grid must be ordered lexicographically by a blockvector structure
   as created by 'CreateBVStripe2D'/'3D'.

   The grid should be regular (i.e. each line/plane contains the same number of
   vectors) and part of the unit square/cube.

   SEE ALSO:
   FFConstructTestvector_loc, CreateBVStripe2D, CreateBVStripe3D

   RETURN VALUE:
   .n   void
   D*/
/*************************************************************************/
static void MeshwidthForFFConstructTestvector_loc( const VECTOR *v, const VECTOR *vn, DOUBLE *meshwidth, DOUBLE *coord )
{
  DOUBLE pos1[DIM], pos2[DIM];

  VectorPosition(v,pos1);
  VectorPosition(vn,pos2);

  *coord = pos1[0];
  /* on the same horizontal line? */
  if( (*meshwidth = fabs( *coord - pos2[0] )) > 1e-6 )
  {
    /* assert that the 2 points are on the same vertical line */
    /*		ASSERT( fabs(pos1[1] - pos2[1]) <= 1e-6 );*/ /* comment out because of cross points */
  }
  else
  {
    *coord = pos1[1];
    if( (*meshwidth = fabs( *coord - pos2[1] )) <= 1e-6 )
    {
      /* if not on the same vertical line then an error */
      /*			ASSERT( FALSE );*/ /* comment out because of cross points */
    }
  }

  /*if( *coord > 0.5 )
   *coord = 1.0 - *coord;*/	/* warum war das drin ???? */
}


void CalculateTv_loc_linesegments( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr )
/* isolates the single linesegments (leaf blockvectors) and calculates for
   each of these segments an ordinary (1D) testvector */
{
  register DOUBLE pos, hkpi;
  register VECTOR *v, *end_v;
  register BLOCKVECTOR *bv_i, *bv_end;
  DOUBLE meshwidth, coord, line_pos;

  if ( BV_IS_EMPTY( bv ) )
    return;

#ifdef QQQTTTT /* auskommentiert fuer Variante, die sich an den Vector-Koordinaten oreintiert */
  if ( BV_IS_LEAF_BV(bv) )
  {             /* go simply over the vectors of the blockvector */
    v = BVFIRSTVECTOR(bv);
    MeshwidthForFFConstructTestvector_loc( v, SUCCVC(v), &meshwidth, &coord );
    hkpi = PI * wavenr * meshwidth;
    pos = PI * coord * wavenr;

    end_v = BVENDVECTOR( bv );
    BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )                /* over all points in the line */
    {
      /*printf("pos %g  sin %g\n", pos, sin(pos) );*/
      VVALUE( v, tv_comp ) = sin( pos );
      pos += hkpi;
    }

    /*printf("leaf tv ; wavenr %g    mesh %g    coord %g\n", wavenr, meshwidth, coord ); printvBS( bv, tv_comp );*/
  }
#endif
  if ( BV_IS_LEAF_BV(bv) )
  {
    DOUBLE_VECTOR position;
    register DOUBLE kpi;
    register const INT orientation = (BVORIENTATION(bv)==BVHORIZONTAL ? 0 : 1);

    /* go simply over the vectors of the blockvector */
    kpi = PI * wavenr;
    end_v = BVENDVECTOR( bv );
    BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )                /* over all points in the line */
    {
      VectorPosition( v, position );
      /*printf("pos %g  sin %g\n", pos, sin(pos) );*/
      VVALUE( v, tv_comp ) = sin( position[orientation] * kpi );
    }

    /*printf("leaf tv ; wavenr = %g,  orientation = %d\n", wavenr, orientation ); printvBS( bv, tv_comp );*/
  }
  else
  {
    bv_i = BVDOWNBV( bv );
    bv_end = BVDOWNBVEND( bv );

    for ( ; bv_i != bv_end; bv_i = BVSUCC( bv_i) )
    {
      ASSERT( BVTVTYPE(bv_i) == BV1DTV );
      CalculateTv_loc_linesegments( bv_i, tv_comp, wavenr );
    }
    /*printf("inner tv \n"); printvBS( bv, tv_comp );*/
  }
}


/****************************************************************************/
/*D
   FFConstructTestvector_loc - construct a sine-shaped testvector (2D/3D)

   SYNOPSIS:
   void FFConstructTestvector_loc( const BLOCKVECTOR *bv, INT tv_comp,
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
   FFConstructTestvector, CreateBVStripe2D, CreateBVStripe3D

   RETURN VALUE:
   .n   void
   D*/
/*************************************************************************/

void FFConstructTestvector_loc( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr, DOUBLE wavenr_3D )
{
  register DOUBLE tensor, pos, hkpi, plane_hkpi, plane_pos;
  register VECTOR *v, *end_v;
  register BLOCKVECTOR *bv_i, *bv_end;
  DOUBLE meshwidth, coord, line_pos;

  if ( BVTVTYPE( bv ) != BV2DTV )
  /* line(s) block */
  {
    /* zick-zack instead of sin */
    /*DOUBLE incr = (2.0 * wavenr) / (BVNUMBEROFVECTORS(bv) + 1.0 );
            DOUBLE pos = incr;

            end_v = BVENDVECTOR( bv );
            BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )
            {
                    VVALUE( v, tv_comp ) = pos;
                    if ( fabs(pos) > 0.999 ) incr *= -1.0;
                    pos += incr;
            }
            return;
     */

    CalculateTv_loc_linesegments( bv, tv_comp, wavenr );

#ifdef QQQTTT
    if ( BV_IS_LEAF_BV(bv) )
    {                   /* go simply over the vectors of the blockvector */
      v = BVFIRSTVECTOR(bv);
      MeshwidthForFFConstructTestvector_loc( v, SUCCVC(v), &meshwidth, &coord );
      hkpi = PI * wavenr * meshwidth;
      pos = PI * coord * wavenr;

      end_v = BVENDVECTOR( bv );
      BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )                      /* over all points in the line */
      {
        /*printf("pos %g  sin %g\n", pos, sin(pos) );*/
        VVALUE( v, tv_comp ) = sin( pos );
        pos += hkpi;
      }
    }
    else
    {}
#endif

  }
  else
  /* plane block; tensor product testvector */
  {
    ASSERT( !BV_IS_LEAF_BV(bv) );
    bv_i = BVDOWNBV(bv);
    bv_end = BVDOWNBVEND(bv);
    v = BVFIRSTVECTOR(bv_i);

    /* planewise */
    MeshwidthForFFConstructTestvector_loc( v, BVFIRSTVECTOR(BVSUCC(bv_i)), &meshwidth, &coord );
    plane_hkpi = PI * wavenr * meshwidth;
    plane_pos = PI * coord * wavenr;
    tensor = sin( plane_pos );

    /* linewise */
    MeshwidthForFFConstructTestvector_loc( v, SUCCVC(v), &meshwidth, &coord );
    hkpi = PI * wavenr * meshwidth;
    line_pos = PI * coord * wavenr;

    for ( ; bv_i != bv_end; bv_i = BVSUCC(bv_i) )               /* over all lines */
    {
      pos = line_pos;
      end_v = BVENDVECTOR( bv_i );
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
   FFMultWithM - for a frequency filtered matrix M calculate y := M * x

   SYNOPSIS:
   INT FFMultWithM( const BLOCKVECTOR *bv,
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
   'M' is an frequency filtered decomposed matrix as produced by 'TFFDecomp'
   or 'FFDecomp', i.e. M = ( L + T ) * T^-1 * ( T + U ).
   This function calculates the product 'y := M * x'. Do not mix up this
   function with 'FFMultWithMInv'!

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
   TFFDecomp, FFDecomp, FFMultWithMInv

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   D*/
/****************************************************************************/

INT FFMultWithM( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT y_comp, INT x_comp )
{
  register BLOCKVECTOR *bv_i, *bv_ip1, *bv_stop;
  register BV_DESC *bvd_i, *bvd_ip1, *bvd_temp;
  BV_DESC bvd1, bvd2;
  INT aux_comp, L_comp, T_comp;

  aux_comp = GET_AUX_VEC;
  L_comp = STIFFMAT_ON_LEVEL(bv);
  T_comp = DECOMPMAT_ON_LEVEL(bv);

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
    dmatmul_addBS( bv_i, bvd_ip1, bvdf, aux_comp, L_comp, x_comp );

    /* aux_i := (T_i)^-1 * aux_i */
#ifdef ModelP
    FFMultWithMInv( bv_i, bvd_i, bvdf, aux_comp, aux_comp, NULL, NULL );
#else
    FFMultWithMInv( bv_i, bvd_i, bvdf, aux_comp, aux_comp );
#endif

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
    dmatmul_addBS( bv_ip1, bvd_ip1, bvdf, y_comp, T_comp, aux_comp );

    /* y_i+1 += L_(i+1,i) * aux_i */
    dmatmul_addBS( bv_ip1, bvd_i, bvdf, y_comp, L_comp, aux_comp );

    /* prepare BVDs for next loop */
    SWAP( bvd_i, bvd_ip1, bvd_temp );
    BVD_DEC_LAST_ENTRY( bvd_i, 2, bvdf );
  }
  /* y_0 := 0 */
  dsetBS( bv_ip1, y_comp, 0.0 );

  /* y_0 += T_0 * aux_0 */
  dmatmul_addBS( bv_ip1, bvd_ip1, bvdf, y_comp, T_comp, aux_comp );

  FREE_AUX_VEC(aux_comp);

  return NUM_OK;
}


/****************************************************************************/
/*D
   FFMultWithMInv - for a frequency filtered matrix M calculate v := M^-1 * b

   SYNOPSIS:
   INT FFMultWithMInv( const BLOCKVECTOR *bv,
                                                const BV_DESC *bvd,
                                                const BV_DESC_FORMAT *bvdf,
                                                INT v_comp,
                                                INT b_comp );

   PARAMETERS:
   .  bv - blockvector covering M
   .  bvd - description of the blockvector
   .  bvdf - format to interpret the 'bvd'
   .  v_comp - position of the result in the VECTOR-data
   .  b_comp - position of the right hand side vector in the VECTOR-data

   DESCRIPTION:
   'M' is an frequency filtered decomposed matrix as produced by 'TFFDecomp'
   or 'FFDecomp', i.e. M = ( L + T ) * T^-1 * ( T + U ).
   This function calculates 'v := M^-1 * b', i.e. solves 'M * v = b' for v.
   Do not mix up this function with 'FFMultWithM'!

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
   TFFDecomp, FFDecomp, FFMultWithM

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    error code from 'solveLUMatBS'
   D*/
/****************************************************************************/

INT FFMultWithMInv(
  const BLOCKVECTOR *bv,
  const BV_DESC *bvd,
  const BV_DESC_FORMAT *bvdf,
  INT v_comp,
  INT b_comp
#ifdef ModelP
  ,const VECDATA_DESC *v,                                       /* braucht man nicht mehr! entfernen !!!!*/
  GRID *grid
#endif
  )
{
  register BLOCKVECTOR *bv_i, *bv_ip1, *bv_stop;
  register BV_DESC *bvd_i=NULL, *bvd_ip1, *bvd_temp;
  BV_DESC bvd1, bvd2;
  BLOCKVECTOR *bv_first;
  INT aux_comp, auxsub_comp, L_comp;
#ifdef MINV_2D_EXACT
  INT auxA_comp, auxB_comp;
#endif

  ASSERT( !BV_IS_EMPTY(bv) );

  if ( BV_IS_LEAF_BV(bv) )
  {
    solveLUMatBS( bv, bvd, bvdf, v_comp, DECOMPMAT_ON_LEVEL(bv), b_comp );
    return NUM_OK;
  }

  if ( BV_IS_DIAG_BV(bv) )
  {
    bvd1 = *bvd;
    bv_stop = BVDOWNBVEND(bv);
    for ( bv_i = BVDOWNBV( bv ); bv_i != bv_stop; bv_i = BVSUCC( bv_i ) )
      if( !BV_IS_EMPTY(bv_i) )
      {
        BVD_PUSH_ENTRY( &bvd1, BVNUMBER(bv_i), bvdf );
#ifdef ModelP
        FFMultWithMInv( bv_i, &bvd1, bvdf, v_comp, b_comp, NULL, NULL );
#else
        FFMultWithMInv( bv_i, &bvd1, bvdf, v_comp, b_comp );
#endif
        BVD_DISCARD_LAST_ENTRY(&bvd1);
      }

    return(NUM_OK);
  }

#ifdef MINV_2D_EXACT
  auxA_comp = GET_AUX_VEC;
  ASSERT( auxA_comp>-1 );
  auxB_comp = GET_AUX_VEC;
  ASSERT( auxB_comp>-1 );
#endif

  aux_comp = GET_AUX_VEC;
  ASSERT( aux_comp>-1 );
  L_comp = STIFFMAT_ON_LEVEL( bv );

  ASSERT( v_comp != aux_comp );
  ASSERT( b_comp != aux_comp );

  /* To minimize the incrementation of BVDs there are used two (one for
     index i and the other for index i+1) which are swapped in the loop;
     thus only one incrementation by 2 is necessary */
  bvd1 = bvd2 = *bvd;
  bvd_i = &bvd1;
  bvd_ip1 = &bvd2;

  /* solve lower triangular matrix; except the last block in this loop */
  /* aux := ( L + T )^-1 * b */

  /* set up stop-block (last non-empty block) */
  bv_stop = BVDOWNBVLAST(bv);
  while( BV_IS_EMPTY( bv_stop ) && (bv_stop != BVDOWNBV(bv)) )          /* search last nonempty bv */
    bv_stop = BVPRED( bv_stop );

  /* set up first block */
  bv_i = BVDOWNBV(bv);
  while( BV_IS_EMPTY( bv_i ) && (bv_i != BVDOWNBVEND(bv)) )             /* search first nonempty bv */
    bv_i = BVSUCC( bv_i );
  ASSERT( bv_i != BVDOWNBVEND(bv) );
  BVD_PUSH_ENTRY( bvd_i, BVNUMBER(bv_i), bvdf );
  bv_first = bv_i;

  /* set up second block */
  bv_ip1 = BVSUCC( bv_i );
  while( (bv_ip1 != BVDOWNBVEND(bv)) && BV_IS_EMPTY( bv_ip1 ) )         /* search first nonempty bv */
    bv_ip1 = BVSUCC( bv_ip1 );
  if ( bv_ip1 != BVDOWNBVEND(bv) )
    BVD_PUSH_ENTRY( bvd_ip1, BVNUMBER(bv_ip1), bvdf );
  /* else: bv_ip1 and bvd_ip1 are never used;
     thus the content of them have not to be updated */

  for ( ; bv_i != bv_stop; )
  {
    /* aux_i := (T_i)^-1 * b_i */
#ifdef MINV_2D_EXACT
    if ( BV_IS_LEAF_BV(bv_i) )
      solveLUMatBS( bv_i, bvd_i, bvdf, aux_comp, DECOMPMAT_ON_LEVEL( bv_i ), b_comp );
    else
      gs_solveBS ( bv_i, bvd_i, bvdf, 5e-14, 1000, L_comp, aux_comp, b_comp, auxA_comp, TRUE, TRUE );
#else
                #ifdef ModelP
    FFMultWithMInv( bv_i, bvd_i, bvdf, aux_comp, b_comp, NULL, NULL );

    /* in the case of ModelP: make the new solution on the interface consistent */
    if( BVNUMBER(bv_i) == -100 )             /* lines blockvector */
#ifdef FFCOMM
      FFVectorConsistent( bv_i, aux_comp );
#else
      if( l_vector_consistentBS( grid, bvd_i, bvdf, aux_comp )!=NUM_OK ) REP_ERR_RETURN (1);
#endif
                #else
    FFMultWithMInv( bv_i, bvd_i, bvdf, aux_comp, b_comp );
                #endif
#endif
    /* b_i+1 -= L_(i+1,i) * aux_i */
    dmatmul_minusBS( bv_ip1, bvd_i, bvdf, b_comp, L_comp, aux_comp );

    /* prepare BVDs for next loop */
    bv_i = bv_ip1;
    SWAP( bvd_i, bvd_ip1, bvd_temp );
    bv_ip1 = BVSUCC( bv_ip1 );
    while( (bv_ip1 != BVDOWNBVEND(bv)) && BV_IS_EMPTY( bv_ip1 ) )               /* search first nonempty bv */
      bv_ip1 = BVSUCC( bv_ip1 );
    if( bv_ip1 != BVDOWNBVEND(bv) )
    {
      BVD_DISCARD_LAST_ENTRY( bvd_ip1 );
      BVD_PUSH_ENTRY( bvd_ip1, BVNUMBER(bv_ip1), bvdf );
    }
    /* else: bv_ip1 and bvd_ip1 are never used;
       thus the content of them have not to be updated */
  }
  /* special treatment: v_last = (T_last)^-1 * b_last */
#ifdef MINV_2D_EXACT
  /*if ( BV_IS_LEAF_BV(bv_i) )
          solveLUMatBS( bv_i, bvd_i, bvdf, v_comp, DECOMPMAT_ON_LEVEL( bv_i ), b_comp );
     else*/
  {
    dcopyBS( bv_i, auxB_comp, b_comp );                 /* necessary if b_comp == v_comp */
    gs_solveBS ( bv_i, bvd_i, bvdf, 5e-14, 1000, L_comp, v_comp, auxB_comp, auxA_comp, TRUE, TRUE );
  }
#else
        #ifdef ModelP
  FFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, b_comp, NULL, NULL );

  /* lines or cross blockvector */
  if( (BVNUMBER(bv_i) == -100) || (BVNUMBER(bv_i) == -101) )
  {
#ifdef FFCOMM
    if ( BVNUMBER(bv_i) == -100 )
      FFVectorConsistent( bv_i, v_comp );
    else
    if( l_vector_consistentBS( grid, bvd_i, bvdf, v_comp )!=NUM_OK ) REP_ERR_RETURN (1);
#else
    if( l_vector_consistentBS( grid, bvd_i, bvdf, v_comp )!=NUM_OK ) REP_ERR_RETURN (1);
#endif
  }
#else
  FFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, b_comp );
        #endif
#endif

  /* solve upper triangular matrix; the last block is already calculated */
  /* v := (T^-1*U + I )^-1 * aux */
  ASSERT( bv_i == bv_stop );

  /* set up stop block */
  bv_stop = BVPRED( bv_first );

  /* set up last block (bvd_ip1, bv_ip1 is not used)*/
  SWAP( bvd_i, bvd_ip1, bvd_temp );

  /* set up last but one block (bv(d)_i )*/
  bv_i = BVPRED( bv_i );
  while( (bv_i != bv_stop) && BV_IS_EMPTY( bv_i ) )             /* search first nonempty bv */
    bv_i = BVPRED( bv_i );
  if ( bv_i != bv_stop )
  {
    BVD_DISCARD_LAST_ENTRY( bvd_i );
    BVD_PUSH_ENTRY( bvd_i, BVNUMBER(bv_i), bvdf );
  }

  for ( ; bv_i != bv_stop; )
  {
    /* v_i := L_(i,i+1) * v_i+1 */
    dsetBS( bv_i, v_comp, 0.0 );
    dmatmul_addBS( bv_i, bvd_ip1, bvdf, v_comp, L_comp, v_comp );

    /* v_i := (T_i)^-1 * v_i */
#ifdef MINV_2D_EXACT
    /*if ( BV_IS_LEAF_BV(bv_i) )
            solveLUMatBS( bv_i, bvd_i, bvdf, v_comp, DECOMPMAT_ON_LEVEL( bv_i ), b_comp );
       else*/
    {
      dcopyBS( bv_i, auxB_comp, v_comp );
      gs_solveBS ( bv_i, bvd_i, bvdf, 5e-14, 1000, L_comp, v_comp, auxB_comp, auxA_comp, TRUE, TRUE );
    }
#else
                #ifdef ModelP
    FFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, v_comp, NULL, NULL );

    /* in the case of ModelP: make the new solution on the interface consistent */
    if( BVNUMBER(bv_i) == -100 )             /* lines blockvector */
#ifdef FFCOMM
      FFVectorConsistent( bv_i, v_comp );
#else
      if( l_vector_consistentBS( grid, bvd_i, bvdf, v_comp )!=NUM_OK ) REP_ERR_RETURN (1);
#endif
                #else
    FFMultWithMInv( bv_i, bvd_i, bvdf, v_comp, v_comp );
                #endif
#endif

    /* v_i := aux_i - v_i */
    dminusaddBS( bv_i, v_comp, aux_comp );

    /* prepare BVDs for next loop */
    SWAP( bvd_i, bvd_ip1, bvd_temp );
    bv_i = BVPRED( bv_i );
    while( (bv_i != bv_stop) && BV_IS_EMPTY( bv_i ) )                   /* search first nonempty bv */
      bv_i = BVPRED( bv_i );
    if( bv_i != bv_stop )
    {
      BVD_DISCARD_LAST_ENTRY( bvd_i );
      BVD_PUSH_ENTRY( bvd_i, BVNUMBER(bv_i), bvdf );
    }
    /* else: bv_i and bvd_i are never used;
       thus the content of them have not to be updated */
  }

  FREE_AUX_VEC( aux_comp );

#ifdef MINV_2D_EXACT
  FREE_AUX_VEC( auxB_comp );
  FREE_AUX_VEC( auxA_comp );
#endif

  return NUM_OK;
}

#ifdef QQQQQQQQQQQQ
INT FFMultWithMInvDD( GRID *grid,
                      const BV_DESC_FORMAT *bvdf,
                      const VECDATA_DESC *v,
                      const VECDATA_DESC *b )
{
  register BLOCKVECTOR *bv_subdom, *bv_lines;
  BV_DESC bvd_subdom, bvd_lines;
  INT v_comp, b_comp, aux_comp, L_comp;
  bv_subdom = GFIRSTBV(grid);
  bv_lines = BVSUCC(bv_subdom);

  ASSERT( !BV_IS_LEAF_BV(bv_subdom) );
  ASSERT( !BV_IS_LEAF_BV(bv_lines) );
  ASSERT( VD_IS_SCALAR(v) );
  ASSERT( VD_IS_SCALAR(b) );

  v_comp = VD_SCALCMP(v);
  b_comp = VD_SCALCMP(b);
  L_comp = STIFFMAT_ON_LEVEL( bv_subdom );
  ASSERT( L_comp>-1 );

  aux_comp = GET_AUX_VEC;
  ASSERT( aux_comp>-1 );

  ASSERT( v_comp != aux_comp );
  ASSERT( b_comp != aux_comp );

  BVD_INIT( &bvd_subdom );
  BVD_PUSH_ENTRY( &bvd_subdom, BVNUMBER(bv_subdom), bvdf );
  BVD_INIT( &bvd_lines );
  BVD_PUSH_ENTRY( &bvd_lines, BVNUMBER(bv_lines), bvdf );

  /* solve lower triangular matrix; except the last block */
  /* aux := ( L + T )^-1 * b */

  /* aux_subdom := (T_subdom)^-1 * b_subdomi */
  FFMultWithMInv( bv_subdom, &bvd_subdom, bvdf, aux_comp, b_comp );

  /* b_lines -= L_(lines,subdom) * aux_subdom */
  dmatmul_minusBS( bv_lines, &bvd_subdom, bvdf, b_comp, L_comp, aux_comp );

  /* special treatment: v_last = (T_last)^-1 * b_last */
  /* aux_lines := (T_lines)^-1 * b_lines */
  FFMultWithMInv( bv_lines, &bvd_lines, bvdf, v_comp, b_comp );

  /* in the case of ModelP: make the new solution on the interface consistent */
  if( l_vector_consistentBS( grid, v )!=NUM_OK ) REP_ERR_RETURN (1);

  /* cross points still missing */


  /* solve upper triangular matrix; the last block is already calculated */
  /* v := (T^-1*U + I )^-1 * aux */

  /* v_subdom := L_(subdom,lines) * v_lines */
  dsetBS( bv_subdom, v_comp, 0.0 );
  dmatmul_addBS( bv_subdom, bvd_lines, bvdf, v_comp, L_comp, v_comp );

  /* v_subdom := (T_subdom)^-1 * v_subdom */
  FFMultWithMInv( bv_subdom, bvd_subdom, bvdf, v_comp, v_comp );

  /* v_subdom := aux_subdom - v_subdom */
  dminusaddBS( bv_subdom, v_comp, aux_comp );

  FREE_AUX_VEC( aux_comp );

  return NUM_OK;
}
#endif /* QQQQQQQQQQQ */

#endif /* __BLOCK_VECTOR_DESC__ */
