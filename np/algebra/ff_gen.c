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

#include "compiler.h"
#include "gm.h"          /* for data structure               */
#include "devices.h"     /* for UserWrite, PrintErrorMessage */
#include "commands.h"    /* for GetCurrentMultigrid              */
#include "cmdint.h"      /* for CreateCommand                */

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

/* mute level for FFs */
INT mute_level  = 0;

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
      return GM_OUT_OF_MEM;

    }
    if ( (mem = GetMem( MGHEAP(MYMG(grid)), BVNUMBEROFVECTORS(bv)*sizeof(DOUBLE), FROM_BOTTOM )) == NULL )
    {
      PrintErrorMessage( 'E', "storeVectorBS", "Not enough memory to store the vector" );
      return GM_OUT_OF_MEM;
    }
    BVUSERDATA( bv ) = mem;
  }

  end_v = BVENDVECTOR( bv );
  BLOCK_L_VLOOP( v, BVFIRSTVECTOR(bv), end_v )
  *mem++ = VVALUE( v, x_comp );

  /* the number of copied elements must be equal the number of vectors
          stored in the blockvector for them the memory is allocated */
  assert( ( mem - BACKUP_MEM(bv) ) == BVNUMBEROFVECTORS( bv ) );

  return NUM_OK;
}


INT     restoreVectorBS( BLOCKVECTOR *bv, INT x_comp )
/* restore the given vector from an allocated buffer in the blockvector */
{
  register DOUBLE *mem;
  register VECTOR *v, *end_v;

  mem = BACKUP_MEM( bv );

  assert( mem != NULL );

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

  assert(FALSE); /* neues boundary konzept einbauen */
  return 1;

  register VSEGMENT *vseg;
  register BndCondProcPtr bnd_cond;
  register BNDSEGDESC *bnd_desc;

  DOUBLE value;
  INT type;

  /* only node vectors are considered */
  assert( VTYPE(v) == NODEVECTOR );

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

  printf("comp (%d)\n",m_nr);
  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
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

  printf("comp (%d)\n",m_nr);
  for (v=FIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
  {
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

  printf("comp (%d)\n",m_nr);

  if ( BVNUMBEROFVECTORS( bv_row ) == 0 || BVNUMBEROFVECTORS( bv_col ) == 0 )
  {
    printf( "empty\n" );
    return;
  }

  for (v=BVFIRSTVECTOR(bv_row); v!= BVENDVECTOR(bv_row); v=SUCCVC(v))
  {
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

  printf("comp (%d)\n",m_nr);

  if ( BVNUMBEROFVECTORS( bv_row ) == 0 )
  {
    printf( "empty\n" );
    return;
  }

  for (v=BVFIRSTVECTOR(bv_row); v!= BVENDVECTOR(bv_row); v=SUCCVC(v))
  {
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
    printf( "%s Nr. %d ", indent, BVNUMBER(bv) );
    if ( BVNUMBEROFVECTORS(bv) != 0 )
    {
      printf( "number of vectors %d ", BVNUMBEROFVECTORS(bv) );
      printf( "first vector %d ", VINDEX(BVFIRSTVECTOR(bv)) );
      printf( "last vector %d\n", VINDEX(BVLASTVECTOR(bv)) );

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
    assert( link != NULL );
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
    PrintErrorMessage( 'E', "FF_PrepareSolver", "grid is not a square!" );
    return 1;
  }

  /* remove the boundary nodes from the calculations */
  n = n - 2;

  if ( CreateBVStripe2D( grid, n*n, n ) != GM_OK )
  {
    PrintErrorMessage( 'F', "FF_PrepareSolver", "can not build blockvector structure" );
    return 1;
  }

#endif

#ifdef __THREEDIM__
  n = (INT)(pow( (double)nn, 1.0/3.0 ) + 0.00001 );       /* rounding */

  if ( n * n * n != nn )
  {
    PrintErrorMessage( 'E', "FF_PrepareSolver", "grid is not a cube!" );
    return 1;
  }

  /* remove the boundary nodes from the calculations */
  n = n - 2;

  if ( CreateBVStripe3D( grid, n*n*n, n, n ) != GM_OK )
  {
    PrintErrorMessage( 'F', "FF_PrepareSolver", "can not build blockvector structure" );
    return 1;
  }

#endif

  /* check consistency of the blockvector structure */
  bv_inner = GFIRSTBV( grid );
  assert( !BV_IS_LEAF_BV( bv_inner ) );
  assert( BV_IS_LEAF_BV( BVSUCC(bv_inner) ) );      /* blockvector for the boundary vectors */
  assert( BVSUCC(BVSUCC(bv_inner)) == NULL );       /* the boundary vectors are the last */
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
            PrintErrorMessage( 'E', "FF_PrepareSolver", "error in disposing connection ############\n" );
      }
  }

  TFFmuchBigger = 100.0;
  TFFsmallTV = 1e-3;
  TFFaccuracy = 1e-10;

  mute_level = GetMuteLevel();

#ifdef weglassen

#ifdef __TWODIM__
  UserWrite( "TFF2D prepared.\n" );
#endif

#ifdef __THREEDIM__
  UserWrite( "TFF3D prepared.\n" );
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

  return 0;
}
#endif /* __BLOCK_VECTOR_DESC__ */
