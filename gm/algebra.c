// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  algebra.c                                                                                                     */
/*																			*/
/* Purpose:   management for algebraic structures							*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 294										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de                                               */
/*																			*/
/*			  blockvector data structure:									*/
/*			  Christian Wrobel                                                                              */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:    1.12.93 begin, ug 3d                                                                             */
/*			  26.10.94 begin combination 2D/3D version						*/
/*			  28.09.95 blockvector implemented (Christian Wrobel)			*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/*		defines to exclude functions										*/
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "compiler.h"
#include "heaps.h"
#include "fifo.h"
#include "ugenv.h"
#include "debug.h"
#include "general.h"

#include "devices.h"

#include "dlmgr.h"
#include "gm.h"
#include "refine.h"
#include "ugm.h"
#include "evm.h"
#include "misc.h"
#include "algebra.h"
#include "dlmgr.h"

#ifdef ModelP
#include "parallel.h"
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

/* for LexAlgDep */
#ifdef __SGI__
/* qsort on SGI cannot handle nontransitive compare routines,               */
/* resulting by rounding errors                                             */
#define ORDERRES                0       /* resolution for LexAlgDep					*/
#else
#define ORDERRES                1e-3    /* resolution for LexAlgDep					*/
#endif

/* for GetDomainPart indicating an element is meant rather than an element side */
#define NOSIDE          -1

/* for LexOrderVectorsInGrid */
#define MATTABLESIZE    32

/* constants for the direction of domain halfening */
#define BV_VERTICAL     0
#define BV_HORIZONTAL   1

/* for ordering of matrices */
#define MATHPOS                 +1
#define MATHNEG                 -1

/* temp vector flag for LineOrderVectors */
static INT ce_VCSTRONG;
#define VCSTRONG_LEN                                    1
#define VCSTRONG(p)                                             CW_READ(p,ce_VCSTRONG)
#define SETVCSTRONG(p,n)                                CW_WRITE(p,ce_VCSTRONG,n)

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

INT MatrixType[MAXVECTORS][MAXVECTORS];
const char *ObjTypeName[MAXVOBJECTS];

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT theAlgDepDirID;                      /* env type for Alg Dep dir                     */
static INT theAlgDepVarID;                      /* env type for ALG_DEP vars                    */

static INT theFindCutDirID;                     /* env type for FindCut dir                     */
static INT theFindCutVarID;                     /* env type for FIND_CUT vars                   */

static FindCutProcPtr FindCutSet;       /* pointer to find cut procedure		*/

static MULTIGRID *GBNV_mg;                      /* multigrid							*/
static INT GBNV_n;                                      /* list items							*/
static INT GBNV_curr;                           /* curr pos								*/
static VECTOR **GBNV_list=NULL;         /* list pointer							*/
static INT GBNV_MarkKey;                        /* key for Mark/Release					*/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/* 2 often used blockvector description formats */
const BV_DESC_FORMAT DH_bvdf =
{ 2, 16,
  { (BVD_ENTRY_TYPE)0x03,        (BVD_ENTRY_TYPE)0x0f,
            (BVD_ENTRY_TYPE)0x03f,       (BVD_ENTRY_TYPE)0x0ff,
            (BVD_ENTRY_TYPE)0x03ff,      (BVD_ENTRY_TYPE)0x0fff,
            (BVD_ENTRY_TYPE)0x03fff,     (BVD_ENTRY_TYPE)0x0ffff,
            (BVD_ENTRY_TYPE)0x03ffff,    (BVD_ENTRY_TYPE)0x0fffff,
            (BVD_ENTRY_TYPE)0x03fffff,   (BVD_ENTRY_TYPE)0x0ffffff,
            (BVD_ENTRY_TYPE)0x03ffffff,  (BVD_ENTRY_TYPE)0x0fffffff,
            (BVD_ENTRY_TYPE)0x3fffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff,  (BVD_ENTRY_TYPE)0xffffffff},
  { (BVD_ENTRY_TYPE)0xfffffffc, (BVD_ENTRY_TYPE)0xfffffff3,
            (BVD_ENTRY_TYPE)0xffffffcf, (BVD_ENTRY_TYPE)0xffffff3f,
            (BVD_ENTRY_TYPE)0xfffffcff, (BVD_ENTRY_TYPE)0xfffff3ff,
            (BVD_ENTRY_TYPE)0xffffcfff, (BVD_ENTRY_TYPE)0xffff3fff,
            (BVD_ENTRY_TYPE)0xfffcffff, (BVD_ENTRY_TYPE)0xfff3ffff,
            (BVD_ENTRY_TYPE)0xffcfffff, (BVD_ENTRY_TYPE)0xff3fffff,
            (BVD_ENTRY_TYPE)0xfcffffff, (BVD_ENTRY_TYPE)0xf3ffffff,
            (BVD_ENTRY_TYPE)0xcfffffff, (BVD_ENTRY_TYPE)0x3fffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
            (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff}};

const BV_DESC_FORMAT one_level_bvdf =
{ 32, 1,
  {     (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff},
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

const BV_DESC_FORMAT two_level_bvdf =
{ 16, 2,
  {     (BVD_ENTRY_TYPE)0x0000ffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff},
  { (BVD_ENTRY_TYPE)0xffff0000, (BVD_ENTRY_TYPE)0x0000ffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

const BV_DESC_FORMAT three_level_bvdf =
{ 10, 3,
  {     (BVD_ENTRY_TYPE)0x000003ff, (BVD_ENTRY_TYPE)0x000fffff,
                (BVD_ENTRY_TYPE)0x3fffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff,
                (BVD_ENTRY_TYPE)0xffffffff, (BVD_ENTRY_TYPE)0xffffffff},
  { (BVD_ENTRY_TYPE)0xfffffc00, (BVD_ENTRY_TYPE)0xfff003ff,
            (BVD_ENTRY_TYPE)0xc00fffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* for LexOrderVectorsInGrid */
static const INT *Order,*Sign;
static INT SkipV,SignRad;
static DOUBLE InvMeshSize;

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static INT ReallocateVector (GRID *g, VECTOR *v, INT part, INT side);

static INT CreateBVPlane( BLOCKVECTOR **bv_plane, const BV_DESC *bvd_plane, const BV_DESC_FORMAT *bvdf, VECTOR **v, INT stripes, INT vectors_per_stripe, GRID *grid );
static INT BlockHalfening( GRID *grid, BLOCKVECTOR *bv, INT left, INT bottom, INT width, INT height, INT side, INT orientation, INT leaf_size );
static INT DisposeIMatrices (GRID *theGrid, MATRIX *theMatrix);

/****************************************************************************/
/*D
   CheckMatrixList - Check matrix list

   SYNOPSIS:
   static void CheckMatrixList (VECTOR *theVector);

   PARAMETERS:
   .  theVector - pointer to a vector list

   DESCRIPTION:
   This function checks matrix list and prints an error message
   if the list is inconsistent.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

static void CheckMatrixList (VECTOR *theVector)
{
  MATRIX *theMatrix;

  for (theMatrix=VSTART(theVector); theMatrix!= NULL; theMatrix=MNEXT(theMatrix))
  {
    if (MDIAG(theMatrix))
    {
      if (theVector != MDEST(theMatrix))
        UserWriteF("matrix ??? should be diagonal matrix but dest. ptr does not point back\n");
      if (VSTART(theVector) != theMatrix)
        UserWriteF("matrix ??? is diagonal matrix but is not at first place in the list\n");
    }
    else
    if (MDEST(MADJ(theMatrix)) != theVector)
      UserWriteF("adj of matrix ??? does not point back\n");
  }
}

/****************************************************************************/
/*D
   InitBVDF - Initialize a blockvector description format

   SYNOPSIS:
   INT InitBVDF( BV_DESC_FORMAT *bvdf, BLOCKNUMBER max_blocks )

   PARAMETERS:
   .  bvdf - blockvector description format to be initialized
   .  max_blocks - maximum number of blocks provided for each blockvector level

   DESCRIPTION:
   All routines working with a 'blockvector description' need a explanation
   of the internal encoding used in the description. This routine generates
   this necessary explanation. Since the encoding depends on the number
   of blockvectors you want to use on each blockvector level, you have to
   specify this quantity. The background is, that a sequence of blockvector
   numbers is encoded in a bitstring and the number of bits used for each
   level has to be fixed before the first use. In order to get fast access
   to the stored information, a couple of bitmasks are necessary and are
   precalculated and stored in the 'blockvector description format' by this
   function. The memory to store the information is provided by the caller.

   HINTS:
   Since the hole amount of bits available for encoding a sequence of
   blockvector numbers is constant, you may use only few blocks on each
   level and for that have many levels available; or if you have many
   blocks on each level you can work only with few levels at all. Thus you
   have to find a trade off of your requirements.

   Since the relevant information of the 'max_blocks' parameter is only
   the number of bits, necessary to encode this number, all numbers up to
   the next power of 2 will have the same result.

   PREDEFINED FORMATS:
   There are predefined 'blockvector description formats', which are
   often used. This formats you can use without initialization.
   .  DH_bvdf  - useful for domain halfening methods (up to 4 blocks per level)
   .  one_level_bvdf - max. number of blocks in a level but only 1 level
   .  two_level_bvdf - bvdf for 2 blocklevels
   .  three_level_bvdf - bvdf for 3 blocklevels

   EXAMPLE:
   To create a format for managing a octree-like structure (that is
   subdividing a 3D cube into 8 subcubes) you can use the following
   created format:
   .vb
   BV_DESC_FORMAT my_bvdf;

   InitBVDF( &my_bvdf, 8 );
   .ve

   RETURN VALUE:
   INT
   .n GM_OK if ok
   .n GM_OUT_OF_RANGE if 'max_blocks' < 2 or if 'max_blocks' greater than possible

   SEE ALSO:
   BLOCKVECTOR, BV_DESC, BV_DESC_FORMAT

   D*/
/****************************************************************************/

INT InitBVDF( BV_DESC_FORMAT *bvdf, BLOCKNUMBER max_blocks )
{
  INT bits;
  INT level;
  BVD_ENTRY_TYPE mask;

  /* calculate number of bits to represent the numbers 0..max_blocks-1,
     i.e. the next integer greater or equal to the logarithmus dualis
     of max_blocks-1 */
  if ( max_blocks < (BLOCKNUMBER)2 )
    return GM_OUT_OF_RANGE;
  max_blocks--;
  bits = 0;
  while( max_blocks > (BLOCKNUMBER)0 )
  {
    max_blocks >>= 1;
    bits++;
  }


  bvdf->bits = bits;
  bvdf->max_level = BVD_MAX_ENTRIES / bits;
  if ( bvdf->max_level == 0 )
    return GM_OUT_OF_RANGE;

  /* calculate the necessary bitmasks */
  mask = ((BVD_ENTRY_TYPE)1 << bits) - (BVD_ENTRY_TYPE)1;
  bvdf->level_mask[0] = mask;
  bvdf->neg_digit_mask[0] = ~mask;
  for ( level = 1; level < BVD_MAX_ENTRIES; level ++ )
  {
    mask <<= bits;
    bvdf->level_mask[level] = bvdf->level_mask[level-1] | mask;
    bvdf->neg_digit_mask[level] = ~mask;
  }

  return GM_OK;
}

/****************************************************************************/
/*D
   PushEntry - Adds a further blockvector number to the existing sequence

   SYNOPSIS:
   INT PushEntry( BV_DESC *bvd, BLOCKNUMBER bnr, const BV_DESC_FORMAT *bvdf )

   PARAMETERS:
   .  bvd - blockvector description to which a further blocknumber is added
   .  bnr - blocknumber to be added
   .  bvdf - blockvector description format to work on the 'bvd'

   DESCRIPTION:
   In a blockvector description the position of a blockvector in a blockvector
   hierarchy is stored as a sequence of blocknumbers. With this function you
   can add a further number to the existing sequence, i.e. to describe
   a subblock of the smallest already existing block in the 'bvd'.

   This function is embedded in the macro 'BVD_PUSH_ENTRY' and provides
   together with the macro 'BVD_DISCARD_LAST_ENTRY' a stack-mechanism
   for storing a sequence of blockvector numbers.

   RETURN VALUE:
   INT
   .n GM_OK if ok
   .n GM_OUT_OF_RANGE if the maximum number of blockvector levels is exeeded

   SEE ALSO:
   BLOCKVECTOR, BV_DESC

   D*/
/****************************************************************************/

INT PushEntry( BV_DESC *bvd, BLOCKNUMBER bnr, const BV_DESC_FORMAT *bvdf )
{
  /* exist, if there is no space for a further number */
  if ( bvd->current >= bvdf->max_level )
    return GM_OUT_OF_RANGE;

  /* store the number */
  BVD_SET_ENTRY( bvd, bvd->current, bnr, bvdf );
  bvd->current++;

  return GM_OK;
}

/****************************************************************************/
/*D
   FindBV - Find the described blockvector in the grid

   SYNOPSIS:
   BLOCKVECTOR *FindBV( const GRID *grid, BV_DESC *bvd, const BV_DESC_FORMAT *bvdf )

   PARAMETERS:
   .  grid - the grid, which in the blockvector is to be searched
   .  bvd - blockvector description specifying the blockvector
   .  bvdf - blockvector description format to read the 'bvd'

   DESCRIPTION:
   Returns the blockvector specified by 'bvd'.

   RETURN VALUE:
   BLOCKVECTOR *
   .n pointer to the specified blockvector
   .n NULL if the specified blockvector does not exist

   SEE ALSO:
   BLOCKVECTOR, BV_DESC, BVD_PUSH_ENTRY

   D*/
/****************************************************************************/

BLOCKVECTOR *FindBV( const GRID *grid, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf )
{
  register BLOCKVECTOR *bv;
  register BLOCKNUMBER nr;
  BV_DESC copy_bvd = *bvd;

  /* initialize the search with the topmost blockvector in the grid */
  bv = GFIRSTBV( grid );
  BVD_INIT_SEQ_READ( &copy_bvd );
  nr = BVD_READ_NEXT_ENTRY( &copy_bvd, bvdf );

  /* search trough all blockvector levels defined in bvd */
  while ( TRUE )
  {
    while ( BVNUMBER(bv) != nr )             /* search block with number nr in block list*/
    {
      bv = BVSUCC( bv );
      if( bv == NULL )
        return bv;                              /* desired number not found, return NULL */
    }

    /* block with number nr found, prepare search on next block level */
    nr = BVD_READ_NEXT_ENTRY( &copy_bvd, bvdf  );
    if ( nr != NO_BLOCKVECTOR )
      /* there is a further blockvector level specified */
      if ( BV_IS_LEAF_BV( bv ) )
        return NULL;                            /* more block level specified than present */
      else
        bv = BVDOWNBV( bv );                            /* go to next blockvector level */
    else
      /* no more block level specified, current block is the result */
      return bv;
  }
}

/****************************************************************************/
/*D
   GetDomainPart - compute part information of geometrical object

   SYNOPSIS:
   INT GetDomainPart (const INT s2p[], const GEOM_OBJECT *obj, INT side)

   PARAMETERS:
   .  s2p - table translating subdomain to domain part
   .  obj - geometric object (node, element or edge)
   .  side - if element side is meant for obj==element side has to be >=0, negative else

   DESCRIPTION:
   Compute the part info for a geometrical object including element sides.

   RETURN VALUE:
   INT
   .n    part if ok
   .n    -n else
   D*/
/****************************************************************************/

INT GetDomainPart (const INT s2p[], const GEOM_OBJECT *obj, INT side)
{
  NODE *nd,*n0,*n1;
  EDGE *ed;
  ELEMENT *elem;
  VERTEX *v0,*v1;
  BNDS *bs;
  INT part=-1,subdom,move,left,right;

  switch (OBJT(obj))
  {
  case NDOBJ :
    nd = (NODE*)obj;
    v0 = MYVERTEX(nd);
    if (OBJT(v0)==IVOBJ) {
      subdom = NSUBDOM(nd);
      ASSERT(subdom>0);
      part = s2p[subdom];
    }
    else
    {
      /* get part info from domain module */
      if (BNDP_BndPDesc(V_BNDP(v0),&move,&part))
        REP_ERR_RETURN(-2);
      ASSERT(NSUBDOM(nd) == 0);
    }
    break;

  case IEOBJ :
  case BEOBJ :
    elem = (ELEMENT*)obj;
    if (side==NOSIDE)
    {
      /* get info for element */

      subdom = SUBDOMAIN(elem);
      ASSERT(subdom>0);
      part = s2p[subdom];
    }
    else
    {
      /* get info for element side */

      ASSERT(side<SIDES_OF_ELEM(elem));
      ASSERT(side>=0);

      if ((OBJT(elem)==BEOBJ) && ((bs = ELEM_BNDS(elem,side))!=NULL))
      {
        /* this is a boundary side: ask domain module */
        if (BNDS_BndSDesc(ELEM_BNDS(elem,side),&left,&right,&part))
          REP_ERR_RETURN(-3);
      }
      else
      {
        /* there still is the possibility that the side is a boundary side
           because VECTORs are created while CreateElement but boundary sides later
           by CreateSonElementSide.
           The vector eventually will be reallocated by ReinspectSonSideVector later */
        subdom = SUBDOMAIN(elem);
        ASSERT(subdom>0);
        part = s2p[subdom];
      }
    }
    break;

  case EDOBJ :
    ed = (EDGE*)obj;
    n0 = NBNODE(LINK0(ed));
    n1 = NBNODE(LINK1(ed));
    v0 = MYVERTEX(n0);
    v1 = MYVERTEX(n1);
    if ((OBJT(v0)==BVOBJ) && (OBJT(v1)==BVOBJ))
      if (BNDP_BndEDesc(V_BNDP(v0),V_BNDP(v1),&part) == 0)
        return(part);
    subdom = EDSUBDOM(ed);
    if (subdom > 0)
      return(s2p[subdom]);
    subdom = NSUBDOM(n0);
    if (subdom > 0)
      return(s2p[subdom]);
    subdom = NSUBDOM(n1);
    if (subdom > 0)
      return(s2p[subdom]);
    REP_ERR_RETURN(-4);

  default : REP_ERR_RETURN (-5);
  }
  return (part);
}

#ifdef ModelP
INT GetVectorSize (GRID *theGrid, INT VectorObjType, GEOM_OBJECT *object)
{
  MULTIGRID *mg;
  INT part,vtype;

  mg = MYMG(theGrid);
  part = GetDomainPart(BVPD_S2P_PTR(MG_BVPD(mg)),object,NOSIDE);
  if (part < 0)
    REP_ERR_RETURN(-1);
  vtype = FMT_PO2T(MGFORMAT(mg),part,VectorObjType);

  return(FMT_S_VEC_TP(MGFORMAT(mg),vtype));
}
#endif

/****************************************************************************/
/*D
   CreateVectorInPart -  Return pointer to a new vector structure

   SYNOPSIS:
   INT CreateVectorInPart (GRID *theGrid, INT DomPart, INT ObjType, GEOM_OBJECT *object, VECTOR **vHandle);

   PARAMETERS:
   .  theGrid - grid where vector should be inserted
   .  DomPart - part of the domain where vector is created
   .  ObjType - one of the types defined in gm.h
   .  object  - associate vector with this object
   .  vHandle - handle of new vector, i.e. a pointer to a pointer where
                                a pointer to the new vector is placed.

   DESCRIPTION:
   This function returns a pointer to a new vector structure.
   The vector type is determined by 'DomPart' and 'ObjType'
   First the free list is checked for a free entry, if none
   is available, a new structure is allocated from the heap.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occured.
   D*/
/****************************************************************************/

static INT CreateVectorInPart (GRID *theGrid, INT DomPart, INT VectorObjType,
                               GEOM_OBJECT *object, VECTOR **vHandle)
{
  MULTIGRID *theMG;
  VECTOR *pv;
  INT ds, Size, vtype;

  *vHandle = NULL;

  theMG = MYMG(theGrid);
  vtype = FMT_PO2T(MGFORMAT(theMG),DomPart,VectorObjType);
  ds = FMT_S_VEC_TP(MGFORMAT(theMG),vtype);
  if (ds == 0)
    return (0);                         /* HRR: this is ok now, no XXXXVEC in part of the domain */

  Size = sizeof(VECTOR)-sizeof(DOUBLE)+ds;
  pv = (VECTOR *)GetMemoryForObject(theMG,Size,VEOBJ);
  if (pv==NULL)
    REP_ERR_RETURN(1);

  /* initialize data */
  SETOBJT(pv,VEOBJ);
  SETVTYPE(pv,vtype);
  SETVPART(pv,DomPart);
  ds = VPART(pv);
  if (ds!=DomPart)
    return (1);
  SETVDATATYPE(pv,BITWISE_TYPE(vtype));
  SETVOTYPE(pv,VectorObjType);
  SETVCLASS(pv,3);
  SETVNCLASS(pv,0);
  SETVBUILDCON(pv,1);
  SETVNEW(pv,1);
  /*SETPRIO(pv,PrioMaster);*/

    #ifdef __BLOCK_VECTOR_DESC__
  BVD_INIT( &VBVD( pv ) );
    #endif

        #ifdef ModelP
  DDD_AttrSet(PARHDR(pv),GRID_ATTR(theGrid));
        #endif

  VOBJECT(pv) = object;
  VINDEX(pv) = (long)theGrid->nVector;
  SUCCVC(pv) = FIRSTVECTOR(theGrid);

  GRID_LINK_VECTOR(theGrid,pv,PrioMaster);

  *vHandle = pv;

  PRINTDEBUG(gm,1,("%s-vector created (%d): p=%d, t=%d\n",ObjTypeName[VOTYPE(pv)],ID(VOBJECT(pv)),VPART(pv),VTYPE(pv)));

  return (0);
}

/****************************************************************************/
/*D
   CreateVector -  Return pointer to a new vector structure

   SYNOPSIS:
   INT CreateVector (GRID *theGrid, INT VectorObjType, GEOM_OBJECT *object, VECTOR **vHandle);

   PARAMETERS:
   .  theGrid - grid where vector should be inserted
   .  VectorObjType - one of the types defined in gm.h
   .  object  - associate vector with this object
   .  vHandle - handle of new vector, i.e. a pointer to a pointer where
                                a pointer to the new vector is placed.

   DESCRIPTION:
   This function returns a pointer to a new vector structure. First the free list
   is checked for a free entry, if none is available, a new structure is allocated
   from the heap.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occured.
   D*/
/****************************************************************************/

INT CreateVector (GRID *theGrid, INT VectorObjType, GEOM_OBJECT *object, VECTOR **vHandle)
{
  MULTIGRID *mg;
  INT part;

  *vHandle = NULL;
  mg = MYMG(theGrid);
  part = GetDomainPart(BVPD_S2P_PTR(MG_BVPD(mg)),object,NOSIDE);
  if (part < 0)
    REP_ERR_RETURN(1);
  if (CreateVectorInPart(theGrid,part,VectorObjType,object,vHandle)) {
    REP_ERR_RETURN(1);
  }
  else
    return (0);
}

INT CreateSideVector (GRID *theGrid, INT side, GEOM_OBJECT *object, VECTOR **vHandle)
{
  MULTIGRID *mg;
  INT part;

  *vHandle = NULL;
  mg = MYMG(theGrid);
  part = GetDomainPart(BVPD_S2P_PTR(MG_BVPD(mg)),object,side);
  if (part<0)
    REP_ERR_RETURN(1);

  if (CreateVectorInPart(theGrid,part,SIDEVEC,object,vHandle))
    REP_ERR_RETURN(1);

  SETVECTORSIDE(*vHandle,side);
  SETVCOUNT(*vHandle,1);

  return (0);
}

/****************************************************************************/
/*D
   CreateBlockvector -  Return pointer to a new blockvector structure

   SYNOPSIS:
   INT CreateBlockvector( GRID *theGrid, BLOCKVECTOR **BVHandle )

   PARAMETERS:
   .  theGrid - grid where blockvector should be inserted
   .  BVHandle - handle of new blockvector, i.e. a pointer to a pointer where
   a pointer to the new blockvector is placed. Is NULL if there was no more
   memory.

   DESCRIPTION:
   This function returns pointer to a new blockvector structure. First
   the free list is checked for a free entry, if none is available,
   a new structure is allocated from the heap.

   RETURN VALUE:
   INT
   .n      GM_OK if ok
   .n      GM_OUT_OF_MEM if there is not not enough memory available
   D*/
/****************************************************************************/

INT CreateBlockvector( GRID *theGrid, BLOCKVECTOR **BVHandle )
{
  MULTIGRID *theMG;
  BLOCKVECTOR *bv;

  theMG = MYMG(theGrid);

  /* try to find an entry in the free list */
  *BVHandle = NULL;
  bv = (BLOCKVECTOR*)GetMemoryForObject(theMG,sizeof(BLOCKVECTOR),BLOCKVOBJ);
  if (bv == NULL)
    return GM_OUT_OF_MEM;

  SETOBJT(bv,BLOCKVOBJ);
  *BVHandle = bv;

  return GM_OK;
}

static INT InsertBlockvector_l0 (GRID *theGrid, BLOCKVECTOR *insertBV, BLOCKVECTOR *theBV, INT after, INT makeVC)
{
  if (GFIRSTBV(theGrid)==NULL)
  {
    GFIRSTBV(theGrid)=GLASTBV(theGrid)=insertBV;
    BVPRED(insertBV)=BVSUCC(insertBV)=NULL;

    if (makeVC==YES)
    {
      PREDVC(BVFIRSTVECTOR(insertBV)) = NULL;
      SUCCVC(BVLASTVECTOR(insertBV)) = NULL;
    }
  }
  else if (after)
  {
    if (theBV==NULL)
    {
      BVSUCC(insertBV) = GFIRSTBV(theGrid);
      BVPRED(GFIRSTBV(theGrid)) = insertBV;
      GFIRSTBV(theGrid) = insertBV;
      BVPRED(insertBV) = NULL;

      if (makeVC==YES)
      {
        PREDVC(BVFIRSTVECTOR(insertBV)) = NULL;
        SUCCVC(BVLASTVECTOR(insertBV)) = BVFIRSTVECTOR(BVSUCC(insertBV));
        PREDVC(BVFIRSTVECTOR(BVSUCC(insertBV))) = BVLASTVECTOR(insertBV);
        SFIRSTVECTOR(theGrid) = BVFIRSTVECTOR(insertBV);
      }
    }
    else
    {
      BVPRED(insertBV) = theBV;
      BVSUCC(insertBV) = BVSUCC(theBV);
      if (theBV==GLASTBV(theGrid)) GLASTBV(theGrid)=insertBV;
      else BVPRED(BVSUCC(theBV))=insertBV;
      BVSUCC(theBV) = insertBV;

      if (makeVC==YES)
      {
        SUCCVC(BVLASTVECTOR(insertBV)) = SUCCVC(BVLASTVECTOR(theBV));
        if (SUCCVC(BVLASTVECTOR(insertBV))!=NULL) PREDVC(SUCCVC(BVLASTVECTOR(insertBV))) = BVLASTVECTOR(insertBV);
        else LASTVECTOR(theGrid) = BVLASTVECTOR(insertBV);
        SUCCVC(BVLASTVECTOR(theBV)) = BVFIRSTVECTOR(insertBV);
        PREDVC(BVFIRSTVECTOR(insertBV)) = BVLASTVECTOR(theBV);
      }
    }
  }
  else
  {
    if (theBV==NULL)
      return (InsertBlockvector_l0(theGrid,insertBV,GLASTBV(theGrid),1,makeVC));
    else
      return (InsertBlockvector_l0(theGrid,insertBV,BVPRED(theBV),1,makeVC));
  }

  return (GM_OK);
}

INT CreateBlockvector_l0 (GRID *theGrid, BLOCKVECTOR **BVHandle, BLOCKVECTOR *insertBV, INT after)
{
  BLOCKVECTOR *theBV;

  if (CreateBlockvector (theGrid, &theBV)!=GM_OK) return (GM_OUT_OF_MEM);
  if (InsertBlockvector_l0(theGrid,theBV,insertBV,after,NO)) return (GM_OUT_OF_MEM);

  *BVHandle = theBV;
  return (GM_OK);
}

static INT CutBlockvector_l0 (GRID *theGrid, BLOCKVECTOR *theBV, INT makeVC)
{
  if (theBV==NULL) return (GM_ERROR);

  switch ((theBV==GFIRSTBV(theGrid)) | ((theBV==GLASTBV(theGrid))<<1))
  {
  case 0 :
    BVSUCC(BVPRED(theBV)) = BVSUCC(theBV);
    BVPRED(BVSUCC(theBV)) = BVPRED(theBV);

    if (makeVC==YES)
    {
      SUCCVC(BVLASTVECTOR(BVPRED(theBV))) = BVFIRSTVECTOR(BVSUCC(theBV));
      PREDVC(BVFIRSTVECTOR(BVSUCC(theBV))) = BVLASTVECTOR(BVPRED(theBV));
    }
    break;
  case 1 :
    BVPRED(BVSUCC(theBV)) = NULL;
    GFIRSTBV(theGrid) = BVSUCC(theBV);

    if (makeVC==YES)
      PREDVC(BVFIRSTVECTOR(GFIRSTBV(theGrid))) = NULL;
    break;
  case 2 :
    BVSUCC(BVPRED(theBV)) = NULL;
    GLASTBV(theGrid) = BVPRED(theBV);

    if (makeVC==YES)
      SUCCVC(BVLASTVECTOR(BVPRED(theBV))) = NULL;
    break;
  case 3 :
    GFIRSTBV(theGrid) = GLASTBV(theGrid) = NULL;
    break;
  }

  return (GM_OK);
}

/****************************************************************************/
/*D
   CreateConnection -  Return pointer to a new connection structure

   SYNOPSIS:
   CONNECTION *CreateConnection (GRID *theGrid, VECTOR *from, VECTOR *to);

   PARAMETERS:
   .  TheGrid - grid where matrix should be inserted
   .  from - source vector
   .  to - destination vector

   DESCRIPTION:
   This function allocates a new 'CONNECTION' and inserts the two
   'MATRIX' structures in the lists of 'from' and 'to' vector.
   Since the operation is symmetric, the order of 'from' and 'to'
   is not important.

   RETURN VALUE:
   INT
   .n    NULL if error occured.
   .n    else a pointer to the new 'CONNECTION' is returned.
   D*/
/****************************************************************************/

CONNECTION *CreateConnection (GRID *theGrid, VECTOR *from, VECTOR *to)
{
  MULTIGRID *theMG;
  HEAP *theHeap;
  CONNECTION *pc;
  MATRIX *pm;
  INT RootType, DestType, MType, ds, Diag, Size;

  /* set Diag, RootType and DestType	*/
  Diag = ((from == to) ? 1 : 0);
  RootType = VTYPE(from);
  DestType = VTYPE(to);
  MType = MatrixType[RootType][DestType];

  /* check expected size */
  theMG = MYMG(theGrid);
  theHeap = theMG->theHeap;
  ds = FMT_S_MAT_TP(MGFORMAT(theMG),MType);
  if (ds == 0)
    return (NULL);
  Size = sizeof(MATRIX)-sizeof(DOUBLE)+ds;
  if (MSIZEMAX<Size) return (NULL);

  /* is there already the desired connection ? */
  pc = GetConnection(from,to);
  if (pc != NULL)
  {
    SETCEXTRA(pc,0);
    return (pc);
  }

  if (Diag)
    pc = (CONNECTION*)GetMemoryForObject(theMG,Size,MAOBJ);
  else
    pc = (CONNECTION*)GetMemoryForObject(theMG,2*Size,COOBJ);
  if (pc==NULL) return (NULL);

  /* initialize data */
  pm = CMATRIX0(pc);
  SETOBJT(pm,MAOBJ);
  SETMROOTTYPE(pm,RootType);
  SETMDESTTYPE(pm,DestType);
  SETMDIAG(pm,Diag);
  SETMOFFSET(pm,0);
  SETMSIZE(pm,Size);
  SETMNEW(pm,1);
  MDEST(pm) = to;
  if (!Diag)
  {
    pm = CMATRIX1(pc);
    CTRL(pm) = 0;
    SETOBJT(pm,MAOBJ);
    SETMROOTTYPE(pm,DestType);
    SETMDESTTYPE(pm,RootType);
    SETMDIAG(pm,Diag);
    SETMOFFSET(pm,1);
    SETMSIZE(pm,Size);
    SETMNEW(pm,1);
    MDEST(pm) = from;
  }

  /* set sizes */
  if (!Diag)
  {
    Size = (unsigned INT)pm - (unsigned INT)pc;
    SETMSIZE(pc,Size);
    SETMSIZE(pm,Size);
  }

  /* put in matrix list */
  if (Diag)
  {
    /* insert at first place in the list (only one matrix) */
    MNEXT(CMATRIX0(pc)) = VSTART(from);
    VSTART(from) = CMATRIX0(pc);
  }
  else
  {
    /* insert at second place in the list (both matrices) */
    pm = VSTART(from);
    if (pm == NULL)
    {
      MNEXT(CMATRIX0(pc)) = NULL;
      VSTART(from) = CMATRIX0(pc);
    }
    else
    {
      MNEXT(CMATRIX0(pc)) = MNEXT(pm);
      MNEXT(pm) = CMATRIX0(pc);
    }

    pm = VSTART(to);
    if (pm == NULL)
    {
      MNEXT(CMATRIX1(pc)) = NULL;
      VSTART(to) = CMATRIX1(pc);
    }
    else
    {
      MNEXT(CMATRIX1(pc)) = MNEXT(pm);
      MNEXT(pm) = CMATRIX1(pc);
    }
  }

  /* counters */
  theGrid->nCon++;

  return(pc);
}

/****************************************************************************/
/*D
   CreateExtraConnection -  Return pointer to a new matrix structure with extra flag set.

   SYNOPSIS:
   CONNECTION *CreateExtraConnection (GRID *theGrid, VECTOR *from,
   VECTOR *to);

   PARAMETERS:
   .  theGrid - grid level where connection will be inserted.
   .  from,to - Pointers to vectors where connection is inserted.

   DESCRIPTION:
   This function returns a pointer to a new 'CONNECTION'
   structure with extra flag set. This e.g. for a direct solver
   or ILU with fill in. The new connections can be distinguished
   from the connections necessary for the stiffness matrix.

   RETURN VALUE:
   INT
   .n    NULL if error occured.
   .n    else pointer to new 'CONNECTION'
   D*/
/****************************************************************************/

CONNECTION      *CreateExtraConnection  (GRID *theGrid, VECTOR *from, VECTOR *to)
{
  CONNECTION *pc;

  pc = CreateConnection(theGrid,from,to);
  if (pc==NULL) return(NULL);
  SETCEXTRA(pc,1);
  return(pc);
}

INT CreateElementList (GRID *theGrid, NODE *theNode, ELEMENT *theElement)
{
  ELEMENTLIST *pel;

  for (pel=NODE_ELEMENT_LIST(theNode); pel!=NULL; pel=NEXT(pel))
    if(pel->el==theElement)
      return(0);

  pel = (ELEMENTLIST *)GetMemoryForObject(theGrid->mg,sizeof(ELEMENTLIST),-1);
  if (pel == NULL)
    return(1);

  pel->next = NODE_ELEMENT_LIST(theNode);
  pel->el = theElement;
  NDATA(theNode) = (void *) pel;

  return(0);
}

/****************************************************************************/
/*D
   DisposeVector - Remove vector from the data structure

   SYNOPSIS:
   INT DisposeVector (GRID *theGrid, VECTOR *theVector);

   PARAMETERS:
   .  theGrid - grid level where theVector is in.
   .  theVector - 'VECTOR' to be disposed.

   DESCRIPTION:
   This function removes vector from the data structure and places
   it in the free list.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT DisposeVector (GRID *theGrid, VECTOR *theVector)
{
  MATRIX *theMatrix, *next;
  INT Size;

  if (theVector == NULL)
    return(0);

  HEAPFAULT(theVector);

  /* remove all connections concerning the vector */
  for (theMatrix=VSTART(theVector); theMatrix!=NULL; theMatrix=next)
  {
    next = MNEXT(theMatrix);
    if (DisposeConnection(theGrid,MMYCON(theMatrix)))
      RETURN (1);
  }

#ifdef __INTERPOLATION_MATRIX__
  if (DisposeIMatrices(theGrid,VISTART(theVector)))
    RETURN (1);
#endif

  /* now remove vector from vector list */
  GRID_UNLINK_VECTOR(theGrid,theVector);

  /* reset count flags */
  SETVCOUNT(theVector,0);


  /* delete the vector itself */
  Size = sizeof(VECTOR)-sizeof(DOUBLE)
         + FMT_S_VEC_TP(MGFORMAT(MYMG(theGrid)),VTYPE(theVector));
  if (PutFreeObject(theGrid->mg,theVector,Size,VEOBJ))
    RETURN(1);

  return(0);
}

/****************************************************************************/
/*D
   ReallocateVector - change vector allocated with wrong part

   SYNOPSIS:
   INT ReinspectSonSideVector (GRID *g, ELEMENT *elem, INT side, VECTOR **vHandle)

   PARAMETERS:
   .  g - grid level where the element is in.
   .  elem - element of side vector
   .  side - element side
   .  vHandle - handle to side vector (inialized with old, may be changed)

   DESCRIPTION:
   This changes a side vector which was allocated with the wrong part (maybe unknown for side vectors
   when creating an element).

   RETURN VALUE:
   INT
   .n    GM_OK if ok
   .n    GM_ERROR if error occured.
   D*/
/****************************************************************************/

INT ReinspectSonSideVector (GRID *g, ELEMENT *elem, INT side, VECTOR **vHandle)
{
  MULTIGRID *mg;
  VECTOR *vold,*vnew;
  FORMAT *fmt;
  INT partnew,partold,vtnew,vtold,dsnew,dsold;

  mg  = MYMG(g);
  fmt = MGFORMAT(mg);

  vold = *vHandle;

  /* check whether part has actually changed */
  partold = (vold!=NULL) ? VPART(vold) : BVPD_S2P(MG_BVPD(mg),SUBDOMAIN(elem));
  partnew = GetDomainPart(BVPD_S2P_PTR(MG_BVPD(mg)),(GEOM_OBJECT*)elem,side);
  if (partnew<0)
    REP_ERR_RETURN(GM_ERROR);
  if (partnew==partold)
    return (GM_OK);

  /* check whether vtype has actually changed */
  vtold = (vold!=NULL) ? VTYPE(vold) : FMT_PO2T(fmt,partold,SIDEVEC);
  vtnew = FMT_PO2T(fmt,partnew,SIDEVEC);
  if (vtnew==vtold)
  {
    if (vold!=NULL)
    {
      /* just change part */
      SETVPART(vold,partnew);
    }
    PRINTDEBUG(gm,1,("SIDEVEC (%d,%d): part\n",ID(elem),side));
    return (GM_OK);
  }

  /* check whether size has actually changed */
  dsold = FMT_S_VEC_TP(fmt,vtold);
  dsnew = FMT_S_VEC_TP(fmt,vtnew);
  if (dsold==dsnew)
  {
    if (vold!=NULL)
    {
      /* just change part and type */
      SETVTYPE(vold,vtnew);
      SETVPART(vold,partnew);

      DisposeConnectionFromVector(g,vold);

      SETVBUILDCON(vold,1);
    }
    PRINTDEBUG(gm,1,("SIDEVEC (%d,%d): part and type\n",ID(elem),side));
    return (GM_OK);
  }

  PRINTDEBUG(gm,1,("SIDEVEC (%d,%d): part, type and size\n",ID(elem),side));

  /* create new vector */
  if (CreateVectorInPart(g,partnew,SIDEVEC,(GEOM_OBJECT*)elem,&vnew))
    REP_ERR_RETURN(GM_ERROR);

  /* now we can dispose the old vector */
  if (DisposeVector(g,vold))
    REP_ERR_RETURN(GM_ERROR);

  *vHandle = vnew;

  return (GM_OK);
}

/****************************************************************************/
/*D
   DisposeBlockvector - Dispose blockvector

   SYNOPSIS:
   INT DisposeBlockvector( GRID *theGrid, BLOCKVECTOR *bv )

   PARAMETERS:
   .  theGrid - grid level where bv is in.
   .  bv - 'BLOCKVECTOR' to be disposed.

   DESCRIPTION:
   This function places the blockvector in the free list. The data structure
   must be kept consistent by the caller! If bv is the NULL-vector nothing
   will happen.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT DisposeBlockvector( GRID *theGrid, BLOCKVECTOR *bv )
{
  if ( bv == NULL ) return 0;
  return PutFreeObject( MYMG(theGrid), bv,sizeof(BLOCKVECTOR),BLOCKVOBJ);
}

/****************************************************************************/
/*D
   FreeBVList - Frees blockvector-list and all its sons

   SYNOPSIS:
   static void FreeBVList (GRID *grid, BLOCKVECTOR *bv);

   PARAMETERS:
   .  grid - the grid from which the blockvectors are to be removed
   .  bv - start of the blockvector-list to be removed

   DESCRIPTION:
   Frees double linked blockvector-list, starting with 'bv', and all
   its sons recursively, i.e. all the
   subblockvectors. The vectors belonging to the blockvectors
   are not influenced.

   RETURN VALUE:
   void

   SEE ALSO:
   BLOCKVECTOR, FreeAllBV

   D*/
/****************************************************************************/

static void FreeBVList (GRID *grid, BLOCKVECTOR *bv)
{
  register BLOCKVECTOR *bv_h;

  if ( bv == NULL )
    return;

  /* free all blocks in the current block list */
  for ( ; bv != NULL; )
  {
    /* free rekursively all higher block levels */
    if ( !BV_IS_LEAF_BV( bv ) )
      FreeBVList( grid, BVDOWNBV( bv ) );

    /* free the current block and go to the next */
    bv_h = BVSUCC( bv );
    DisposeBlockvector( grid, bv );
    bv = bv_h;
  }
}

/****************************************************************************/
/*D
   FreeAllBV - Frees all allocated blockvectors in the grid

   SYNOPSIS:
   void FreeAllBV( GRID *grid )

   PARAMETERS:
   .  grid - the grid from which all blockvectors are to be removed

   DESCRIPTION:
   Frees all allocated blockvectors in the grid and resets GFIRSTBV(grid)
   and GLASTBV(grid) to 'NULL' and gives the memory back to the memory
   of the 'MULTIGRID'.

   RETURN VALUE:
   void

   SEE ALSO:
   BLOCKVECTOR, FreeBVList

   D*/
/****************************************************************************/

void FreeAllBV (GRID *grid)
{
  FreeBVList( grid, GFIRSTBV( grid ) );
  GFIRSTBV( grid ) = NULL;
  GLASTBV( grid ) = NULL;
}

/****************************************************************************/
/*D
   DisposeConnection - Remove connection from the data structure

   SYNOPSIS:
   INT DisposeConnection (GRID *theGrid, CONNECTION *theConnection);

   PARAMETERS:
   .  theGrid - the grid to remove from
   .  theConnection - connection to dispose

   DESCRIPTION:
   This function removes a connection from the data structure. The connection
   is removed from the list of the two vectors and is placed in the
   free list.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT DisposeConnection (GRID *theGrid, CONNECTION *theConnection)
{
  VECTOR *from, *to;
  MATRIX *Matrix, *ReverseMatrix, *SearchMatrix;

  HEAPFAULT(theConnection);

  /* remove matrix(s) from their list(s) */
  Matrix = CMATRIX0(theConnection);
  to = MDEST(Matrix);
  if (MDIAG(Matrix))
  {
    from = to;
    VSTART(to) = MNEXT(Matrix);
  }
  else
  {
    ReverseMatrix = CMATRIX1(theConnection);
    from = MDEST(ReverseMatrix);
    if (VSTART(from) == Matrix)
      VSTART(from) = MNEXT(Matrix);
    else
      for (SearchMatrix=VSTART(from); SearchMatrix!=NULL; SearchMatrix=MNEXT(SearchMatrix))
        if (MNEXT(SearchMatrix) == Matrix)
          MNEXT(SearchMatrix) = MNEXT(Matrix);
    if (VSTART(to) == ReverseMatrix)
      VSTART(to) = MNEXT(ReverseMatrix);
    else
      for (SearchMatrix=VSTART(to); SearchMatrix!=NULL; SearchMatrix=MNEXT(SearchMatrix))
        if (MNEXT(SearchMatrix) == ReverseMatrix)
          MNEXT(SearchMatrix) = MNEXT(ReverseMatrix);
  }

  /* free connection object */
  if (MDIAG(Matrix))
    PutFreeObject(MYMG(theGrid),Matrix,MSIZE(Matrix),MAOBJ);
  else
    PutFreeObject(MYMG(theGrid),Matrix,2*MSIZE(Matrix),COOBJ);

  theGrid->nCon--;

  /* return ok */
  return(0);
}

/****************************************************************************/
/*D
   DisposeDoubledSideVector - Dispose vector associated with two sides in 3D

   SYNOPSIS:
   INT DisposeDoubledSideVector (GRID *theGrid, ELEMENT *Elem0, INT Side0,
   ELEMENT *Elem1, INT Side1);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  Elem0,side0 - first element and side
   .  Elem1,side1 - second element and side

   DESCRIPTION:
   This function disposes a vector associated with two sides in 3D.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

#ifdef __THREEDIM__
INT DisposeDoubledSideVector (GRID *theGrid, ELEMENT *Elem0, INT Side0, ELEMENT *Elem1, INT Side1)
{
  VECTOR *Vector0, *Vector1;
  int vc0,vc1;

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
  {
    assert(NBELEM(Elem0,Side0)==Elem1 && NBELEM(Elem1,Side1)==Elem0);
    Vector0 = SVECTOR(Elem0,Side0);
    Vector1 = SVECTOR(Elem1,Side1);
    if (Vector0 == Vector1)
      return (0);
    if ((Vector0==NULL) || (Vector1==NULL))
      /* this is the case at boundaries between different parts
         the part not using vectors in SIDEs will not need a pointer
         to the side vector */
      return (0);
    vc0=VCOUNT(Vector0);
    vc1=VCOUNT(Vector1);
    assert(VCOUNT(Vector0)==1 && VCOUNT(Vector1)==1);
    assert(VSTART(Vector0)==NULL || VSTART(Vector1)==NULL);
    if (VSTART(Vector0)==NULL)
    {
      SET_SVECTOR(Elem0,Side0,Vector1);
      SETVCOUNT(Vector1,2);
      if (DisposeVector (theGrid,Vector0))
        RETURN (1);
    }
    else
    {
      SET_SVECTOR(Elem1,Side1,Vector0);
      SETVCOUNT(Vector0,2);
      if (DisposeVector (theGrid,Vector1))
        RETURN (1);
    }
    return (0);

  }

  RETURN (1);
}
#endif

/****************************************************************************/
/*D
   DisposeConnectionFromVector - Remove all connections associated with a vector

   SYNOPSIS:
   INT DisposeConnectionFromVector (GRID *theGrid, VECTOR *theVector);

   PARAMETERS:
   .  theGrid - grid level where vector belongs to
   .  theVector - vector where connections are disposed from

   DESCRIPTION:
   This function removes all connections from a vector.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT DisposeConnectionFromVector (GRID *theGrid, VECTOR *theVector)
{
  while(VSTART(theVector) != NULL)
    if (DisposeConnection (theGrid,MMYCON(VSTART(theVector))))
      return (1);

  return (0);
}


/****************************************************************************/
/*D
   DisposeConnectionFromElement	- Removes all connections from all vectors
   associated with an element

   SYNOPSIS:
   INT DisposeConnectionFromElement (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - grid level where element is on
   .  theElement - element from which to dispose connections

   DESCRIPTION:
   This function removes all connections from all vectors
   associated with an element.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT DisposeConnectionFromElement (GRID *theGrid, ELEMENT *theElement)
{
  INT i;
  VECTOR *vList[20];
  INT cnt;

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
  {
    GetVectorsOfElement(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
    {
      if (DisposeConnectionFromVector(theGrid,vList[i])) RETURN(GM_ERROR);
      SETVBUILDCON(vList[i],1);
    }
  }
    #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
  {
    GetVectorsOfSides(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
    {
      if (DisposeConnectionFromVector(theGrid,vList[i])) RETURN(GM_ERROR);
      SETVBUILDCON(vList[i],1);
    }
  }
    #endif
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
  {
    GetVectorsOfEdges(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
    {
      if (DisposeConnectionFromVector(theGrid,vList[i])) RETURN(GM_ERROR);
      SETVBUILDCON(vList[i],1);
    }
  }
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
  {
    GetVectorsOfNodes(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
    {
      if (DisposeConnectionFromVector(theGrid,vList[i])) RETURN(GM_ERROR);
      SETVBUILDCON(vList[i],1);
    }
  }

  return(GM_OK);
}

/****************************************************************************/
/*D
   DisposeConnectionFromElementInNeighborhood -  Remove matrices

   SYNOPSIS:
   static INT DisposeConnectionFromElementInNeighborhood (GRID *theGrid, ELEMENT *theElement,
   INT Depth);

   PARAMETERS:
   .  theGrid - the grid to remove from
   .  theElement - that element
   .  Depth -  that many slices around the element

   DESCRIPTION:
   This function removes connections concerning an element from the data structure
   and stores flags saying: "connection has to be rebuild",
   it does this in a neighborhood of the elem of depth Depth, where depth
   is the distance in the element-neighborship-graph (see also 'FORMAT').

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT DisposeConnectionFromElementInNeighborhood (GRID *theGrid, ELEMENT *theElement, INT Depth)
{
  INT i;

  if (Depth < 0) RETURN (GM_ERROR);

  /* create connection at that depth */
  if(!EBUILDCON(theElement))
    if (DisposeConnectionFromElement(theGrid,theElement))
      RETURN (GM_ERROR);
  SETEBUILDCON(theElement,1);

  /* dispose connection in neighborhood */
  if (Depth > 0)
  {
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if (DisposeConnectionFromElementInNeighborhood(theGrid,NBELEM(theElement,i),Depth-1))
        RETURN (GM_ERROR);
  }

  RETURN (GM_OK);
}

INT DisposeConnectionsInNeighborhood (GRID *theGrid, ELEMENT *theElement)
{
  INT Depth;
  Depth = (INT)(floor(0.5*(double)FMT_CONN_DEPTH_MAX(MGFORMAT(MYMG(theGrid)))));
  return(DisposeConnectionFromElementInNeighborhood(theGrid,theElement,Depth));
}

INT     DisposeElementFromElementList (GRID *theGrid, NODE *theNode,
                                       ELEMENT *theElement)
{
  ELEMENTLIST *pel,*next;

  pel = NODE_ELEMENT_LIST(theNode);
  if (pel == NULL) return(0);
  if (pel->el == theElement) {
    NDATA(theNode) = (void *) pel->next;
    return(PutFreeObject(theGrid->mg,pel,sizeof(ELEMENTLIST),-1));
  }
  next = pel->next;
  while (next != NULL) {
    if (next->el == theElement) {
      pel->next = next->next;
      return(PutFreeObject(theGrid->mg,next,sizeof(ELEMENTLIST),-1));
    }
    pel = next;
    next = pel->next;
  }

  return(0);
}

INT     DisposeElementList (GRID *theGrid, NODE *theNode)
{
  ELEMENTLIST *pel,*next;

  pel = NODE_ELEMENT_LIST(theNode);
  while (pel != NULL) {
    next = pel->next;
    if (PutFreeObject(theGrid->mg,pel,sizeof(ELEMENTLIST),-1))
      return(1);
    pel = next;
  }

  return(0);
}

/****************************************************************************/
/*D
   GetMatrix - Return pointer to matrix if it exists

   SYNOPSIS:
   MATRIX *GetMatrix (const VECTOR *FromVector, const VECTOR *ToVector);

   PARAMETERS:
   .  FromVector - starting vector of the Matrix
   .  ToVector - destination vector of the Matrix

   DESCRIPTION:
   This function returns pointer to matrix if it exists. The function
   runs through the single linked list, since the list is
   assumed to be small (sparse matrix!) the cost is assumed to be negligible.

   RETURN VALUE:
   MATRIX *
   .n       pointer to Matrix,
   .n       NULL if Matrix does not exist.
   D*/
/****************************************************************************/

MATRIX *GetMatrix (const VECTOR *FromVector, const VECTOR *ToVector)
{
  MATRIX *theMatrix;

  for (theMatrix=VSTART(FromVector); theMatrix!=NULL; theMatrix = MNEXT(theMatrix))
    if (MDEST(theMatrix)==ToVector)
      return (theMatrix);

  /* return not found */
  return (NULL);
}

MATRIX *GetOrderedMatrix (const VECTOR *FromVector, const VECTOR *ToVector)
{
  MATRIX *theMatrix;

  if (FromVector == ToVector)
    return(VSTART(ToVector));
  else if (VINDEX(FromVector) > VINDEX(ToVector)) {
    for (theMatrix=MNEXT(VSTART(FromVector));
         theMatrix!=NULL; theMatrix = MNEXT(theMatrix))
      if (MDEST(theMatrix)==ToVector)
        return (theMatrix);
  }
  else {
    for (theMatrix=MNEXT(VSTART(ToVector));
         theMatrix!=NULL; theMatrix = MNEXT(theMatrix))
      if (MDEST(theMatrix)==FromVector)
        return (MADJ(theMatrix));
  }

  /* return not found */
  return (NULL);
}

/****************************************************************************/
/*D
   GetConnection - Return pointer to connection if it exists

   SYNOPSIS:
   CONNECTION *GetConnection (const VECTOR *FromVector, const VECTOR *ToVector);

   PARAMETERS:
   .  FromVector - starting vector of the con
   .  ToVector - destination vector of the con

   DESCRIPTION:
   This function returns pointer to connection if it exists.

   RETURN VALUE:
   CONNECTION *
   .n       pointer to
   .n       NULL if connection does not exist.
   D*/
/****************************************************************************/

CONNECTION *GetConnection (const VECTOR *FromVector, const VECTOR *ToVector)
{
  MATRIX *Matrix;

  Matrix = GetMatrix(FromVector,ToVector);
  if (Matrix != NULL)
    return (MMYCON(Matrix));

  /* return not found */
  return (NULL);
}

/****************************************************************************/
/*D
   GetVectorsOfElement - Get a pointer list to all element data

   SYNOPSIS:
   INT GetVectorsOfElement (const ELEMENT *theElement, INT *cnt, VECTOR **vList);

   PARAMETERS:
   .  theElement - that element
   .  cnt - how many vectors
   .  vList - array to store vector list

   DESCRIPTION:
   This function returns a pointer to the 'VECTOR' associated with the
   element (if the element is allowed to have one). 'cnt' will either
   be 0 or 1.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

INT GetVectorsOfElement (const ELEMENT *theElement, INT *cnt, VECTOR **vList)
{
  *cnt = 0;
  if (EVECTOR(theElement) != NULL)
    vList[(*cnt)++] = EVECTOR(theElement);

  IFDEBUG(gm,0)
  if (*cnt)
  {
    assert(vList[0] != NULL);
    assert(VOTYPE(vList[0]) == ELEMVEC);
  }
  ENDDEBUG

  return(GM_OK);
}

/****************************************************************************/
/*D
   GetVectorsOfSides - Get a pointer list to all side data

   SYNOPSIS:
   INT GetVectorsOfSides (const ELEMENT *theElement, INT *cnt, VECTOR **vList);

   PARAMETERS:
   .  theElement - that element
   .  cnt - how many vectors
   .  vList - array to store vector list

   DESCRIPTION:
   This function gets a pointer array to all 'VECTOR's in sides of the given element.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

#ifdef __THREEDIM__
INT GetVectorsOfSides (const ELEMENT *theElement, INT *cnt, VECTOR **vList)
{
  INT i;

  *cnt = 0;

  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    if (SVECTOR(theElement,i) != NULL)
      vList[(*cnt)++] = SVECTOR(theElement,i);

  IFDEBUG(gm,0)
  for (i=0; i<(*cnt); i++)
  {
    assert(vList[i] != NULL);
    assert(VOTYPE(vList[i]) == SIDEVEC);
  }
  ENDDEBUG

  return(GM_OK);
}
#endif

/****************************************************************************/
/*D
   GetVectorsOfEdges - Get a pointer list to all edge data

   SYNOPSIS:
   INT GetVectorsOfEdges (const ELEMENT *theElement, INT *cnt, VECTOR **vList);

   PARAMETERS:
   .  theElement -  that element
   .  cnt - how many vectors
   .  vList - array to store vector list

   DESCRIPTION:
   This function gets a pointer array to all 'VECTOR's in edges of the element.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

INT GetVectorsOfEdges (const ELEMENT *theElement, INT *cnt, VECTOR **vList)
{
  EDGE *theEdge;
  INT i;

  *cnt = 0;
  for (i=0; i<EDGES_OF_ELEM(theElement); i++)
    if ((theEdge =
           GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
                   CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)))) != NULL)
      if (EDVECTOR(theEdge) != NULL)
        vList[(*cnt)++] = EDVECTOR(theEdge);

  IFDEBUG(gm,0)
  for (i=0; i<(*cnt); i++)
  {
    assert(vList[i] != NULL);
    assert(VOTYPE(vList[i]) == EDGEVEC);
  }
  ENDDEBUG

  return(GM_OK);
}

/****************************************************************************/
/*D
   GetVectorsOfNodes - Get a pointer list to all node data

   SYNOPSIS:
   INT GetVectorsOfNodes (const ELEMENT *theElement, INT *cnt, VECTOR **vList);

   PARAMETERS:
   .  theElement -  that element
   .  cnt - how many vectors
   .  vList - array to store vector list

   DESCRIPTION:
   This function gets a pointer array to all 'VECTOR's in nodes of the element.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

INT GetVectorsOfNodes (const ELEMENT *theElement, INT *cnt, VECTOR **vList)
{
  INT i;

  *cnt = 0;
  for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    if (NVECTOR(CORNER(theElement,i)) != NULL)
      vList[(*cnt)++] = NVECTOR(CORNER(theElement,i));

  IFDEBUG(gm,0)
  for (i=0; i<(*cnt); i++)
  {
    assert(vList[i] != NULL);
    assert(VOTYPE(vList[i]) == NODEVEC);
  }
  ENDDEBUG

  return (GM_OK);
}

/****************************************************************************/
/*D
   GetVectorsOfOType - Get a pointer list to all vector data of specified object type

   SYNOPSIS:
   INT GetVectorsOfOType (const ELEMENT *theElement, INT type, INT *cnt, VECTOR **vList)

   PARAMETERS:
   .  theElement -  that element
   .  type - fill array with vectors of this type
   .  cnt - how many vectors
   .  vList - array to store vector list

   DESCRIPTION:
   This function returns a pointer array to all 'VECTOR's of the element of the specified type.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

INT GetVectorsOfOType (const ELEMENT *theElement, INT type, INT *cnt, VECTOR **vList)
{
  switch (type)
  {
  case NODEVEC :
    return (GetVectorsOfNodes(theElement,cnt,vList));
  case ELEMVEC :
    return (GetVectorsOfElement(theElement,cnt,vList));
  case EDGEVEC :
    return (GetVectorsOfEdges(theElement,cnt,vList));
        #ifdef __THREEDIM__
  case SIDEVEC :
    return (GetVectorsOfSides(theElement,cnt,vList));
        #endif
  }
  RETURN (GM_ERROR);
}

/****************************************************************************/
/*D
   DataTypeFilterVList - remove vectors with non-matching type from list

   SYNOPSIS:
   INT DataTypeFilterVList (INT dt, VECTOR **vec, INT *cnt)

   PARAMETERS:
   .  dt - data type including all vtypes needed
   .  vec - vector list
   .  cnt  - number of vectors return in the list

   DESCRIPTION:
   This function removes vectors with non-matching types from the list.

   RETURN VALUE:
   INT
   .n    number of remaining vectors in the list
   D*/
/****************************************************************************/

INT DataTypeFilterVList (INT dt, VECTOR **vec, INT *cnt)
{
  INT i,n;

  n = *cnt;
  *cnt = 0;
  for (i=0; i<n; i++)
    if (VDATATYPE(vec[i]) & dt)
      vec[(*cnt)++] = vec[i];

  return (*cnt);
}

/****************************************************************************/
/*D
   GetVectorsOfDataTypesInObjects - get vector list including all vectors of specified vtypes

   SYNOPSIS:
   INT GetVectorsOfDataTypesInObjects (const ELEMENT *theElement, INT dt, INT obj, INT *cnt, VECTOR *VecList[])

   PARAMETERS:
   .  theElement - pointer to an element
   .  dt  - data type including all vtypes needed
   .  obj - flags for objects types needed
   .  cnt - number of vectors return in the list
   .  vec - vector list

   DESCRIPTION:
   This function gets a list of vectors of the specified vtypes corresponding to an element.

   RETURN VALUE:
   INT
   .n    GM_OK if ok
   .n    GM_ERROR if error occured
   D*/
/****************************************************************************/

INT GetVectorsOfDataTypesInObjects (const ELEMENT *theElement, INT dt, INT obj, INT *cnt, VECTOR *VecList[])
{
  INT i,n;

  *cnt = n = 0;

  if (obj & BITWISE_TYPE(NODEVEC))
  {
    if (GetVectorsOfNodes(theElement,&i,VecList) != GM_OK)
      return(GM_ERROR);
    n += i;
  }

  if (obj & BITWISE_TYPE(EDGEVEC))
  {
    if (GetVectorsOfEdges(theElement,&i,VecList+n) != GM_OK)
      return(GM_ERROR);
    n += i;
  }

  if (obj & BITWISE_TYPE(ELEMVEC))
  {
    if (GetVectorsOfElement(theElement,&i,VecList+n) != GM_OK)
      return(GM_ERROR);
    n += i;
  }

    #ifdef __THREEDIM__
  if (obj & BITWISE_TYPE(SIDEVEC))
  {
    if (GetVectorsOfSides(theElement,&i,VecList+n) != GM_OK)
      return(GM_ERROR);
    n += i;
  }
    #endif

  *cnt = n;

  /* remove vectors of types not needed */
  DataTypeFilterVList(dt,VecList,cnt);

  return (GM_OK);
}

/****************************************************************************/
/*D
   PrepareGetBoundaryNeighbourVectors - prepare 'GetBoundaryNeighbourVectors'

   SYNOPSIS:
   INT PrepareGetBoundaryNeighbourVectors (GRID *theGrid, INT *MaxListLen)

   PARAMETERS:
   .  theGrid - grid level
   .  MaxListLen - max size of neighbourhood and therefore max list lenght of the
                                VecList array of 'GetBoundaryNeighbourVectors'

   DESCRIPTION:
   This function stores lists of boundary vector neighbourhoods in temp storage
   of the multigrid. The neighborhoods of the boundary vectors can be received
   via a call of 'GetBoundaryNeighbourVectors' one by one.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 else

   SEE ALSO:
   'GetBoundaryNeighbourVectors', 'ResetGetBoundaryNeighbourVectors', 'FinishBoundaryNeighbourVectors'
   D*/
/****************************************************************************/

INT PrepareGetBoundaryNeighbourVectors (GRID *theGrid, INT *MaxListLen)
{
#ifdef __TWODIM__
  ELEMENT *elem;
  VECTOR *vec,*v0,*v1;
  INT i;

  if (GBNV_list!=NULL)
    /* last process not finished properly by call of 'GetBoundaryNeighbourVectors' */
    REP_ERR_RETURN(1);

  /* count boundary node vectors */
  GBNV_n = 0;
  for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
    if (VOTYPE(vec)==NODEVEC)
      if (OBJT(MYVERTEX((NODE*)VOBJECT(vec)))==BVOBJ)
        GBNV_n++;

  /* allocate list storage: 3 pointers each */
  GBNV_mg = MYMG(theGrid);
  MarkTmpMem(MGHEAP(GBNV_mg),&GBNV_MarkKey);
  GBNV_list = (VECTOR **) GetTmpMem(MGHEAP(GBNV_mg),3*GBNV_n*sizeof(VECTOR *),GBNV_MarkKey);
  if (GBNV_list==NULL)
    REP_ERR_RETURN(1);

  /* store offset in vector index field */
  i = 0;
  for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
    if (VOTYPE(vec)==NODEVEC)
      if (OBJT(MYVERTEX((NODE*)VOBJECT(vec)))==BVOBJ)
      {
        VINDEX(vec) = i;
        GBNV_list[i] = vec;
        i += 3;
      }

  /* loop elements and fill neighbours */
  /* TODO: parallel also orphan(?) elements to be complete */
  for (elem=FIRSTELEMENT(theGrid); elem!=NULL; elem=SUCCE(elem))
    if (OBJT(elem)==BEOBJ)
      for (i=0; i<SIDES_OF_ELEM(elem); i++)
        if (ELEM_BNDS(elem,i)!=NULL)
        {
          /* 2D: two corners */
          v0 = NVECTOR(CORNER(elem,CORNER_OF_SIDE(elem,i,0)));
          v1 = NVECTOR(CORNER(elem,CORNER_OF_SIDE(elem,i,1)));

          ASSERT(GBNV_list[VINDEX(v0)]==v0);
          ASSERT(GBNV_list[VINDEX(v1)]==v1);
          GBNV_list[VINDEX(v0)+2] = v1;
          GBNV_list[VINDEX(v1)+1] = v0;
        }
  GBNV_curr = 0;

  /* this is simple in 2D: center, pred and succ in positive sense */
  *MaxListLen = 3;

  return (0);

#endif /* __TWODIM__ */

#ifdef __THREEDIM__
  /* 3D requires somewhat more work! */

  /* but it should be possible to create an oriented list
     of neighbours for each boundary vector */

  REP_ERR_RETURN (1);
#endif /* __THREEDIM__ */
}

/****************************************************************************/
/*D
   ResetGetBoundaryNeighbourVectors - reset current neighbourhood to begin of list

   SYNOPSIS:
   INT ResetGetBoundaryNeighbourVectors (void)

   PARAMETERS:
   .  void - none

   DESCRIPTION:
   This function resets the pointer to the current neighbourhood to the beginning of
   the list. 'GetBoundaryNeighbourVectors' will start again with the first one.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 else

   SEE ALSO:
   'GetBoundaryNeighbourVectors', 'PrepareGetBoundaryNeighbourVectors', 'FinishBoundaryNeighbourVectors'
   D*/
/****************************************************************************/

INT ResetGetBoundaryNeighbourVectors (void)
{
  if (GBNV_list==NULL)
    REP_ERR_RETURN(1);

  GBNV_curr = 0;
  return (0);
}

/****************************************************************************/
/*D
   GetBoundaryNeighbourVectors - get a neighbourhood of boundary vectors

   SYNOPSIS:
   INT GetBoundaryNeighbourVectors (INT dt, INT obj, INT *cnt, VECTOR *VecList[], INT *end)

   PARAMETERS:
   .  dt - datatypes of center vectors (bitwise)
   .  obj - object types of center vectors (bitwise)
   .  cnt - vector list length
   .  VecList - vector list
   .  end - if YES the currently returned list was the last one

   DESCRIPTION:
   This function returns a neighbourhood of boundary vectors, center first and the
   remainder oriented in positiv sense. The boundary vector heighbourhood lists
   are created by a call of 'PrepareGetBoundaryNeighbourVectors'. Use
   'ResetGetBoundaryNeighbourVectors' to begin again with the first
   neighbourhood and finish processing the boundary vectors with a call of
   'FinishBoundaryNeighbourVectors' which releases the temporary storage
   occupied by 'PrepareGetBoundaryNeighbourVectors'..

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 else

   SEE ALSO:
   'PrepareGetBoundaryNeighbourVectors', 'ResetGetBoundaryNeighbourVectors', 'FinishBoundaryNeighbourVectors'
   D*/
/****************************************************************************/

INT GetBoundaryNeighbourVectors (INT dt, INT obj, INT *cnt, VECTOR *VecList[])
{
  VECTOR *vec;

  *cnt = 0;

  if (GBNV_list==NULL)
    REP_ERR_RETURN(1);

  /* find next center vec matching data type */
  for (vec=GBNV_list[GBNV_curr]; GBNV_curr<3*GBNV_n; GBNV_curr+=3, vec=GBNV_list[GBNV_curr])
    if (BITWISE_TYPE(VTYPE(vec)) & dt)
      break;
  if (GBNV_curr>=3*GBNV_n)
    /* no (more) vector found */
    return (0);

  if (VOTYPE(vec)!=NODEVEC)
    REP_ERR_RETURN(1);

  /* vector, pre and succ in positive sense */
  VecList[(*cnt)++] = GBNV_list[GBNV_curr];
  VecList[(*cnt)++] = GBNV_list[GBNV_curr+1];
  VecList[(*cnt)++] = GBNV_list[GBNV_curr+2];

  /* move on to next position */
  GBNV_curr += 3;

  return (0);
}

/****************************************************************************/
/*D
   FinishBoundaryNeighbourVectors - finish processing of boundary vectors

   SYNOPSIS:
   INT FinishBoundaryNeighbourVectors (void)

   PARAMETERS:
   .  void - nonw

   DESCRIPTION:
   This function releases the temporary memory allocated by 'PrepareGetBoundaryNeighbourVectors'
   from the multigrid heap.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 else

   SEE ALSO:
   'PrepareGetBoundaryNeighbourVectors', 'ResetGetBoundaryNeighbourVectors', 'GetBoundaryNeighbourVectors'
   D*/
/****************************************************************************/

INT FinishBoundaryNeighbourVectors (void)
{
        #ifdef __TWODIM
  if (ReleaseTmpMem(MGHEAP(GBNV_mg)),GBNV_MarkKey)
    REP_ERR_RETURN(1);
        #endif

  GBNV_list = NULL;

  return (0);
}

/****************************************************************************/
/*D
   GetAllVectorsOfElement - get vector list

   SYNOPSIS:
   INT GetAllVectorsOfElement (GRID *theGrid, ELEMENT *theElement, VECTOR **vec);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  theElement - pointer to an element
   .  vec - vector list

   DESCRIPTION:
   This function gets a list of vectors corresponding to an element.

   RETURN VALUE:
   INT
   .n    number of components
   .n    -1 if error occured
   D*/
/****************************************************************************/

INT GetAllVectorsOfElement (GRID *theGrid, ELEMENT *theElement, VECTOR **vec)
{
  INT i;
  INT cnt;

  cnt = 0;
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
  {
    if (GetVectorsOfNodes(theElement,&i,vec) == GM_ERROR)
      RETURN(-1);
    cnt += i;
  }
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
  {
    if (GetVectorsOfEdges(theElement,&i,vec+cnt) == GM_ERROR)
      RETURN(-1);
    cnt += i;
  }
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
  {
    if (GetVectorsOfElement(theElement,&i,vec+cnt) == GM_ERROR)
      RETURN(-1);
    cnt += i;
  }
    #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
  {
    if (GetVectorsOfSides(theElement,&i,vec+cnt) == GM_ERROR)
      RETURN(-1);
    cnt += i;
  }
    #endif

  return (cnt);
}

/****************************************************************************/
/*D
   DisposeExtraConnections - Remove all extra connections from the grid

   SYNOPSIS:
   INT DisposeExtraConnections (GRID *theGrid);

   PARAMETERS:
   .  theGrid - grid to remove from

   DESCRIPTION:
   This function removes all extra connections from the grid, i.e. those
   that have been allocated with the 'CreateExtraConnection' function.

   RETURN VALUE:
   INT
   .n       0 if ok
   .n       1 if error occured
   D*/
/****************************************************************************/

INT DisposeExtraConnections (GRID *theGrid)
{
  VECTOR *theVector;
  MATRIX *theMatrix, *nextMatrix;
  CONNECTION *theCon;

  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    theMatrix = VSTART(theVector);
    while (theMatrix!=NULL)
    {
      nextMatrix = MNEXT(theMatrix);
      theCon = MMYCON(theMatrix);
      if (CEXTRA(theCon)) DisposeConnection(theGrid,theCon);
      theMatrix = nextMatrix;
    }
  }
  return(GM_OK);
}

INT DisposeConnectionsInGrid (GRID *theGrid)
{
  VECTOR *theVector;
  MATRIX *theMatrix, *nextMatrix;
  CONNECTION *theCon;

  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    theMatrix = VSTART(theVector);
    while (theMatrix!=NULL)
    {
      nextMatrix = MNEXT(theMatrix);
      theCon = MMYCON(theMatrix);
      DisposeConnection(theGrid,theCon);
      theMatrix = nextMatrix;
    }
  }
  return(GM_OK);
}

/****************************************************************************/
/*D
   ElementElementCreateConnection - Create connections of two elements

   SYNOPSIS:
   static INT ElementElementCreateConnection (GRID *theGrid, ELEMENT *Elem0,
   ELEMENT *Elem1, INT ActDepth, INT *ConDepth, INT* MatSize);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  Elem0,Elem1 - elements to be connected
   .  ActDepth - distance of the two elements in the element neighborship graph
   .  ConDepth - Array containing the connection depth desired. This is constructed
   .  MatSize - Array containing the connection size. This is constructed from the information in the FORMAT.

   DESCRIPTION:
   This function creates connections between all the 'VECTOR's associated
   with two given elements according to the specifications in the 'ConDepth' array.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT ElementElementCreateConnection (GRID *theGrid, ELEMENT *Elem0, ELEMENT *Elem1, INT ActDepth, INT *ConDepth, INT *MatSize)
{
  INT cnt0,cnt1,i,j,mtype;
  VECTOR *vec0[MAX_SIDES_OF_ELEM+MAX_EDGES_OF_ELEM+MAX_CORNERS_OF_ELEM+1];
  VECTOR *vec1[MAX_SIDES_OF_ELEM+MAX_EDGES_OF_ELEM+MAX_CORNERS_OF_ELEM+1];

  cnt0 = GetAllVectorsOfElement(theGrid,Elem0,vec0);
  if (Elem0 == Elem1)
  {
    for (i=0; i<cnt0; i++)
      for (j=i; j<cnt0; j++)
      {
        mtype = MatrixType[VTYPE(vec0[i])][VTYPE(vec0[j])];
        if (MatSize[mtype] > 0)
          if (ActDepth <= ConDepth[mtype])
            if (CreateConnection(theGrid,vec0[i],vec0[j])==NULL)
              RETURN(GM_ERROR);
      }
    if (NELIST_DEF_IN_GRID(theGrid))
      for (i=0; i<CORNERS_OF_ELEM(Elem0); i++)
        if (CreateElementList(theGrid,CORNER(Elem0,i),Elem0))
          RETURN(GM_ERROR);

    return (0);
  }

  cnt1 = GetAllVectorsOfElement(theGrid,Elem1,vec1);
  for (i=0; i<cnt0; i++)
    for (j=0; j<cnt1; j++)
    {
      mtype = MatrixType[VTYPE(vec0[i])][VTYPE(vec1[j])];
      if (MatSize[mtype] > 0)
        if (ActDepth <= ConDepth[mtype])
          if (CreateConnection(theGrid,vec0[i],vec1[j])==NULL)
            RETURN(GM_ERROR);
    }

  return (0);
}

/****************************************************************************/
/*D
   GetElementInfoFromSideVector	- Get pointers to elements having a common side

   SYNOPSIS:
   INT GetElementInfoFromSideVector (const VECTOR *theVector,
   ELEMENT **Elements, INT *Sides);

   PARAMETERS:
   .  theVector - given vector associated with a side of an element in 3D
   .  Elements - array to be filled with two pointers of elements
   .  Sides - array to be filled with side number within respective element

   DESCRIPTION:
   Given a 'VECTOR' associated with the side of an element, this function
   retrieves pointers to at most two elements having this side in common.
   If the side is part of the exterior boundary of the domain, then
   'Elements[1]' will be the 'NULL' pointer.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

#ifdef __THREEDIM__
INT GetElementInfoFromSideVector (const VECTOR *theVector, ELEMENT **Elements, INT *Sides)
{
  INT i;
  ELEMENT *theNeighbor;

  if (VOTYPE(theVector) != SIDEVEC)
    RETURN (1);
  Elements[0] = (ELEMENT *)VOBJECT(theVector);
  Sides[0] = VECTORSIDE(theVector);

  /* find neighbor */
  Elements[1] = theNeighbor = NBELEM(Elements[0],Sides[0]);
  if (theNeighbor == NULL) return (0);

  /* search the side */
  for (i=0; i<SIDES_OF_ELEM(theNeighbor); i++)
    if (NBELEM(theNeighbor,i) == Elements[0])
      break;

  /* found ? */
  if (i<SIDES_OF_ELEM(theNeighbor))
  {
    Sides[1] = i;
    return (0);
  }
  RETURN (1);
}
#endif

/****************************************************************************/
/*
   ResetUsedFlagInNeighborhood - Reset all 'USED' flags in neighborhood of an element

   SYNOPSIS:
   static INT ResetUsedFlagInNeighborhood (ELEMENT *theElement,
   INT ActDepth, INT MaxDepth);

   PARAMETERS:
   .  theElement - given element
   .  ActDepth - recursion depth
   .  MaxDepth - end of recursion

   DESCRIPTION:
   This function calls itself recursively and resets all 'USED' flags in the
   neighborhood of depth 'MaxDepth' of 'theElement'. For the first call
   'ActDepth' should be zero.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ResetUsedFlagInNeighborhood (ELEMENT *theElement, INT ActDepth, INT MaxDepth)
{
  int i;

  /* is anything to do ? */
  if (theElement==NULL) return (0);

  /* action */
  if (ActDepth>=0) SETUSED(theElement,0);

  /* call all neighbors recursively */
  if (ActDepth<MaxDepth)
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if (ResetUsedFlagInNeighborhood(NBELEM(theElement,i),ActDepth+1,MaxDepth)) RETURN (1);

  return (0);
}

static INT ConnectWithNeighborhood (ELEMENT *theElement, GRID *theGrid, ELEMENT *centerElement, INT *ConDepth, INT *MatSize, INT ActDepth, INT MaxDepth)
{
  int i;

  /* is anything to do ? */
  if (theElement==NULL) return (0);

  /* action */
  if (ActDepth>=0)
    if (ElementElementCreateConnection(theGrid,centerElement,theElement,
                                       ActDepth,ConDepth,MatSize))
      RETURN (1);

  /* call all neighbors recursively */
  if (ActDepth<MaxDepth)
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if (ConnectWithNeighborhood(NBELEM(theElement,i),theGrid,
                                  centerElement,ConDepth,MatSize,
                                  ActDepth+1,MaxDepth)) RETURN (1);

  return (0);
}

/****************************************************************************/
/*D
   CreateConnectionsInNeighborhood - Create connection of an element

   SYNOPSIS:
   INT CreateConnectionsInNeighborhood (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  theElement - pointer to an element

   DESCRIPTION:
   This function creates connection for all 'VECTOR's of an element
   with the depth specified in the 'FORMAT' structure.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT CreateConnectionsInNeighborhood (GRID *theGrid, ELEMENT *theElement)
{
  FORMAT *theFormat;
  INT MaxDepth;
  INT *ConDepth;
  INT *MatSize;

  /* set pointers */
  theFormat = GFORMAT(theGrid);
  MaxDepth = FMT_CONN_DEPTH_MAX(theFormat);
  ConDepth = FMT_CONN_DEPTH_PTR(theFormat);
  MatSize = FMT_S_MATPTR(theFormat);

  /* reset used flags in neighborhood */
  if (ResetUsedFlagInNeighborhood(theElement,0,MaxDepth))
    RETURN (1);

  /* create connection in neighborhood */
  if (ConnectWithNeighborhood(theElement,theGrid,theElement,ConDepth,
                              MatSize,0,MaxDepth))
    RETURN (1);

  return (0);
}

/****************************************************************************/
/*
   ConnectInsertedWithNeighborhood - Create connection of an inserted element

   SYNOPSIS:
   static INT ConnectInsertedWithNeighborhood (ELEMENT *theElement, GRID *theGrid,
   INT ActDepth, INT MaxDepth);

   PARAMETERS:
   .  theElement -
   .  theGrid -
   .  ActDepth -
   .  MaxDepth -

   DESCRITION:
   This function creates connection of an inserted element ,
   i.e. an element is inserted interactively by the user.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT ConnectInsertedWithNeighborhood (ELEMENT *theElement, GRID *theGrid, INT ActDepth, INT MaxDepth)
{
  int i;

  /* is anything to do ? */
  if (theElement==NULL) return (0);

  /* action */
  if (ActDepth>=0)
    if (CreateConnectionsInNeighborhood(theGrid,theElement))
      RETURN (1);

  /* call all neighbors recursively */
  if (ActDepth<MaxDepth)
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if (ConnectInsertedWithNeighborhood(NBELEM(theElement,i),theGrid,ActDepth+1,MaxDepth)) RETURN (1);

  return (0);
}

/****************************************************************************/
/*D
   InsertedElementCreateConnection -  Create connection of an inserted element

   SYNOPSIS:
   INT InsertedElementCreateConnection (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - grid level
   .  theElement -  pointer to an element

   DESCRITION:
   This function creates connections of an inserted element ,
   i.e. when an element is inserted interactively by the user. This function
   is inefficient and only intended for that purpose.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

INT InsertedElementCreateConnection (GRID *theGrid, ELEMENT *theElement)
{
  FORMAT *theFormat;
  INT MaxDepth;

  if (!MG_COARSE_FIXED(MYMG(theGrid)))
    RETURN (1);

  /* set pointers */
  theFormat = GFORMAT(theGrid);
  MaxDepth = (INT)(floor(0.5*(double)FMT_CONN_DEPTH_MAX(theFormat)));

  /* reset used flags in neighborhood */
  if (ResetUsedFlagInNeighborhood(theElement,0,MaxDepth))
    RETURN (1);

  /* call 'CreateConnectionsInNeighborhood'
     in a neighborhood of theElement */
  if (ConnectInsertedWithNeighborhood (theElement,theGrid,0,MaxDepth))
    RETURN (1);

  return (0);
}

/****************************************************************************/
/*D
   GridCreateConnection	- Create all connections needed on a grid level
   SYNOPSIS:
   INT GridCreateConnection (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function creates all connections needed on the grid.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT GridCreateConnection (GRID *theGrid)
{
  ELEMENT *theElement;
  VECTOR *vList[20];
  INT i,cnt;

  if (!MG_COARSE_FIXED(MYMG(theGrid)))
    RETURN (1);

  /* lets see if there's something to do */
  if (theGrid == NULL)
    return (0);

    #ifdef ModelP
        #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    DDD_XferBegin();
        #endif
        #endif

  /* set EBUILDCON-flags also in elements accessing a vector with VBUILDCON true */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    {
      EDGE *ed;

      for (i=0; i<EDGES_OF_ELEM(theElement); i++) {
        ed = GetEdge(CORNER(theElement,
                            CORNER_OF_EDGE(theElement,i,0)),
                     CORNER(theElement,
                            CORNER_OF_EDGE(theElement,i,1)));
        ASSERT (ed != NULL);
        if (EDVECTOR(ed) == NULL) {
          CreateVector(theGrid,EDGEVEC,(GEOM_OBJECT *)ed,vList);
          EDVECTOR(ed) = vList[0];
                    #ifdef ModelP
                        #ifdef __THREEDIM__
          SETPRIO(EDVECTOR(ed),PRIO(ed));
                    #endif
                        #endif
        }
      }
    }

            #ifdef ModelP
    if (!EMASTER(theElement)) continue;
                #endif

    /* see if it is set */
    if (EBUILDCON(theElement)) continue;

    /* check flags in vectors */
                #ifdef __THREEDIM__
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
    {
      GetVectorsOfSides(theElement,&cnt,vList);
      for (i=0; i<cnt; i++)
        if (VBUILDCON(vList[i])) {SETEBUILDCON(theElement,1); break;}
    }
    if (EBUILDCON(theElement)) continue;
        #endif
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    {
      GetVectorsOfEdges(theElement,&cnt,vList);
      for (i=0; i<cnt; i++)
        if (VBUILDCON(vList[i])) {SETEBUILDCON(theElement,1); break;}
    }
    if (EBUILDCON(theElement)) continue;
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
    {
      GetVectorsOfNodes(theElement,&cnt,vList);
      for (i=0; i<cnt; i++)
        if (VBUILDCON(vList[i])) {SETEBUILDCON(theElement,1); break;}
    }
  }

    #ifdef ModelP
        #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    DDD_XferEnd();
        #endif
        #endif

  /* run over all elements with EBUILDCON true and build connections */
  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
    if (EBUILDCON(theElement))             /* this is the trigger ! */
      if (CreateConnectionsInNeighborhood(theGrid,theElement))
        RETURN (1);

  return(GM_OK);
}

/****************************************************************************/
/*D
   CreateAlgebra - Creates the algebra for a grid

   SYNOPSIS:
   INT CreateAlgebra (MULTIGRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function allocates VECTORs in all geometrical objects of the grid
   (where described by format) and then calls 'GridCreateConnection'.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   D*/
/****************************************************************************/

#ifdef ModelP
static int Gather_VectorVNew (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  ((INT *)data)[0] = VNEW(theVector);

  return(0);
}

static int Scatter_VectorVNew (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  SETVNEW(theVector,MAX(VNEW(theVector),((INT *)data)[0]));

  return(0);
}

static int Scatter_GhostVectorVNew (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  SETVNEW(theVector,((INT *)data)[0]);

  return(0);
}
#endif

INT SetSurfaceClasses (MULTIGRID *theMG)
{
  VECTOR *v;
  INT level,fullrefine;

  fullrefine = TOPLEVEL(theMG);
  for (level=TOPLEVEL(theMG); level>=BOTTOMLEVEL(theMG); level--)
    for (v=PFIRSTVECTOR(GRID_ON_LEVEL(theMG,level)); v!= NULL; v=SUCCVC(v)) {
      SETNEW_DEFECT(v,(VCLASS(v)>=2));
      SETFINE_GRID_DOF(v,((VCLASS(v)>=2)&&(VNCLASS(v)<=1)));
      if (FINE_GRID_DOF(v))
        fullrefine = level;
    }
        #ifdef ModelP
  fullrefine = UG_GlobalMinINT(fullrefine);
        #endif

  FULLREFINELEVEL(theMG) = fullrefine;

  return(0);
}

INT CreateAlgebra (MULTIGRID *theMG)
{
  GRID *g;
  FORMAT *fmt;
  VECTOR *vec;
  ELEMENT *elem;
  NODE *nd;
  LINK *li;
  EDGE *ed;
  INT side,i;

  if (MG_COARSE_FIXED(theMG) == FALSE) {
    for (i=0; i<=TOPLEVEL(theMG); i++) {
      g = GRID_ON_LEVEL(theMG,i);

      if (NVEC(g)>0) return(0);

      fmt = MGFORMAT(MYMG(g));

      /* loop nodes and edges */
      for (nd=PFIRSTNODE(g); nd!=NULL; nd=SUCCN(nd)) {
        /* node vector */
        if (FMT_USES_OBJ(fmt,NODEVEC))
        {
          ASSERT(NVECTOR(nd)==NULL);
          if (CreateVector (g,NODEVEC,(GEOM_OBJECT *)nd,&vec))
            REP_ERR_RETURN (GM_ERROR);
          NVECTOR(nd) = vec;
        }

        /* edge vectors */
        if (FMT_USES_OBJ(fmt,EDGEVEC))
          for (li=START(nd); li!=NULL; li=NEXT(li)) {
            ed = MYEDGE(li);
            if (li==LINK0(ed))                                     /* to avoid double access of edges*/
            {
              ASSERT(EDVECTOR(ed)==NULL);
              if (CreateVector(g,EDGEVEC,(GEOM_OBJECT *)ed,&vec))
                REP_ERR_RETURN (GM_ERROR);
              EDVECTOR(ed) = vec;
            }
          }
      }

      /* loop elements and element sides */
      for (elem=PFIRSTELEMENT(g); elem!=NULL; elem=SUCCE(elem)) {
        /* to tell GridCreateConnection to build connections */
        if (EMASTER(elem)) SETEBUILDCON(elem,1);

        /* element vector */
        if (FMT_USES_OBJ(fmt,ELEMVEC)) {
          ASSERT(EVECTOR(elem)==NULL);
          if (CreateVector (g,ELEMVEC,(GEOM_OBJECT *)elem,&vec))
            REP_ERR_RETURN (GM_ERROR);
          SET_EVECTOR(elem,vec);
        }

        /* side vectors */
        if (FMT_USES_OBJ(fmt,SIDEVEC))
          for (side=0; side<SIDES_OF_ELEM(elem); side++)
            if (SVECTOR(elem,side)==NULL) {
              if (CreateSideVector (g,side,
                                    (GEOM_OBJECT *)elem,&vec))
                REP_ERR_RETURN (GM_ERROR);
              SET_SVECTOR(elem,side,vec);
            }
      }
      MG_COARSE_FIXED(theMG) = TRUE;

      /* now connections */
      if (GridCreateConnection(g))
        REP_ERR_RETURN (1);
    }
  }

    #ifdef ModelP
    #ifndef __EXCHANGE_CONNECTIONS__
  MGCreateConnection(theMG);
    #endif
  /* update VNEW-flags */
  DDD_IFExchange(BorderVectorSymmIF,sizeof(INT),
                 Gather_VectorVNew,Scatter_VectorVNew);
  DDD_IFOneway(VectorAllIF,IF_FORWARD,sizeof(INT),
               Gather_VectorVNew,Scatter_GhostVectorVNew);
    #endif

  SetSurfaceClasses(theMG);

  return(GM_OK);
}

/****************************************************************************/
/*D
   MGCreateConnection - Create all connections in multigrid

   SYNOPSIS:
   INT MGCreateConnection (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - pointer to mulrigrid

   DESCRIPTION:
   This function creates all connections in multigrid.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT MGCreateConnection (MULTIGRID *theMG)
{
  INT i;
  GRID *theGrid;
  ELEMENT *theElement;

  if (!MG_COARSE_FIXED(theMG))
    RETURN (1);

  for (i=0; i<=theMG->topLevel; i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      SETEBUILDCON(theElement,1);
    if (GridCreateConnection(theGrid)) RETURN (1);
  }

  return (0);
}

/****************************************************************************/
/*D
   ResetAlgebraInternalFlags - Reset USED and EBUILDCON flags in elements

   SYNOPSIS:
   INT PrepareAlgebraModification (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - pointer to multigrid

   DESCRIPTION:
   This function resets  USED and EBUILDCON flags in elements.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT PrepareAlgebraModification (MULTIGRID *theMG)
{
  int j,k;
  ELEMENT *theElement;
  VECTOR *theVector;
  MATRIX *theMatrix;

  j = theMG->topLevel;
  for (k=0; k<=j; k++)
  {
    for (theElement=PFIRSTELEMENT(GRID_ON_LEVEL(theMG,k)); theElement!=NULL; theElement=SUCCE(theElement))
    {
      SETUSED(theElement,0);
      SETEBUILDCON(theElement,0);
    }
    for (theVector=PFIRSTVECTOR(GRID_ON_LEVEL(theMG,k)); theVector!= NULL; theVector=SUCCVC(theVector))
      SETVBUILDCON(theVector,0);
    for (theVector=PFIRSTVECTOR(GRID_ON_LEVEL(theMG,k)); theVector!= NULL; theVector=SUCCVC(theVector))
    {
      SETVNEW(theVector,0);
      for (theMatrix=VSTART(theVector); theMatrix!=NULL; theMatrix = MNEXT(theMatrix))
        SETMNEW(theMatrix,0);
    }
  }

  return (0);
}

/****************************************************************************/
/*D
   ElementElementCheck - Check connection of two elements

   SYNOPSIS:
   static INT ElementElementCheck (GRID *theGrid,
   ELEMENT *Elem0, ELEMENT *Elem1, INT ActDepth, INT *ConDepth, INT *MatSize);

   PARAMETERS:
   .  theGrid -  grid level to check
   .  Elem0,Elem1 - elements to be checked.
   .  ActDepth - recursion depth
   .  ConDepth - connection depth as provided by format.
   .  MatSize - matrix sizes as provided by format

   DESCRIPTION:
   This function recursively checks connection of two elements.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     != 0 if errors occured.
   D*/
/****************************************************************************/

static INT ElementElementCheck (GRID *theGrid, ELEMENT *Elem0, ELEMENT *Elem1, INT ActDepth, INT *ConDepth, INT *MatSize)
{
  INT cnt0,cnt1,i,j,mtype;
  VECTOR *vec0[MAX_SIDES_OF_ELEM+MAX_EDGES_OF_ELEM+MAX_CORNERS_OF_ELEM+1];
  VECTOR *vec1[MAX_SIDES_OF_ELEM+MAX_EDGES_OF_ELEM+MAX_CORNERS_OF_ELEM+1];
  CONNECTION *theCon;
  char msg[128];
  INT errors = 0;

  sprintf(msg,PFMT " ERROR: missing connection between elem0=" EID_FMTX " elem1=" EID_FMTX,
          me,EID_PRTX(Elem0),EID_PRTX(Elem1));

  cnt0 = GetAllVectorsOfElement(theGrid,Elem0,vec0);
  if (Elem0 == Elem1)
  {
    for (i=0; i<cnt0; i++)
      for (j=0; j<cnt0; j++)
      {
        mtype = MatrixType[VTYPE(vec0[i])][VTYPE(vec0[j])];
        if (MatSize[mtype] > 0)
          if (ActDepth <= ConDepth[mtype])
          {
            /* check connection in both directions */
            theCon = GetConnection(vec0[i],vec0[j]);
            if (theCon==NULL)
            {
              errors++;
              UserWriteF("%s vec0[%d]=" VINDEX_FMTX
                         " to vec0[%d]=" VINDEX_FMTX "\n",
                         msg,i,VINDEX_PRTX(vec0[i]),
                         j,VINDEX_PRTX(vec0[j]));
            }
            else
            {
              theCon = GetConnection(vec0[j],vec0[i]);
              if (theCon==NULL)
              {
                errors++;
                UserWriteF("%s vec0[%d]=" VINDEX_FMTX
                           " to vec0[%d]=" VINDEX_FMTX "\n",
                           msg,j,VINDEX_PRTX(vec0[j]),
                           i,VINDEX_PRTX(vec0[i]));
              }
              else
              {
                SETCUSED(theCon,1);
              }
            }
          }
      }
    return (errors);
  }

  cnt1 = GetAllVectorsOfElement(theGrid,Elem1,vec1);
  for (i=0; i<cnt0; i++)
    for (j=0; j<cnt1; j++)
    {
      mtype = MatrixType[VTYPE(vec0[i])][VTYPE(vec1[j])];
      if (MatSize[mtype] > 0)
        if (ActDepth <= ConDepth[mtype])
        {
          /* check connection in both directions */
          theCon = GetConnection(vec0[i],vec1[j]);
          if (theCon==NULL)
          {
            errors++;
            UserWriteF("%s vec0[%d]=" VINDEX_FMTX
                       " to vec1[%d]=" VINDEX_FMTX "\n",
                       msg,i,VINDEX_PRTX(vec0[i]),
                       j,VINDEX_PRTX(vec1[j]));
          }
          else
          {
            theCon = GetConnection(vec1[j],vec0[i]);
            if (theCon == NULL)
            {
              errors++;
              UserWriteF("%s vec1[%d]=" VINDEX_FMTX
                         " to vec0[%d]=%x/" VINDEX_FMTX "\n",
                         msg,j,VINDEX_PRTX(vec1[j]),
                         i,VINDEX_PRTX(vec0[i]));
            }
            else
              SETCUSED(theCon,1);
          }
        }
    }
  return (errors);
}

static INT CheckNeighborhood (GRID *theGrid, ELEMENT *theElement, ELEMENT *centerElement, INT *ConDepth, INT ActDepth, INT MaxDepth, INT *MatSize)
{
  INT i,errors = 0;

  /* is anything to do ? */
  if (theElement==NULL) return (0);

  /* action */
  if (ActDepth>=0)
    if ((errors+=ElementElementCheck(theGrid,centerElement,theElement,ActDepth,ConDepth,MatSize))!=0)
      return (errors);

  /* call all neighbors recursively */
  if (ActDepth<MaxDepth)
    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      if ((errors+=CheckNeighborhood(theGrid,NBELEM(theElement,i),centerElement,ConDepth,ActDepth+1,MaxDepth,MatSize))!=0)
        return(errors);

  return (0);
}

/****************************************************************************/
/*D
   ElementCheckConnection - Check connection of the element

   SYNOPSIS:
   INT   ElementCheckConnection (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  theElement - pointer to element

   DESCRIPTION:
   This function checks all connections of the given element.

   RETURN VALUE:
   INT
   .n   number of errors
   .n   0 if ok
   .n   != 0 if errors occured in that connection.
   D*/
/****************************************************************************/

INT ElementCheckConnection (GRID *theGrid, ELEMENT *theElement)
{
  FORMAT *theFormat;
  INT MaxDepth;
  INT *ConDepth;
  INT *MatSize;

  /* set pointers */
  theFormat = GFORMAT(theGrid);
  MaxDepth = FMT_CONN_DEPTH_MAX(theFormat);
  ConDepth = FMT_CONN_DEPTH_PTR(theFormat);
  MatSize = FMT_S_MATPTR(theFormat);

  /* call elements recursivly */
  return(CheckNeighborhood(theGrid,theElement,theElement,ConDepth,0,MaxDepth,MatSize));
}

/****************************************************************************/
/*D
   CheckConnections - Check if connections are correctly allocated

   SYNOPSIS:
   INT CheckConnections (GRID *theGrid);

   PARAMETERS:
   .  theGrid -  grid level to check

   DESCRIPTION:
   This function checks if connections are correctly allocated.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

static INT CheckConnections (GRID *theGrid)
{
  ELEMENT *theElement;
  INT errors;
  INT error;

  errors = 0;

  for (theElement=FIRSTELEMENT(theGrid);
       theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    if ((error=ElementCheckConnection(theGrid,theElement))!=0)
    {
      UserWriteF(PFMT "%d: element=" EID_FMTX " has bad connections\n",
                 me, EID_PRTX(theElement));
      errors+=error;
    }
  }

  return(errors);
}


/****************************************************************************/
/*D
   CheckVector - Checks validity of geom_object	and its vector

   SYNOPSIS:
   INT CheckVector (const FORMAT *fmt, const INT s2p[], GEOM_OBJECT *theObject, char *ObjectString,
                VECTOR *theVector, INT VectorObjType, INT side);

   PARAMETERS:
   .  fmt - FORMAT of associated multigrid
   .  s2p - subdomain to part table
   .  theObject - the object which points to theVector
   .  ObjectString - for error message
   .  theVector - the vector of theObject
   .  VectorObjType - NODEVEC,...
   .  side - element side for SIDEVEC, NOSIDE else

   DESCRIPTION:
   This function checks the consistency between an geom_object and
   its vector.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

static INT CheckVector (const FORMAT *fmt, const INT s2p[], GEOM_OBJECT *theObject, char *ObjectString,
                        VECTOR *theVector, INT VectorObjType, INT side)
{
  GEOM_OBJECT *VecObject;
  MATRIX *mat;
  INT errors = 0,vtype,DomPart,ds;

  if (theVector == NULL)
  {
    /* check if size is really 0 */
    DomPart = GetDomainPart(s2p,theObject,side);
    vtype = FMT_PO2T(fmt,DomPart,VectorObjType);
    ds = FMT_S_VEC_TP(fmt,vtype);
    if (ds>0)
    {
      errors++;
      UserWriteF("%d: %s ID=%ld  has NO VECTOR",me, ObjectString,
                 ID(theObject));
                        #ifdef ModelP
                        #ifdef __THREEDIM__
      if (VectorObjType == EDGEVEC)
        UserWriteF(" prio=%d",PRIO((EDGE*)theObject));
                        #endif
                        #endif
      UserWrite("\n");
    }
  }
  else
  {
    ds = FMT_S_VEC_TP(fmt,VTYPE(theVector));
    if (ds==0)
    {
      errors++;
      UserWriteF("%d: %s ID=%ld  exists but should not\n",me, ObjectString,
                 ID(theObject));
    }
    SETVCUSED(theVector,1);

    VecObject = VOBJECT(theVector);
    if (VecObject == NULL)
    {
      errors++;
      UserWriteF("%d: vector=" VINDEX_FMTX " %s GID=" GIDFMT " has NO BACKPTR\n",
                 me, VINDEX_PRTX(theVector), ObjectString,
                 (OBJT(theObject)==BEOBJ || OBJT(theObject)==IEOBJ) ?
                 EGID(&(theObject->el)) : (OBJT(theObject)==NDOBJ) ?
                 GID(&(theObject->nd)) :
                                #ifdef __TWODIM__
                 ((int *)(theObject))[0]
                                #else
                 GID(&(theObject->ed))
                                #endif
                 );
    }
    else
    {
      if (VOTYPE(theVector) != VectorObjType)
      {
        errors++;
        UserWriteF("%d: %s vector=" VINDEX_FMTX " has incompatible type=%d, "
                   "should be type=%s\n",
                   me, ObjectString, VINDEX_PRTX(theVector), VTYPE(theVector),
                   ObjTypeName[VectorObjType]);
      }

      if (VecObject != theObject)
      {
        if (OBJT(VecObject) != OBJT(theObject))
        {
          int error = 1;

          /* both objects may be elements */
          if ((OBJT(VecObject)==BEOBJ || OBJT(VecObject)==IEOBJ) &&
              (OBJT(theObject)==BEOBJ || OBJT(theObject)==IEOBJ) )
          {
            ELEMENT *theElement = (ELEMENT *)theObject;
            ELEMENT *vecElement = (ELEMENT *)VecObject;
            int i;

                                                #ifdef ModelP
            if (EMASTER(theElement) ||
                EMASTER(vecElement) )
            {
                                                #endif
            for (i=0; i<SIDES_OF_ELEM(theElement); i++)
              if (NBELEM(theElement,i) == vecElement)
              {
                /* they are neighbors -> ok */
                error = 0;
                break;
              }
                                                #ifdef ModelP
          }
                                                #endif
            if (error)
            {
              UserWriteF("%d: vector=" VINDEX_FMTX " has type %s, but points "
                         "to wrong vecobj=" EID_FMTX " NO NB of obj=" EID_FMTX "\n",
                         me,VINDEX_PRTX(theVector),ObjectString,
                         EID_PRTX(vecElement),EID_PRTX(theElement));
            }
          }
          else
          {
            errors++;
            UserWriteF("%d: vector=" VINDEX_FMTX " has type %s, but points "
                       "to wrong obj=%d type OBJT=%d\n",
                       me,VINDEX_PRTX(theVector),ObjectString,ID(VecObject),
                       OBJT(VecObject));
          }
        }
        else
        {
                                        #ifdef __THREEDIM__
          if (VectorObjType == SIDEVEC)
          {
            /* TODO: check side vector */
          }
          else
                                        #endif
          {
            errors++;
            UserWriteF("%d: %s vector=" VINDEX_FMTX " is referenced by "
                       "obj0=%x, but points to wrong obj1=%x\n",
                       me, ObjectString, VINDEX_PRTX(theVector),
                       theObject, VecObject);
                                                #ifdef ModelP
            if (strcmp(ObjectString,"EDGE")==0)
              UserWriteF("%d: obj0: n0=%d n1=%d  obj1: "
                         "n0=%d n1=%d\n",me,
                         GID(NBNODE(LINK0(&(theObject->ed)))),
                         GID(NBNODE(LINK1(&(theObject->ed)))),
                         GID(NBNODE(LINK0(&(VecObject->ed)))),
                         GID(NBNODE(LINK1(&(VecObject->ed)))) );
                                                #endif
          }
        }
      }
    }
    /* check connectivity of matrices */
    for (mat=VSTART(theVector); mat!=NULL; mat=MNEXT(mat))
      if (MDEST(mat)==NULL)
      {
        errors++;
        UserWriteF("%d: %s vector=" VINDEX_FMTX ": matrix dest==NULL\n",
                   me, ObjectString, VINDEX_PRTX(theVector));
      }
      else if (MDEST(MADJ(mat))!=theVector)
      {
        errors++;
        UserWriteF("%d: %s vector=" VINDEX_FMTX ": adj matrix dest does not coincide"
                   " with vector conn=%x mat=%x mdest=%x\n",
                   me, ObjectString, VINDEX_PRTX(theVector),MMYCON(mat),MDEST(mat),MDEST(MADJ(mat)));
      }

  }

  return(errors);
}

/****************************************************************************/
/*D
   CheckAlgebra - Check the algebraic part of the data structure

   SYNOPSIS:
   INT CheckAlgebra (GRID *theGrid);

   PARAMETERS:
   .  theGrid -  grid level to check

   DESCRIPTION:
   This function checks the consistency of the algebraic data structures
   including the interconnections between the geometric part.
   This function assumes a correct geometric data structure.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR	if error occured.
   D*/
/****************************************************************************/

INT CheckAlgebra (GRID *theGrid)
{
  FORMAT *fmt;
  ELEMENT *theElement;
  NODE *theNode;
  VECTOR *theVector;
  EDGE *theEdge;
  LINK *theLink;
  INT errors;
  INT *s2p;
#       ifdef __THREEDIM__
  INT i;
#       endif

  errors = 0;

  if ((GLEVEL(theGrid)==0) && !MG_COARSE_FIXED(MYMG(theGrid)))
  {
    if ((NVEC(theGrid)>0) || (NC(theGrid)>0))
    {
      errors++;
      UserWriteF("coarse grid not fixed but vectors allocated\n");
    }
    return(errors);
  }

  s2p = BVPD_S2P_PTR(MG_BVPD(MYMG(theGrid)));
  fmt = MGFORMAT(MYMG(theGrid));

  /* reset USED flag */
  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
  {
    SETVCUSED(theVector,0);
  }

  /* check pointers to element, side, edge vector */
  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {

    /* check element vectors */
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
    {
      theVector = EVECTOR(theElement);
      errors += CheckVector(fmt,s2p,(GEOM_OBJECT *) theElement, "ELEMENT",
                            theVector, ELEMVEC,NOSIDE);
    }

                #ifdef __THREEDIM__
    /* check side vectors */
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
    {
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        theVector = SVECTOR(theElement,i);
        errors += CheckVector(fmt,s2p,(GEOM_OBJECT *) theElement, "ELEMSIDE",
                              theVector, SIDEVEC,i);
      }
    }
                #endif
  }

  for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
  {

    /* check node vectors */
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
    {
      theVector = NVECTOR(theNode);
      errors += CheckVector(fmt,s2p,(GEOM_OBJECT *) theNode, "NODE", theVector,
                            NODEVEC,NOSIDE);
    }

    /* check edge vectors */
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
    {
      for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
      {
        theEdge = GetEdge(theNode,NBNODE(theLink));
        if (theEdge != NULL) {
          theVector = EDVECTOR(theEdge);
          errors += CheckVector(fmt,s2p,(GEOM_OBJECT *) theEdge, "EDGE",
                                theVector, EDGEVEC,NOSIDE);
        }
      }
    }
  }

  /* check USED flag */
  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
  {
    if (VCUSED(theVector) != 1)
    {
      errors++;
      UserWriteF("%d: vector ID=%d NOT referenced by an geom_object: "
                 "vtype=%d, objptr=%x",
                 me, ID(theVector), VTYPE(theVector), VOBJECT(theVector));
      if (VOBJECT(theVector) != NULL)
        UserWriteF(" objtype=%d\n",OBJT(VOBJECT(theVector)));
      else
        UserWrite("\n");
    }
    else
      SETVCUSED(theVector,0);
  }

  /* check validity of all defined connections */
  errors += CheckConnections(theGrid);

  /* reset flags in connections */
  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
  {
    MATRIX  *theMatrix;

    for (theMatrix=VSTART(theVector); theMatrix!=NULL;
         theMatrix = MNEXT(theMatrix)) SETCUSED(MMYCON(theMatrix),0);
  }

  /* set flags in connections */
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
  {
    MATRIX  *theMatrix;

    for (theMatrix=VSTART(theVector); theMatrix!=NULL;
         theMatrix = MNEXT(theMatrix)) SETMUSED(MADJ(theMatrix),1);
  }

  /* check matrices and connections */
  for (theVector=PFIRSTVECTOR(theGrid);
       theVector!=NULL;
       theVector=SUCCVC(theVector))
  {
    MATRIX  *theMatrix;
                #ifdef ModelP
    INT prio = PRIO(theVector);
                #endif

    for (theMatrix=VSTART(theVector);
         theMatrix!=NULL;
         theMatrix = MNEXT(theMatrix))
    {
      MATRIX *Adj = MADJ(theMatrix);

      if (MDEST(theMatrix) == NULL)
      {
        errors++;
        UserWriteF(PFMT "ERROR: matrix %x has no dest, start vec="
                   VINDEX_FMTX "\n",
                   me,theMatrix,VINDEX_PRTX(theVector));
      }
      if (MDEST(Adj) != theVector)
      {
        errors++;
        UserWriteF(PFMT "ERROR: dest=%x of adj matrix "
                   " unequal vec=" VINDEX_FMTX "\n",
                   me, MDEST(Adj),VINDEX_PRTX(theVector));
      }

                        #ifdef ModelP
      if (prio != PrioHGhost)
                        #endif
      if (MUSED(theMatrix)!=1 &&  !CEXTRA(MMYCON(theMatrix)))
      {
        errors++;
        UserWriteF(PFMT "ERROR: connection dead vec=" VINDEX_FMTX
                   " vector=" VINDEX_FMTX " con=%x mat=%x matadj=%x level(vec)=%d is_extra_con %d\n",
                   me,VINDEX_PRTX(theVector),VINDEX_PRTX(MDEST(theMatrix)),
                   MMYCON(theMatrix),MDEST(theMatrix),MDEST(MADJ(theMatrix)),
                   GLEVEL(theGrid),CEXTRA(MMYCON(theMatrix)));
      }

                        #ifdef ModelP
      if (GHOSTPRIO(prio) && !CEXTRA(MMYCON(theMatrix)))
      {
        errors++;
        UserWriteF(PFMT "ERROR: ghost vector has matrix vec="
                   VINDEX_FMTX " con=%x mat=%x\n",
                   me,VINDEX_PRTX(theVector),MMYCON(theMatrix),theMatrix);
      }
                        #endif
    }
                #ifdef __INTERPOLATION_MATRIX__
    for (theMatrix=VISTART(theVector);
         theMatrix!=NULL;
         theMatrix = MNEXT(theMatrix))
      if (MDEST(theMatrix) == NULL) {
        errors++;
        UserWriteF(PFMT "ERROR: interpolation matrix %x has no dest,"
                   " start vec=" VINDEX_FMTX "\n",me,theMatrix,
                   VINDEX_PRTX(theVector));
      }
        #endif
  }

  return(errors);
}

/****************************************************************************/
/*D
   VectorInElement -  Decide whether a vector corresponds to an element or not
   SYNOPSIS:
   INT VectorInElement (ELEMENT *theElement, VECTOR *theVector);

   PARAMETERS:
   .  theElement - pointer to element
   .  theVector - pointer to a vector

   DESCRIPTION:
   This function decides whether a given vector belongs to the given element, or
   one of its sides, edges or nodes.

   RETURN VALUE:
   INT
   .n    0 if does not correspond
   .n    1 if does correspond.
   D*/
/****************************************************************************/

INT VectorInElement (ELEMENT *theElement, VECTOR *theVector)
{
  INT i;
  VECTOR *vList[20];
  INT cnt;

  if (VOTYPE(theVector) == ELEMVEC)
  {
    GetVectorsOfElement(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
      if (vList[i]==theVector) RETURN(1);
  }
        #ifdef __THREEDIM__
  if (VOTYPE(theVector) == SIDEVEC)
  {
    GetVectorsOfSides(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
      if (vList[i]==theVector) RETURN(1);
  }
    #endif
  if (VOTYPE(theVector) == EDGEVEC)
  {
    GetVectorsOfEdges(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
      if (vList[i]==theVector) RETURN(1);
  }
  if (VOTYPE(theVector) == NODEVEC)
  {
    GetVectorsOfNodes(theElement,&cnt,vList);
    for (i=0; i<cnt; i++)
      if (vList[i]==theVector) RETURN(1);
  }

  return (0);
}

/****************************************************************************/
/*D
   VectorPosition - Calc coordinate position of vector

   SYNOPSIS:
   INT VectorPosition (const VECTOR *theVector, DOUBLE *position);

   PARAMETERS:
   .  theVector - a given vector
   .  position - array to be filled

   DESCRIPTION:
   This function calcs physical position of vector. For edge vectors the
   midpoint is returned, and for sides and elements the center of mass.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT VectorPosition (const VECTOR *theVector, DOUBLE *position)
{
  INT i;
  EDGE *theEdge;
        #ifdef __THREEDIM__
  ELEMENT *theElement;
  INT theSide,j;
        #endif

  ASSERT(theVector != NULL);

  switch (VOTYPE(theVector))
  {
  case NODEVEC :
    for (i=0; i<DIM; i++)
      position[i] = CVECT(MYVERTEX((NODE*)VOBJECT(theVector)))[i];
    return (0);

  case EDGEVEC :
    theEdge = (EDGE*)VOBJECT(theVector);
    for (i=0; i<DIM; i++)
      position[i] = 0.5*(CVECT(MYVERTEX(NBNODE(LINK0(theEdge))))[i] +
                         CVECT(MYVERTEX(NBNODE(LINK1(theEdge))))[i]   );
    return (0);
                #ifdef __THREEDIM__
  case SIDEVEC :
    theElement = (ELEMENT *)VOBJECT(theVector);
    theSide = VECTORSIDE(theVector);
    for (i=0; i<DIM; i++)
    {
      position[i] = 0.0;
      for(j=0; j<CORNERS_OF_SIDE(theElement,theSide); j++)
        position[i] += CVECT(MYVERTEX(CORNER(theElement,CORNER_OF_SIDE(theElement,theSide,j))))[i];
      position[i] /= CORNERS_OF_SIDE(theElement,theSide);
    }
    return (0);
                #endif
  case ELEMVEC :
    /* calculate center of mass */
    CalculateCenterOfMass( (ELEMENT *) VOBJECT(theVector), position );
    return (0);

  default : PrintErrorMessage('E',"VectorPosition","unrecognized object type for vector");
    assert(0);
  }

  RETURN (GM_ERROR);
}


/****************************************************************************/
/*D
   SeedVectorClasses -  Initialize vector classes

   SYNOPSIS:
   INT SeedVectorClasses (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - given grid
   .  theElement - given element

   DESCRIPTION:
   Initialize vector class in all vectors associated with given element with 3.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT SeedVectorClasses (GRID *theGrid, ELEMENT *theElement)
{
  INT i;
  VECTOR *vList[20];
  INT cnt;

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
  {
    GetVectorsOfElement(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVCLASS(vList[i],3);
  }
        #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
  {
    GetVectorsOfSides(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVCLASS(vList[i],3);
  }
    #endif
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
  {
    GetVectorsOfEdges(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVCLASS(vList[i],3);
  }
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
  {
    GetVectorsOfNodes(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVCLASS(vList[i],3);
  }
  return (0);
}

/****************************************************************************/
/*D
   ClearVectorClasses - Reset vector classes

   SYNOPSIS:
   INT ClearVectorClasses (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   Reset all vector classes in all vectors of given grid to 0.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT ClearVectorClasses (GRID *theGrid)
{
  VECTOR *theVector;

  /* reset class of each vector to 0 */
  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
    SETVCLASS(theVector,0);

  return(0);
}
/****************************************************************************/
/*D
   PropagateVectorClasses - Compute vector classes after initialization

   SYNOPSIS:
   INT PropagateVectorClasses (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   After vector classes have been reset and initialized, this function
   now computes the class 2 and class 1 vectors.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error occured
   D*/
/****************************************************************************/

#ifdef ModelP
static int Gather_VectorVClass (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  PRINTDEBUG(gm,1,(PFMT "Gather_VectorVClass(): v=" VINDEX_FMTX " vclass=%d\n",
                   me,VINDEX_PRTX(theVector),VCLASS(theVector)))

    ((INT *)data)[0] = VCLASS(theVector);

  return(0);
}

static int Scatter_VectorVClass (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  SETVCLASS(theVector,MAX(VCLASS(theVector),((INT *)data)[0]));

  PRINTDEBUG(gm,2,(PFMT "Scatter_VectorVClass(): v=" VINDEX_FMTX " vclass=%d\n",
                   me,VINDEX_PRTX(theVector),VCLASS(theVector)))

  return(0);
}

static int Scatter_GhostVectorVClass (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  SETVCLASS(theVector,((INT *)data)[0]);

  PRINTDEBUG(gm,1,(PFMT "Scatter_GhostVectorVClass(): v=" VINDEX_FMTX " vclass=%d\n",
                   me,VINDEX_PRTX(theVector),VCLASS(theVector)))

  return(0);
}
#endif

INT PropagateVectorClasses (GRID *theGrid)
{
  VECTOR *theVector;
  MATRIX *theMatrix;

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateVectorClasses():"
                   " 1. communication on level %d\n",me,GLEVEL(theGrid)))
  /* exchange VCLASS of vectors */
  DDD_IFAExchange(BorderVectorSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_VectorVClass, Scatter_VectorVClass);
    #endif

  /* set vector classes in the algebraic neighborhood to 2 */
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
    if ((VCLASS(theVector)==3)&&(VSTART(theVector)!=NULL))
      for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL;
           theMatrix=MNEXT(theMatrix))
        if ((VCLASS(MDEST(theMatrix))<3)&&
            (CEXTRA(MMYCON(theMatrix))!=1))
          SETVCLASS(MDEST(theMatrix),2);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateVectorClasses(): 2. communication\n",
                   me))
  /* exchange VCLASS of vectors */
  DDD_IFAExchange(BorderVectorSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_VectorVClass, Scatter_VectorVClass);
    #endif

  /* set vector classes in the algebraic neighborhood to 1 */
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
    if ((VCLASS(theVector)==2)&&(VSTART(theVector)!=NULL))
      for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL;
           theMatrix=MNEXT(theMatrix))
        if ((VCLASS(MDEST(theMatrix))<2)&&
            (CEXTRA(MMYCON(theMatrix))!=1))
          SETVCLASS(MDEST(theMatrix),1);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateVectorClasses(): 3. communication\n",
                   me))
  /* exchange VCLASS of vectors */
  DDD_IFAExchange(BorderVectorSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_VectorVClass, Scatter_VectorVClass);
  /* send VCLASS to ghosts */
  DDD_IFAOneway(VectorAllIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(INT),
                Gather_VectorVClass, Scatter_GhostVectorVClass);
    #endif

  return(0);
}

/****************************************************************************/
/*D
   ClearNextVectorClasses - Reset class of the vectors on the next level

   SYNOPSIS:
   INT ClearNextVectorClasses (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function clears VNCLASS flag in all vectors. This is the first step to
   compute the class of the dofs on the *NEXT* level, which
   is also the basis for determining copies.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT ClearNextVectorClasses (GRID *theGrid)
{
  VECTOR *theVector;

  /* reset class of each vector to 0 */
  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
    SETVNCLASS(theVector,0);

  /* now the refinement algorithm will initialize the class 3 vectors */
  /* on the *NEXT* level.                                                                                       */
  return(0);
}

/****************************************************************************/
/*D
   SeedNextVectorClasses - Set 'VNCLASS' in all vectors associated with element

   SYNOPSIS:
   INT SeedNextVectorClasses (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - given grid
   .  theElement - pointer to element

   DESCRIPTION:
   Set 'VNCLASS' in all vectors associated with the element to 3.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT SeedNextVectorClasses (GRID *theGrid, ELEMENT *theElement)
{
  INT i;
  VECTOR *vList[20];
  INT cnt;

  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
  {
    GetVectorsOfElement(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVNCLASS(vList[i],3);
  }
        #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
  {
    GetVectorsOfSides(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVNCLASS(vList[i],3);
  }
        #endif
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
  {
    GetVectorsOfEdges(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVNCLASS(vList[i],3);
  }
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
  {
    GetVectorsOfNodes(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) SETVNCLASS(vList[i],3);
  }
  return (0);
}


/****************************************************************************/
/*D
   PropagateNextVectorClasses - Compute 'VNCLASS' in all vectors of a grid level

   SYNOPSIS:
   INT PropagateNextVectorClasses (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   Computes values of 'VNCLASS' field in all vectors after seed.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   D*/
/****************************************************************************/

#ifdef ModelP
static int Gather_VectorVNClass (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  PRINTDEBUG(gm,2,(PFMT "Gather_VectorVNClass(): v=" VINDEX_FMTX " vnclass=%d\n",
                   me,VINDEX_PRTX(theVector),VNCLASS(theVector)))

    ((INT *)data)[0] = VNCLASS(theVector);

  return(GM_OK);
}

static int Scatter_VectorVNClass (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  SETVNCLASS(theVector,MAX(VNCLASS(theVector),((INT *)data)[0]));

  PRINTDEBUG(gm,2,(PFMT "Scatter_VectorVNClass(): v=" VINDEX_FMTX " vnclass=%d\n",
                   me,VINDEX_PRTX(theVector),VNCLASS(theVector)))

  return(GM_OK);
}

static int Scatter_GhostVectorVNClass (DDD_OBJ obj, void *data)
{
  VECTOR *theVector = (VECTOR *)obj;

  SETVNCLASS(theVector,((INT *)data)[0]);

  PRINTDEBUG(gm,2,(PFMT "Scatter_GhostVectorVNClass(): v=" VINDEX_FMTX " vnclass=%d\n",
                   me,VINDEX_PRTX(theVector),VNCLASS(theVector)))

  return(GM_OK);
}
#endif

INT PropagateNextVectorClasses (GRID *theGrid)
{
  VECTOR *theVector;
  MATRIX *theMatrix;

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextVectorClasses(): 1. communication\n",me))
  /* exchange VNCLASS of vectors */
  DDD_IFAExchange(BorderVectorSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_VectorVNClass, Scatter_VectorVNClass);
    #endif

  /* set vector classes in the algebraic neighborhood to 2 */
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
    if ((VNCLASS(theVector)==3)&&(VSTART(theVector)!=NULL))
      for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL;
           theMatrix=MNEXT(theMatrix))
        if ((VNCLASS(MDEST(theMatrix))<3)
            &&(CEXTRA(MMYCON(theMatrix))!=1))
          SETVNCLASS(MDEST(theMatrix),2);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextVectorClasses(): 2. communication\n",me))
  /* exchange VNCLASS of vectors */
  DDD_IFAExchange(BorderVectorSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_VectorVNClass, Scatter_VectorVNClass);
    #endif

  /* set vector classes in the algebraic neighborhood to 1 */
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector))
    if ((VNCLASS(theVector)==2)&&(VSTART(theVector)!=NULL))
      for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL;
           theMatrix=MNEXT(theMatrix))
        if ((VNCLASS(MDEST(theMatrix))<2)
            &&(CEXTRA(MMYCON(theMatrix)) !=1))
          SETVNCLASS(MDEST(theMatrix),1);

    #ifdef ModelP
  PRINTDEBUG(gm,1,("\n" PFMT "PropagateNextVectorClasses(): 3. communication\n",me))
  /* exchange VNCLASS of vectors */
  DDD_IFAExchange(BorderVectorSymmIF,GRID_ATTR(theGrid),sizeof(INT),
                  Gather_VectorVNClass, Scatter_VectorVNClass);
  /* send VCLASS to ghosts */
  DDD_IFAOneway(VectorAllIF,GRID_ATTR(theGrid),IF_FORWARD,sizeof(INT),
                Gather_VectorVNClass, Scatter_GhostVectorVNClass);
    #endif

  return(0);
}

/****************************************************************************/
/*D
   MaxNextVectorClass - Returns highest vector class of a dof on next level

   SYNOPSIS:
   INT MaxNextVectorClass (GRID *theGrid, ELEMENT *theElement);

   PARAMETERS:
   .  theGrid - pointer to a grid
   .  theElement - pointer to a element

   DECRIPTION:
   This function returns highest 'VNCLASS' of a vector associated with the
   element.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT MaxNextVectorClass (GRID *theGrid, ELEMENT *theElement)
{
  INT i,m;
  VECTOR *vList[20];
  INT cnt;

  m = 0;
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
  {
    GetVectorsOfElement(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) m = MAX(m,VNCLASS(vList[i]));
  }
        #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
  {
    GetVectorsOfSides(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) m = MAX(m,VNCLASS(vList[i]));
  }
        #endif
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC))
  {
    GetVectorsOfEdges(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) m = MAX(m,VNCLASS(vList[i]));
  }
  if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,NODEVEC))
  {
    GetVectorsOfNodes(theElement,&cnt,vList);
    for (i=0; i<cnt; i++) m = MAX(m,VNCLASS(vList[i]));
  }
  return (m);
}

/****************************************************************************/
/*                                                                          */
/* Function:  LexCompare													*/
/*                                                                          */
/****************************************************************************/

static INT LexCompare (VECTOR **pvec1, VECTOR **pvec2)
{
  DOUBLE_VECTOR pv1,pv2;
  DOUBLE diff[DIM];

  PRINTDEBUG(gm,1,("%d: LexCompare %4d %4d",
                   me,VINDEX(*pvec1),VINDEX(*pvec2)));

  if (SkipV)
  {
    PRINTDEBUG(gm,1,(" vecskip %4d %4d\n",
                     VECSKIP(*pvec1),VECSKIP(*pvec2)));
    if (VECSKIP(*pvec1) && !VECSKIP(*pvec2))
      if (SkipV==GM_PUT_AT_BEGIN)
        return (-1);
      else
        return ( 1);

    if (VECSKIP(*pvec2) && !VECSKIP(*pvec1))
      if (SkipV==GM_PUT_AT_BEGIN)
        return ( 1);
      else
        return (-1);
  }
  VectorPosition(*pvec1,pv1);
  VectorPosition(*pvec2,pv2);

  V_DIM_SUBTRACT(pv2,pv1,diff);
  V_DIM_SCALE(InvMeshSize,diff);

        #ifdef __THREEDIM__
  PRINTDEBUG(gm,1,(" diff %f %f %f\n",diff[0],diff[1],diff[2]));
        #endif

  if (fabs(diff[Order[DIM-1]])<=ORDERRES)
  {
                #ifdef __THREEDIM__
    if (fabs(diff[Order[DIM-2]])<=ORDERRES)
    {
      if (diff[Order[DIM-3]]>0.0) return (-Sign[DIM-3]);
      else return ( Sign[DIM-3]);
    }
    else
                #endif
    if (diff[Order[DIM-2]]>0.0) return (-Sign[DIM-2]);
    else return ( Sign[DIM-2]);
  }
  else
  {
    if (diff[Order[DIM-1]]>0.0) return (-Sign[DIM-1]);
    else return ( Sign[DIM-1]);
  }
}

static INT PolarLexCompare (VECTOR **pvec1, VECTOR **pvec2)
{
  DOUBLE_VECTOR pv1,pv2;
  DOUBLE norm1,norm2,s1,s2,c1,c2;
  DOUBLE_VECTOR polar;

  PRINTDEBUG(gm,1,("%d: LexCompare %4d %4d",
                   me,VINDEX(*pvec1),VINDEX(*pvec2)));

  if (SkipV)
  {
    PRINTDEBUG(gm,1,(" vecskip %4d %4d\n",
                     VECSKIP(*pvec1),VECSKIP(*pvec2)));
    if (VECSKIP(*pvec1) && !VECSKIP(*pvec2))
      if (SkipV==GM_PUT_AT_BEGIN)
        return (-1);
      else
        return ( 1);

    if (VECSKIP(*pvec2) && !VECSKIP(*pvec1))
      if (SkipV==GM_PUT_AT_BEGIN)
        return ( 1);
      else
        return (-1);
  }
  VectorPosition(*pvec1,pv1);
  VectorPosition(*pvec2,pv2);

  /* first check if one pos is (0,0) */
  V_DIM_EUKLIDNORM(pv1,norm1);
  if (fabs(norm1)<=SMALL_C)
    return (-SignRad);
  V_DIM_EUKLIDNORM(pv2,norm2);
  if (fabs(norm2)<=SMALL_C)
    return ( SignRad);

  /* determine radius difference polar[0] and polar angle difference polar[1] */

  /* radial difference between pv1 and pv2 */
  polar[0] = norm1-norm2;

  /* sin and cos of polar angle between pv1 and unit x-vector */
  V_DIM_SCALE(1./norm1,pv1);
  V_DIM_SCALE(1./norm2,pv2);
  V_DIM_SCALAR_PRODUCT(pv1,unit_vec[1],s1);
  V_DIM_SCALAR_PRODUCT(pv2,unit_vec[1],s2);
  if (s1*s2>=0)
  {
    V_DIM_SCALAR_PRODUCT(pv1,unit_vec[0],c1);
    V_DIM_SCALAR_PRODUCT(pv2,unit_vec[0],c2);
    polar[1] = (c1-c2);
    if (s1==0)
      polar[1] *= SIGNUM(s2);
    else
      polar[1] *= SIGNUM(s1);
  }
  else
    polar[1] = s1-s2;

  if (fabs(polar[Order[DIM-1]])<=ORDERRES)
  {
    if (polar[Order[DIM-2]]>0.0) return (-Sign[DIM-2]);
    else return ( Sign[DIM-2]);
  }
  else
  {
    if (polar[Order[DIM-1]]>0.0) return (-Sign[DIM-1]);
    else return ( Sign[DIM-1]);
  }
}

static int MatrixCompare (MATRIX **MatHandle1, MATRIX **MatHandle2)
{
  INT IND1,IND2;

  IND1 = VINDEX(MDEST(*MatHandle1));
  IND2 = VINDEX(MDEST(*MatHandle2));

  if (IND1>IND2)
    return ( 1);
  else
    return (-1);
}


/****************************************************************************/
/*D
   LexOrderVectorsInGrid - order vectors lexicographically

   SYNOPSIS:
   INT LexOrderVectorsInGrid (GRID *theGrid, const INT *order, const INT *sign,
                                INT SpecSkipVecs, INT AlsoOrderMatrices)

   PARAMETERS:
   .  theGrid - grid level
   .  order - hierarchie of ordering directions (x,y[,z])
   .  sign  - signs for the directions
   .  SpecSkipVecs - if TRUE: GM_PUT_AT_BEGIN or GM_PUT_AT_END of vectors with skip TRUE
   .  AlsoOrderMatrices - if TRUE order matrices in the same sense

   DESCRIPTION:
   This function orders the vectors of one level lexicographically.
   It has the complexity of qsort which is n*log(n).

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

INT LexOrderVectorsInGrid (GRID *theGrid, INT mode, const INT *order, const INT *sign, INT which, INT SpecSkipVecs, INT AlsoOrderMatrices)
{
  MULTIGRID *theMG;
  VECTOR **table,*theVec;
  MATRIX *theMat,**MatTable;
  BVP *theBVP;
  BVP_DESC *theBVPDesc;
  INT i,entries,nm;
  HEAP *theHeap;
  INT takeSkip, takeNonSkip;
  INT MarkKey;

  theMG   = MYMG(theGrid);
  entries = NVEC(theGrid);

  /* calculate the diameter of the bounding rectangle of the domain */
  theBVP = MG_BVP(theMG);
  theBVPDesc = MG_BVPD(MYMG(theGrid));
  InvMeshSize = POW2(GLEVEL(theGrid)) * pow(NN(GRID_ON_LEVEL(theMG,0)),1.0/DIM) / BVPD_RADIUS(theBVPDesc);

  /* allocate memory for the node list */
  if (which==0) return (99);
  takeSkip        = which & GM_TAKE_SKIP;
  takeNonSkip     = which & GM_TAKE_NONSKIP;
  entries = 0;
  for (theVec=FIRSTVECTOR(theGrid); theVec!=NULL; theVec=SUCCVC(theVec))
    if ((takeSkip && VECSKIP(theVec)) || (takeNonSkip && !VECSKIP(theVec)))
      entries++;
  if (entries < 2) return(0);
  theHeap = MGHEAP(theMG);
  MarkTmpMem(theHeap,&MarkKey);
  if ((table = (VECTOR **)GetTmpMem(theHeap,entries*sizeof(VECTOR *),MarkKey))==NULL)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    PrintErrorMessage('E',"LexOrderVectorsInGrid",
                      "could not allocate memory from the MGHeap");
    return (2);
  }
  /* fill array of pointers to nodes */
  entries = 0;
  for (theVec=FIRSTVECTOR(theGrid); theVec!=NULL; theVec=SUCCVC(theVec))
    if ((takeSkip && VECSKIP(theVec))||(takeNonSkip && !VECSKIP(theVec))) {
      VINDEX(theVec) = entries;
      table[entries++] = theVec;
    }
  /* sort array of pointers */
  Order = order;
  Sign  = sign;
  SkipV = SpecSkipVecs;

  if (mode==OV_POLAR)
  {
    if (Order[0]==0)
      SignRad = Sign[0];
    else
      SignRad = Sign[1];
    qsort(table,entries,sizeof(*table),
          (int (*)(const void *, const void *))PolarLexCompare);
  }
  else
  {
    qsort(table,entries,sizeof(*table),
          (int (*)(const void *, const void *))LexCompare);
  }

  for (i=0; i<entries; i++)
    GRID_UNLINK_VECTOR(theGrid,table[i]);

  for (i=0; i<entries; i++) {
    VINDEX(table[i]) = i;
        #ifdef ModelP
    GRID_LINK_VECTOR(theGrid,table[i],PRIO(table[i]));
        #else
    GRID_LINK_VECTOR(theGrid,table[i],PrioMaster);
            #endif
  }
  if (!AlsoOrderMatrices) {
    ReleaseTmpMem(theHeap,MarkKey);
    return (0);
  }
  MatTable = (MATRIX **) table;
  /* now we also order the matrices of each vector the same way
     (just using the new INDICES) */
  /* but leaving the diag at its place of course */
  for (theVec=FIRSTVECTOR(theGrid); theVec!=NULL; theVec=SUCCVC(theVec))
  {
    /* fill array for qsort */
    for (nm=0, theMat=VSTART(theVec); theMat!=NULL; theMat=MNEXT(theMat))
      MatTable[nm++] = theMat;
    if (nm < 2) continue;
    qsort(MatTable+1,nm-1,sizeof(*MatTable),
          (int (*)(const void *, const void *))MatrixCompare);

    /* establish pointer connections */
    MNEXT(MatTable[--nm]) = NULL;
    while (nm>0)
    {
      MNEXT(MatTable[nm-1]) = MatTable[nm];
      --nm;
    }
    VSTART(theVec) = MatTable[0];
  }
  ReleaseTmpMem(theHeap,MarkKey);

  return (0);
}

/****************************************************************************/
/*D
   CreateFindCutProc - Create a new find cut procedure in environement

   SYNOPSIS:
   ALG_DEP *CreateFindCutProc (char *name, FindCutProcPtr FindCutProc);

   PARAMETERS:
   .  name - name
   .  FindCutProc -  the find cut procedure

   DESCRIPTION:
   This function creates a new find cut dependency in environement.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

FIND_CUT *CreateFindCutProc (char *name, FindCutProcPtr FindCutProc)
{
  FIND_CUT *newFindCut;

  if (ChangeEnvDir("/FindCut")==NULL)
  {
    UserWrite("cannot change to dir '/FindCut'\n");
    return(NULL);
  }
  newFindCut = (FIND_CUT *) MakeEnvItem (name,theFindCutVarID,sizeof(FIND_CUT));
  if (newFindCut==NULL) return(NULL);

  newFindCut->FindCutProc = FindCutProc;

  return (newFindCut);
}

/****************************************************************************/
/*D
   CreateAlgebraicDependency - Create a new algebraic dependency in environement

   SYNOPSIS:
   ALG_DEP *CreateAlgebraicDependency (char *name,
   DependencyProcPtr DependencyProc);

   PARAMETERS:
   .  name - name
   .  DependencyProc -  the dependency procedure

   DESCRIPTION:
   This function creates a new algebraic dependency in environement.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

ALG_DEP *CreateAlgebraicDependency (char *name, DependencyProcPtr DependencyProc)
{
  ALG_DEP *newAlgDep;

  if (ChangeEnvDir("/Alg Dep")==NULL)
  {
    UserWrite("cannot change to dir '/Alg Dep'\n");
    return(NULL);
  }
  newAlgDep = (ALG_DEP *) MakeEnvItem (name,theAlgDepVarID,sizeof(ALG_DEP));
  if (newAlgDep==NULL) return(NULL);

  newAlgDep->DependencyProc = DependencyProc;

  return (newAlgDep);
}

/****************************************************************************/
/*
   FeedbackVertexVectors - put `every` vector '!VCUSED' into the new order list

   SYNOPSIS:
   static INT FeedbackVertexVectors (GRID *theGrid, VECTOR **CutVectors,
   INT *nb);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  CutVectors -
   .  nb -

   DESCRIPTION:
   This function can be used as a find-cut-procedure for the streamwise ordering
   algorithm. But it just puts `all` the remaining vectors into the new order
   list instead of removing only as much as necessary to get rid of cyclic
   dependencies.

   It is used as default when no cut procedure is specified.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static VECTOR *FeedbackVertexVectors (GRID *theGrid, VECTOR *LastVector, INT *nb)
{
  VECTOR *theVector;

  /* push all remaining vectors */
  *nb = 0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
    if (!VCUSED(theVector))
    {
      (*nb)++;
      PREDVC(LastVector) = theVector;
      LastVector = theVector;
      SETVCUSED(theVector,1);
    }
  return (LastVector);
}

/****************************************************************************/
/*
   OrderMatrices - reorder the matrix list of a vector in a circular order (2D only)

   SYNOPSIS:
   static INT OrderMatrices (VECTOR *vec, INT Sense)

   PARAMETERS:
   .  vec - order matrix list of this vector
   .  sense - MATHPOS or MATHNEG

   DESCRIPTION:
   This function reorders the matrix list of a vector in a circular order (2D only).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

#ifdef __TWODIM__

static DOUBLE_VECTOR mypos;
static INT sense;

static INT SensCompare (MATRIX **MatHandle1, MATRIX **MatHandle2)
{
  DOUBLE dx1,dy1,dx2,dy2;
  DOUBLE_VECTOR pos1,pos2;

  VectorPosition(MDEST(*MatHandle1),pos1);
  VectorPosition(MDEST(*MatHandle2),pos2);
  dx1 = pos1[0] - mypos[0];
  dy1 = pos1[1] - mypos[1];
  dx2 = pos2[0] - mypos[0];
  dy2 = pos2[1] - mypos[1];

  if (dy1>=0)
  {
    if (dy2<0)
      return (-sense);

    if ((dy1==0) && (dy2==0))
    {
      if (dx1>dx2)
        return (-sense);
      else
        return (sense);
    }

    /* both neighbours in upper half plane */
    if ((dy1*dx2-dx1*dy2)<0)
      return (-sense);
    else
      return (sense);
  }
  else
  {
    if (dy2>=0)
      return (sense);

    /* both neighbours in lower half plane */
    if ((dy1*dx2-dx1*dy2)<0)
      return (-sense);
    else
      return (sense);
  }
}

static INT OrderMatrices (VECTOR *vec, INT Sense)
{
  VECTOR *nbv;
  MATRIX *mat,*MatList[MATTABLESIZE];
  INT i,flag,condition1,condition2,start,nm;

  /* fill matrices without diag into list */
  for (nm=0, mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
  {
    if (nm>=MATTABLESIZE)
      return (1);

    MatList[nm++] = mat;
  }
  VectorPosition(vec,mypos);

  sense = Sense;
  qsort(MatList,nm,sizeof(MATRIX*),(int (*)(const void *, const void *))SensCompare);

  /* find 'begin' */
  flag = FALSE;
  for (start=0; start<nm; start++)
  {
    nbv = MDEST(MatList[start]);
    condition1 = ((VCLASS(nbv)<ACTIVE_CLASS) || VCUSED(nbv));
    condition2 = (condition1 || (OBJT(MYVERTEX((NODE*)VOBJECT(nbv)))==BVOBJ));
    if (flag && !condition1)
      break;
    else if (condition2)
      flag = TRUE;
  }

  /* establish pointer connections */
  for (i=start; i<nm+start-1; i++)
    MNEXT(MatList[i%nm]) = MatList[(i+1)%nm];
  MNEXT(MatList[i%nm]) = NULL;
  MNEXT(VSTART(vec)) = MatList[start%nm];

  return (0);
}

#endif

/****************************************************************************/
/*
   ShellOrderVectors	- Reorder double linked vector list by a shell algorithm

   SYNOPSIS:
   INT ShellOrderVectors (GRID *theGrid, VECTOR *seed)

   PARAMETERS:
   .  theGrid - pointer to grid
   .  seed - start at this vector

   DESCRIPTION:
   This function reorders double linked vector list by a shell algorithm.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

INT ShellOrderVectors (GRID *theGrid, VECTOR *seed)
{
  FIFO myfifo;
  void *buffer;
  VECTOR **vlist;
  VECTOR *theV;
  MATRIX *theM;
  HEAP *theHeap=MGHEAP(MYMG(theGrid));
  INT i,n;
  INT MarkKey;

  /* count vectors */
  n = 0;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    n++;

  if (n == 0)
    return(0);

  /* find new ordering */
  MarkTmpMem(theHeap,&MarkKey);
  buffer=(void *)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);
  vlist = (VECTOR**)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);
  fifo_init(&myfifo,buffer,sizeof(VECTOR*)*n);
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    SETVCUSED(theV,0);
  fifo_in(&myfifo,(void *)seed);
  SETVCUSED(seed,1);
  i = 0;
  while(!fifo_empty(&myfifo))
  {
    theV = (VECTOR *)fifo_out(&myfifo);
    vlist[i++] = theV;

                #ifdef __TWODIM__
    OrderMatrices(theV,MATHPOS);
                #endif

    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      if (!CEXTRA(MMYCON(theM)))
        if (!VCUSED(MDEST(theM)))
        {
          fifo_in(&myfifo,(void *)MDEST(theM));
          SETVCUSED(MDEST(theM),1);
        }
  }
  assert(i==n);

  /* reorder vector list in grid */
  for (i=0; i<n; i++) GRID_UNLINK_VECTOR(theGrid,vlist[i]);
  for (i=0; i<n; i++) GRID_LINK_VECTOR(theGrid,vlist[i],PrioMaster);

  ReleaseTmpMem(theHeap,MarkKey);

  return (0);
}

/****************************************************************************/
/*
   OrderVectorAlgebraic	- Reorder double linked node list by a streamline ordering

   SYNOPSIS:
   static INT OrderVectorAlgebraic (GRID *theGrid, INT mode, INT putSkipFirst, INT skipPat);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  mode
   .  putSkipFirst - if TRUE put vectors with a skip pattern larger than SkipPat to begin of list
   .  skipPat - s.a.

   DESCRIPTION:
   This function reorders double linked node list by a streamline ordering.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT CheckConsistence (GRID *theGrid, char *location)
{
  VECTOR *theVector;
  MATRIX *theMatrix;
  INT up, down, nfound;

  /* check USED, N_INFLOW and N_OUTFLOW */
  nfound = 0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    if (VCUSED(theVector)) continue;

    /* count upward and downward matrices */
    up = down = 0;
    for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      if (VCUSED(MDEST(theMatrix))) continue;
      if (MUP(theMatrix)) up++;
      if (MDOWN(theMatrix)) down++;
    }
    if (VUP(theVector)!=up || VDOWN(theVector)!=down)
    {
      nfound++;
      if (nfound==1 && location!=NULL)
        UserWriteF("Failed in %s\n",location);
      UserWriteF("vector = %d: up: %d %d   down: %d %d\n",(int)VINDEX(theVector),(int)VUP(theVector),(int)up,
                 (int)VDOWN(theVector),(int)down);
    }
  }
  return (nfound);
}

static INT CheckVectorList (GRID *theGrid)
{
  INT i;
  VECTOR *theVector;

  /* check # members of succ list */
  i=0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector)) i++;
  if (NVEC(theGrid) != i)
    RETURN (1);

  /* check # members of pred list */
  i=0;
  for (theVector=LASTVECTOR(theGrid); theVector!= NULL; theVector=PREDVC(theVector)) i++;
  if (NVEC(theGrid) != i)
    RETURN (1);

  return (0);
}

static INT CheckBVList (GRID *theGrid)
{
  BLOCKVECTOR *theBV;

  /* check # members of succ list */
  for (theBV=GFIRSTBV(theGrid); BVSUCC(theBV)!= NULL; theBV=BVSUCC(theBV))
    if (BVENDVECTOR(theBV)!=BVFIRSTVECTOR(BVSUCC(theBV)))
      RETURN (1);

  return (0);
}

static void GetBVNumberInGrid (GRID *theGrid, INT *nb)
{
  VECTOR *theVector;

  nb[0]=0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector)) nb[0]++;
  nb[1]=0;
  for (theVector=LASTVECTOR(theGrid); theVector!= NULL; theVector=PREDVC(theVector)) nb[1]++;
  return;
}

static void GetBVNumber (BLOCKVECTOR *theBV, INT *nb)
{
  VECTOR *theVector;

  nb[0]=nb[1]=0;
  for (theVector=BVFIRSTVECTOR(theBV); theVector!= BVENDVECTOR(theBV); theVector=SUCCVC(theVector)) nb[0]++;
  for (theVector=BVLASTVECTOR(theBV); theVector!= PREDVC(BVFIRSTVECTOR(theBV)); theVector=PREDVC(theVector)) nb[1]++;
  return;
}

static INT OrderVectorAlgebraic (GRID *theGrid, INT mode, INT putSkipFirst, INT skipPat)
{
  BLOCKVECTOR *theFirstBV, *theLastBV, *theCutBV, *theBV, *takeOut, **FBVList, **LBVList, **CBVList;
  VECTOR FIRST_handle,LAST_handle;
  VECTOR *FIRST_next_out,*FIRST_last_in;
  VECTOR *LAST_next_out,*LAST_last_in;
  VECTOR *CUT_begin;
  VECTOR *theVector,*theNbVector,*succVector,*predVector, *theFirstVector;
  MATRIX *theMatrix;
  DOUBLE a;
  INT i,k,cycle,nCutTot;
  INT nFIRST,nCUT,nLAST;
  INT up, down;
  HEAP *theHeap;

        #ifdef ModelP
  ASSERT(FALSE);       /* see TODO below */
        #endif

  /********************************************************************/
  /*	init				                                                                                        */
  /********************************************************************/

  /* cancel all BLOCKVECTORS */
  FreeAllBV(theGrid);
  SETUSED(theGrid,0);

  /* init USED, N_INFLOW and N_OUTFLOW */
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    /* reset used flag (indicating the vector hasn't yet its new place in the list) */
    SETVCUSED(theVector,0);

    /* count upward and downward matrices */
    up = down = 0;
    for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      if (MUP(theMatrix)) up++;
      if (MDOWN(theMatrix)) down++;
    }
    SETVUP(theVector,up);
    SETVDOWN(theVector,down);
  }

  /* in the sequel we use (and therefore destroy) the PREDVC-list for book keeping */

  cycle = nCutTot = 0;

  /* init pointers and set the first FIRST and LAST set */
  FIRST_next_out = FIRST_last_in = &FIRST_handle;
  PREDVC(FIRST_last_in) = NULL; nFIRST = 0;
  LAST_next_out  = LAST_last_in  = &LAST_handle;
  PREDVC(LAST_last_in) = NULL; nLAST  = 0;
  nCUT   = 0;
  theFirstBV = theLastBV = theCutBV = NULL;
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    if (((putSkipFirst) && (VECSKIP(theVector) & skipPat == skipPat)) ||
        (VUP(theVector)==0))
    {
      /* append to FIRST list */
      nFIRST++;
      PREDVC(FIRST_last_in) = theVector;
      FIRST_last_in = theVector;
      SETVCUSED(theVector,1);
      VINDEX(theVector) = 0;
    }
    else if (VDOWN(theVector)==0)
    {
      /* append to last list */
      nLAST++;
      PREDVC(LAST_last_in) = theVector;
      LAST_last_in = theVector;
      SETVCUSED(theVector,1);
      VINDEX(theVector) = 1;
    }
  }
  PREDVC(FIRST_last_in) = PREDVC(LAST_last_in) = NULL;

  /* store in BLOCKVECTOR */
  if (nFIRST>0)
  {
    if (CreateBlockvector_l0(theGrid,&theBV,NULL,1)!=GM_OK) RETURN (1);
    theFirstBV = theBV;
    BVNUMBER(theFirstBV)    = 3*cycle+0;
    BVFIRSTVECTOR(theFirstBV)       = BVPRED(&FIRST_handle);
  }

  if (nLAST>0)
  {
    if (CreateBlockvector_l0(theGrid,&theBV,theFirstBV,1)!=GM_OK) RETURN (1);
    theLastBV = theBV;
    BVNUMBER(theLastBV)             = 3*cycle+1;
    BVFIRSTVECTOR(theLastBV)        = BVPRED(&LAST_handle);
  }

  do
  {
    cycle++;

    /****************************************************************************/
    /*	find next FIRST-set in vectors not used                                                                 */
    /****************************************************************************/

    nFIRST = nLAST = 0;

    theFirstVector  = FIRST_last_in;
    for (FIRST_next_out=PREDVC(FIRST_next_out); FIRST_next_out!=NULL; FIRST_next_out=PREDVC(FIRST_next_out))
      for (theMatrix=MNEXT(VSTART(FIRST_next_out)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
      {
        theNbVector = MDEST(theMatrix);
        if (MDOWN(theMatrix) && !VCUSED(theNbVector))
        {
          k = VUP(theNbVector);
          assert(k>0);                                                  /* if 0 is supposed to be VCUSED already */
          SETVUP(theNbVector,--k);
          if (k==0)
          {
            /* this vector has only matrices going down */

            nFIRST++;
            /* append to FIRST list */
            PREDVC(FIRST_last_in) = theNbVector;
            FIRST_last_in = theNbVector;
            PREDVC(FIRST_last_in) = NULL;
            SETVCUSED(theNbVector,1);
            VINDEX(theNbVector) = 3*cycle;
          }
        }
        if ((nCUT>0) && !VCUSED(theNbVector))
          /* we also have to push LAST vectors */
          if (MUP(theMatrix))
          {
            k = VDOWN(theNbVector);
            /*assert(k>0);*/			/* if 0 is supposed to be VCUSED already */
            if (k<=0)
              RETURN (1);
            SETVDOWN(theNbVector,--k);
            if (k==0)
            {
              /* this vector has only matrices going up */

              nLAST++;
              /* append to last list */
              PREDVC(LAST_last_in) = theNbVector;
              LAST_last_in = theNbVector;
              PREDVC(LAST_last_in) = NULL;
              SETVCUSED(theNbVector,1);
              VINDEX(theNbVector) = 3*cycle+1;
            }
          }
      }
    FIRST_next_out = FIRST_last_in;

    /* create first BLOCKVECTOR */
    if (nFIRST>0)
    {
      if (CreateBlockvector_l0(theGrid,&theBV,theFirstBV,1)!=GM_OK) RETURN (1);
      theFirstBV = theBV;
      BVNUMBER(theFirstBV)            = 3*cycle+0;
      BVFIRSTVECTOR(theFirstBV)       = PREDVC(theFirstVector);
    }

    /********************************************************************/
    /*	find next LAST-set in vectors not used			                        */
    /********************************************************************/

    for (LAST_next_out=PREDVC(LAST_next_out); LAST_next_out!=NULL; LAST_next_out=PREDVC(LAST_next_out))
      for (theMatrix=MNEXT(VSTART(LAST_next_out)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
      {
        theNbVector = MDEST(theMatrix);
        if (!VCUSED(theNbVector) && MUP(theMatrix))
        {
          k = VDOWN(theNbVector);
          assert(k>0);                                                  /* if 0 is supposed to be VCUSED already */
          SETVDOWN(theNbVector,--k);
          if (k==0)
          {
            /* this vector has only matrices going up */

            nLAST++;
            /* append to last list */
            PREDVC(LAST_last_in) = theNbVector;
            LAST_last_in = theNbVector;
            PREDVC(LAST_last_in) = NULL;
            SETVCUSED(theNbVector,1);
            VINDEX(theNbVector) = 3*cycle+1;
          }
        }
      }
    LAST_next_out = LAST_last_in;

    /* create Last BLOCKVECTOR */
    if (nLAST>0)
    {
      if (CreateBlockvector_l0(theGrid,&theBV,theLastBV,0)!=GM_OK) RETURN (1);
      theLastBV = theBV;
      BVNUMBER(theLastBV)             = 3*cycle+1;
      BVFIRSTVECTOR(theLastBV)        = LAST_last_in;
    }

    /****************************************************************************/
    /*	get CUT (or Feedback Vertex)-set and do what needs to be done			*/
    /****************************************************************************/

    CUT_begin = FIRST_last_in;
    FIRST_last_in = (*FindCutSet)(theGrid,FIRST_last_in,&nCUT);
    if (FIRST_last_in==NULL) nCUT = 0;
    else PREDVC(FIRST_last_in) = NULL;

    /* create cut BLOCKVECTOR */
    if (nCUT>0)
    {
      if (CreateBlockvector_l0(theGrid,&theBV,theFirstBV,1)!=GM_OK) RETURN (1);
      theCutBV = theBV;
      theFirstBV = theBV;
      BVNUMBER(theCutBV)                      = 3*(cycle+1)+2;
      BVFIRSTVECTOR(theCutBV)         = PREDVC(CUT_begin);
    }

    /* set index */
    for (theVector=PREDVC(CUT_begin); theVector!=NULL; theVector=PREDVC(theVector))
      VINDEX(theVector) = 3*(cycle+1)+2;

    nCutTot += nCUT;

  } while (nCUT>0);

  UserWriteF("# %d cycles: %d cutted from %d\n",(int)cycle,(int)nCutTot,(int)NVEC(theGrid));
  a = POW((DOUBLE)NVEC(theGrid),(DOUBLE)(DIM-1)/(DOUBLE)DIM);     a = (DOUBLE)nCutTot/a;
  UserWriteF("# corr. to %6.2f hyp. planes\n",(float)a);

  /* TODO: use dlmgr macros for lists after streamwise ordering */

  /* insert FIRST list one-by-one to LASTVECTOR list of grid */
  LASTVECTOR(theGrid) = NULL;
  for (theVector=PREDVC(&FIRST_handle); theVector!=NULL; theVector=predVector)
  {
    predVector                      = PREDVC(theVector);
    PREDVC(theVector)   = LASTVECTOR(theGrid);
    LASTVECTOR(theGrid) = theVector;
  }

  /* insert LAST list as a whole to LASTVECTOR list of grid */
  PREDVC(LAST_last_in) = LASTVECTOR(theGrid);
  LASTVECTOR(theGrid)  = PREDVC(&LAST_handle);

  /* construct SUCCVC list */
  succVector = NULL;
  for (theVector=LASTVECTOR(theGrid); theVector!=NULL; theVector=PREDVC(theVector))
  {
    SUCCVC(theVector) = succVector;
    succVector = theVector;
  }
  SFIRSTVECTOR(theGrid) = succVector;
  PREDVC(succVector)   = NULL;

  /* set pointers in BLOCKVECTORs */
  BVLASTVECTOR(GLASTBV(theGrid)) = LASTVECTOR(theGrid);
  for (theBV=GLASTBV(theGrid); theBV!=NULL; theBV=BVPRED(theBV))
  {
    if (BVSUCC(theBV)!=NULL && BVLASTVECTOR(theBV)==NULL)
      BVLASTVECTOR(theBV)=PREDVC(BVFIRSTVECTOR(BVSUCC(theBV)));
    if (BVFIRSTVECTOR(theBV)==NULL)
    {
      assert(0);
      BVFIRSTVECTOR(theBV)=BVENDVECTOR(theBV);
    }
  }

  if (mode==GM_FFLLCC)
  {
    for (theBV=GLASTBV(theGrid); BVPRED(theBV)!=NULL; theBV=BVPRED(theBV))
    {
      takeOut = BVPRED(theBV);

      if (BVNUMBER(takeOut)%3==2)
      {
        if (CutBlockvector_l0(theGrid,takeOut,YES))
          RETURN (1);

        if (InsertBlockvector_l0(theGrid,takeOut,NULL,0,YES))
          RETURN (1);

        theBV = BVSUCC(theBV);
      }
    }
  }
  else if (mode==GM_FFLCLC)
  {
    INT MarkKey;

    theHeap = MYMG(theGrid)->theHeap;
    MarkTmpMem(theHeap,&MarkKey);
    if ((FBVList=(BLOCKVECTOR **)GetTmpMem(theHeap,3*(cycle+1)*sizeof(BLOCKVECTOR*),MarkKey))==NULL) RETURN (1);
    for (i=0; i<3*cycle+3; i++) FBVList[i]=NULL;
    LBVList = FBVList + cycle + 1;
    CBVList = LBVList + cycle + 1;

    for (theBV=GFIRSTBV(theGrid); theBV!=NULL; theBV=BVSUCC(theBV))
      if (BVNUMBER(theBV)%3==0)
        FBVList[BVNUMBER(theBV)/3] = theBV;
      else if (BVNUMBER(theBV)%3==1)
        LBVList[BVNUMBER(theBV)/3] = theBV;
      else
        CBVList[BVNUMBER(theBV)/3] = theBV;

    for (i=cycle; i>0; i--)
    {
      if (CBVList[i]==NULL) continue;
      if (CutBlockvector_l0(theGrid,CBVList[i],YES)) RETURN (1);
      k=i-1; while (LBVList[k]==NULL && k>0) k--;
      if (InsertBlockvector_l0(theGrid,CBVList[i],LBVList[k],0,YES)) RETURN (1);
    }

    ReleaseTmpMem(theHeap,MarkKey);
  }
  else if (mode==GM_CCFFLL)
  {
    for (theBV=GFIRSTBV(theGrid); BVSUCC(theBV)!=NULL; theBV=BVSUCC(theBV))
    {
      takeOut = BVSUCC(theBV);

      if (BVNUMBER(takeOut)%3==2)
      {
        if (CutBlockvector_l0(theGrid,takeOut,YES))
          RETURN (1);

        if (InsertBlockvector_l0(theGrid,takeOut,NULL,1,YES))
          RETURN (1);

        theBV = BVPRED(theBV);
      }
    }
  }


  /* set VCCUT */
  for (theBV=GFIRSTBV(theGrid); theBV!=NULL; theBV=BVSUCC(theBV))
    if (BV_GEN(theBV)==BV_GEN_C)
      for (theVector=BVFIRSTVECTOR(theBV); theVector!=BVENDVECTOR(theBV); theVector=SUCCVC(theVector))
      {
        SETVCCUT(theVector,1);
      }
    else
      for (theVector=BVFIRSTVECTOR(theBV); theVector!=BVENDVECTOR(theBV); theVector=SUCCVC(theVector))
      {
        SETVCCUT(theVector,0);
      }

  /* check # members of succ list */
  i=0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector)) i++;
  if (NVEC(theGrid) != i)
  {
    UserWrite("vectorstructure corrupted\n");
    RETURN (1);
  }

  /* check # members of pred list */
  i=0;
  for (theVector=LASTVECTOR(theGrid); theVector!= NULL; theVector=PREDVC(theVector)) i++;
  if (NVEC(theGrid) != i)
  {
    UserWrite("vectorstructure corrupted\n");
    RETURN (1);
  }

  /* set index w.r.t. new order, beginning with 1 */
  i = 1;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
    VINDEX(theVector) = i++;


  return (0);
}

/****************************************************************************/
/*D
   OrderVectors	-  Driver for general vector ordering

   SYNOPSIS:
   INT OrderVectors (MULTIGRID *theMG, INT levels, INT mode, INT PutSkipFirst, INT SkipPat,
                                        const char *dependency, const char *dep_options, const char *findcut);

   PARAMETERS:
   .  theMG -  multigrid to order
   .  levels -  GM_ALL_LEVELS or GM_CURRENT_LEVEL
   .  mode - GM_FCFCLL or GM_FFCCLL (see 'orderv' command)
   .  PutSkipFirst - if TRUE put vectors with a skip pattern larger than SkipPat to begin of list
   .  SkipPat - s.a.
   .  dependency - name of user defined dependency item
   .  dep_options - options for user dependency function

   DESCRIPTION:
   This function orders 'VECTOR's in a multigrid according to the dependency
   function provided.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR if error occured.
   D*/
/****************************************************************************/

INT OrderVectors (MULTIGRID *theMG, INT levels, INT mode, INT PutSkipFirst, INT SkipPat, const char *dependency, const char *dep_options, const char *findcut)
{
  INT i, currlevel, baselevel;
  ALG_DEP *theAlgDep;
  FIND_CUT *theFindCut;
  GRID *theGrid;
  DependencyProcPtr DependencyProc;

  /* check mode */
  if ((mode!=GM_FCFCLL)&&
      (mode!=GM_FFLLCC)&&
      (mode!=GM_FFLCLC)&&
      (mode!=GM_CCFFLL)) RETURN(GM_ERROR);

  /* current level */
  currlevel = theMG->currentLevel;

  /* get algebraic dependency */
  theAlgDep=NULL;
  if (dependency!=NULL)
  {
    theAlgDep = (ALG_DEP *) SearchEnv(dependency,"/Alg Dep",theAlgDepVarID,theAlgDepDirID);
    if (theAlgDep==NULL)
    {
      UserWrite("algebraic dependency not found\n");
      RETURN(GM_ERROR);
    }
    DependencyProc = theAlgDep->DependencyProc;
    if (DependencyProc==NULL)
    {
      UserWrite("don't be stupid: implement a dependency!\n");
      RETURN(GM_ERROR);
    }
  }

  /* get find cut dependency */
  if (findcut==NULL)
  {
    FindCutSet = FeedbackVertexVectors;
    UserWrite("default cut set proc:\n    leaving order of cyclic dependencies unchanged\n");
  }
  else
  {
    theFindCut = (FIND_CUT *) SearchEnv(findcut,"/FindCut",theFindCutVarID,theFindCutDirID);
    if (theFindCut==NULL)
    {
      UserWrite("find cut proc not found\n");
      RETURN(GM_ERROR);
    }
    FindCutSet = theFindCut->FindCutProc;
    if (FindCutSet==NULL)
    {
      UserWrite("don't be stupid: implement a find cut proc!\n");
      RETURN(GM_ERROR);
    }
  }

  /* go */
  if (levels==GM_ALL_LEVELS)
    baselevel = 0;
  else
    baselevel = currlevel;
  if (theAlgDep!=NULL)
    for (i=baselevel; i<=currlevel; i++)
    {
      theGrid = GRID_ON_LEVEL(theMG,i);
      if ((*DependencyProc)(theGrid,dep_options)) RETURN(GM_ERROR);
    }
  for (i=baselevel; i<=currlevel; i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if (OrderVectorAlgebraic(theGrid,mode,PutSkipFirst,SkipPat)) RETURN(GM_ERROR);
  }

  return (GM_OK);
}

static VECTOR *CutStrongLine (VECTOR *theVector, VECTOR *AppendVec)
{
  MATRIX *mat;

  do
  {
    /* append to list */
    PREDVC(AppendVec) = theVector;
    AppendVec = theVector;
    PREDVC(AppendVec) = NULL;
    SETVCUSED(theVector,1);

    /* look for strong nb */
    for (mat=MNEXT(VSTART(theVector)); mat!=NULL; mat=MNEXT(mat))
      if (MSTRONG(mat))
      {
        theVector = MDEST(mat);

        if (!VCUSED(theVector))
          break;
      }
  }
  while (mat!=NULL);

  return (AppendVec);
}

static VECTOR *FindOptimalStrong (FIFO *fifo)
{
  VECTOR *first,*vec,*nbvec,*minvec;
  MATRIX *mat;
  INT nStrong,minStrong,nUp,minUp;

  /* drop leading vectors used */
  while ((vec=(VECTOR*)fifo_out(fifo))!=NULL)
    if (!VCUSED(vec))
      break;

  if (vec==NULL)
    return (NULL);

  first = vec;

  /* determine minimal number of !used strong nbs */
  minStrong = MAX_I;
  do
    if (!VCUSED(vec))
    {
      nStrong = 0;
      for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
        if (MSTRONG(mat))
        {
          nbvec = MDEST(mat);
          if (!VCUSED(nbvec))
            nStrong++;
        }
      fifo_in(fifo,vec);

      minStrong = MIN(minStrong,nStrong);
    }
  while ((vec=(VECTOR*)fifo_out(fifo))!=first);

  /* determine vec with nStrong == minStrong and minimal number of up nbs */
  /* NB: fifo contains no VCUSED(vec) anymore */
  minUp = MAX_I;
  do
  {
    nUp = nStrong = 0;
    for (mat=MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat))
      if (MSTRONG(mat))
      {
        nbvec = MDEST(mat);
        if (!VCUSED(nbvec))
          nStrong++;
      }
      else if (MUP(mat))
      {
        nbvec = MDEST(mat);
        if (!VCUSED(nbvec))
          nUp++;
      }

    fifo_in(fifo,vec);

    if (nStrong==minStrong)
      if (nUp<minUp)
      {
        minUp = nUp;
        minvec = vec;
      }
  }
  while ((vec=(VECTOR*)fifo_out(fifo))!=first);

  /* remove minvec from fifo */
  if (vec!=minvec)
  {
    do
    {
      if (vec!=minvec)
        fifo_in(fifo,vec);
    }
    while ((vec=(VECTOR*)fifo_out(fifo))!=first);

    /* push first again */
    fifo_in(fifo,vec);
  }

  return (minvec);
}

/****************************************************************************/
/*
   LineOrderVectorsAlgebraic	- Reorder double linked vector list by a streamline ordering

   SYNOPSIS:
   static INT LineOrderVectorsAlgebraic (GRID *theGrid, INT mode);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function reorders double linked vector list by a streamline ordering.

   CAUTION:
   The 'VCUSED' flag has to be initialized by the AlgDep procedure. If 'VCUSED==1
   then the vector is pushed to the FIRST fifo at initialization.

   The FindCutProc can only use the SUCCVC-list (the other one is destroyed).

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT LineOrderVectorsAlgebraic (GRID *theGrid, INT verboselevel)
{
  BLOCKVECTOR *theFirstBV, *theLastBV, *theCutBV, *theBV;
  FIFO FFifo,LFifo;
  VECTOR FIRST_handle,LAST_handle;
  VECTOR *FIRST_next_out,*FIRST_last_in;
  VECTOR *FIRST_line_start,*FIRST_line_end;
  VECTOR *LAST_line_start,*LAST_line_end;
  VECTOR *LAST_next_out,*LAST_last_in;
  VECTOR *theVector,*theNbVector,*succVector,*predVector;
  MATRIX *theMatrix;
  DOUBLE a;
  INT i,k,init,cycle,line,nCutTot,fifosize;
  INT nCUT;
  INT up, down, strong;
  INT StrongInflow,nInflow,StrongOutflow,nOutflow,pushInflow,pushOutflow;
  INT bvn;
  char gen_label[3];
  INT MarkKey;

        #ifdef ModelP
  ASSERT(FALSE);       /* see TODO below */
        #endif

  gen_label[GM_GEN_FIRST] = 'F';
  gen_label[GM_GEN_LAST]  = 'L';
  gen_label[GM_GEN_CUT]   = 'C';

  /********************************************************************/
  /*	init				                                                                                        */
  /********************************************************************/

  /* dispose all BLOCKVECTORS */
  FreeAllBV(theGrid);

  /* init fifos */
  MarkTmpMem(MGHEAP(MYMG(theGrid)),&MarkKey);
  fifosize = 30*floor(sqrt(NVEC(theGrid)))*sizeof(VECTOR*);
  fifo_init(&FFifo,GetTmpMem(MGHEAP(MYMG(theGrid)),fifosize,MarkKey),fifosize);
  fifo_init(&LFifo,GetTmpMem(MGHEAP(MYMG(theGrid)),fifosize,MarkKey),fifosize);

  /* init USED, N_INFLOW and N_OUTFLOW */
  StrongInflow = StrongOutflow = 0;
  nInflow = nOutflow = 0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    SETVCSTRONG(theVector,0);

    /* reorder matrix lists in math pos sense */
                #ifdef __TWODIM__
    OrderMatrices(theVector,MATHNEG);
                #endif

    /* count upward, downward and strong matrices */
    up = down = strong = 0;
    for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      if (MUP(theMatrix)) up++;
      if (MDOWN(theMatrix)) down++;
      if (MSTRONG(theMatrix)) strong++;
    }
    SETVUP(theVector,up);
    SETVDOWN(theVector,down);
    if (strong)
      SETVCSTRONG(theVector,1);

    if (VCUSED(theVector))
    {
      nInflow++;
      if (strong)
        StrongInflow++;
    }
    if (VCFLAG(theVector))
    {
      nOutflow++;
      if (strong)
        StrongOutflow++;
    }
  }
  pushInflow  = (nInflow ==StrongInflow);
  pushOutflow = (nOutflow==StrongOutflow);

  /* in the sequel we use (and therefore destroy) the PREDVC-list for book keeping */

  cycle = nCutTot = 0;
  init  = 1;

  /* init pointers and set the first FIRST and LAST set */
  FIRST_next_out = FIRST_last_in = &FIRST_handle;
  PREDVC(FIRST_last_in) = NULL;
  LAST_next_out  = LAST_last_in  = &LAST_handle;
  PREDVC(LAST_last_in) = NULL;
  nCUT   = 0;
  theFirstBV = theLastBV = theCutBV = NULL;
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    if (VCUSED(theVector) && pushInflow)
    {
      if (fifo_in(&FFifo,theVector)!=0)
      {
        PrintErrorMessage('E',"LineOrderVectorsAlgebraic","fifo full");
        return (__LINE__);
      }
      SETVCFLAG(theVector,0);
      SETVCUSED(theVector,0);
      continue;
    }
    if (VCFLAG(theVector) && pushOutflow)
    {
      if (fifo_in(&LFifo,theVector)!=0)
      {
        PrintErrorMessage('E',"LineOrderVectorsAlgebraic","fifo full");
        return (__LINE__);
      }
      SETVCFLAG(theVector,0);
      SETVCUSED(theVector,0);
      continue;
    }
    SETVCFLAG(theVector,VCSTRONG(theVector));
    SETVCUSED(theVector,0);
    if (VUP(theVector)==0)
    {
      /* append to FIRST list */
      PREDVC(FIRST_last_in) = theVector;
      FIRST_last_in = theVector;
      SETVCUSED(theVector,1);
    }
    else if (VDOWN(theVector)==0)
    {
      /* append to LAST list */
      PREDVC(LAST_last_in) = theVector;
      LAST_last_in = theVector;
      SETVCUSED(theVector,1);
    }
  }
  PREDVC(FIRST_last_in) = PREDVC(LAST_last_in) = NULL;

  do
  {
    cycle++;

    /****************************************************************************/
    /*	find next FIRST-set in vectors not used                                                                 */
    /****************************************************************************/

    line = 0;
    do
    {
      FIRST_line_start = PREDVC(FIRST_next_out);
      FIRST_line_end   = FIRST_last_in;
      for (FIRST_next_out=PREDVC(FIRST_next_out); FIRST_next_out!=NULL; FIRST_next_out=PREDVC(FIRST_next_out))
      {
        for (theMatrix=MNEXT(VSTART(FIRST_next_out)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
        {
          theNbVector = MDEST(theMatrix);
          if (VCUSED(theNbVector)) continue;

          if (VCFLAG(theNbVector))
          {
            if (fifo_in(&FFifo,theNbVector)!=0)
            {
              PrintErrorMessage('E',"LineOrderVectorsAlgebraic","fifo full");
              return (__LINE__);
            }
            SETVCFLAG(theNbVector,0);
          }
          else if (MDOWN(theMatrix))
          {
            k = VUP(theNbVector);
            assert(k>0);                                                        /* if 0 is supposed to be VCUSED already */
            SETVUP(theNbVector,--k);
            if (k==0)
            {
              /* this vector has only matrices going down */

              /* append to FIRST list */
              PREDVC(FIRST_last_in) = theNbVector;
              FIRST_last_in = theNbVector;
              PREDVC(FIRST_last_in) = NULL;
              SETVCUSED(theNbVector,1);
            }
          }
          else if (MUP(theMatrix))
          {
            k = VDOWN(theNbVector);
            /*assert(k>0);*/			/* if 0 is supposed to be VCUSED already */
            if (k<=0)
            {
              PrintErrorMessage('E',"LineOrderVectorsAlgebraic","DOWN counter < 0");
              return (__LINE__);
            }
            SETVDOWN(theNbVector,--k);
            if (k==0)
            {
              /* this vector has only matrices going up */

              /* append to last list */
              PREDVC(LAST_last_in) = theNbVector;
              LAST_last_in = theNbVector;
              PREDVC(LAST_last_in) = NULL;
              SETVCUSED(theNbVector,1);
            }
          }
        }
        if (FIRST_next_out==FIRST_line_end)
        {
          /* create block line */
          if (CreateBlockvector_l0(theGrid,&theBV,theFirstBV,1)!=GM_OK)
          {
            PrintErrorMessage('E',"LineOrderVectorsAlgebraic","CreateBlockvector_l0 failed");
            return (__LINE__);
          }
          theFirstBV = theBV;
          if (nCUT>0)
          {
            SET_ORD_GCL(BVNUMBER(theFirstBV),GM_GEN_CUT,cycle,line);
            nCUT = 0;
          }
          else if (init && (line==0))
          {
            SET_ORD_GCL(BVNUMBER(theFirstBV),GM_GEN_FIRST,0,line);
          }
          else
          {
            SET_ORD_GCL(BVNUMBER(theFirstBV),GM_GEN_FIRST,cycle,line);
          }
          BVFIRSTVECTOR(theFirstBV)       = FIRST_line_start;

          line++;
          FIRST_line_start = PREDVC(FIRST_next_out);
          FIRST_line_end   = FIRST_last_in;
        }
      }
      FIRST_next_out = FIRST_last_in;

      /* check if there is a !used vec in the strong fifo */
      if ((theVector=FindOptimalStrong(&FFifo))!=NULL)
        FIRST_last_in = CutStrongLine(theVector,FIRST_last_in);
    }
    while (theVector!=NULL);

    /********************************************************************/
    /*	find next LAST-set in vectors not used			                        */
    /********************************************************************/

    line = 0;
    do
    {
      LAST_line_start = PREDVC(LAST_next_out);
      LAST_line_end   = LAST_last_in;
      for (LAST_next_out=PREDVC(LAST_next_out); LAST_next_out!=NULL; LAST_next_out=PREDVC(LAST_next_out))
      {
        for (theMatrix=MNEXT(VSTART(LAST_next_out)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
        {
          theNbVector = MDEST(theMatrix);
          if (VCUSED(theNbVector)) continue;

          if (VCFLAG(theNbVector))
          {
            /* this is a strongly coupled neighbour */
            if (fifo_in(&LFifo,theNbVector)!=0)
              return (__LINE__);
            SETVCFLAG(theNbVector,0);
          }
          else if (MUP(theMatrix))
          {
            k = VDOWN(theNbVector);
            assert(k>0);                                                        /* if 0 is supposed to be VCUSED already */
            SETVDOWN(theNbVector,--k);
            if (k==0)
            {
              /* this vector has only matrices going up */

              /* append to LAST list */
              PREDVC(LAST_last_in) = theNbVector;
              LAST_last_in = theNbVector;
              PREDVC(LAST_last_in) = NULL;
              SETVCUSED(theNbVector,1);
            }
          }
          else if (MDOWN(theMatrix))
          {
            k = VUP(theNbVector);
            assert(k>0);                                                        /* if 0 is supposed to be VCUSED already */
            SETVUP(theNbVector,--k);
            if (k==0)
            {
              /* this vector has only matrices going down */

              /* append to FIRST list */
              PREDVC(FIRST_last_in) = theNbVector;
              FIRST_last_in = theNbVector;
              PREDVC(FIRST_last_in) = NULL;
              SETVCUSED(theNbVector,1);
            }
          }
        }
        if (LAST_next_out==LAST_line_end)
        {
          if (theLastBV==NULL)
          {
            if (CreateBlockvector_l0(theGrid,&theBV,theFirstBV,1)!=GM_OK) return (__LINE__);
          }
          else
          {
            if (CreateBlockvector_l0(theGrid,&theBV,theLastBV,0)!=GM_OK) return (__LINE__);
          }
          theLastBV = theBV;
          if (init && (line==0))
          {
            SET_ORD_GCL(BVNUMBER(theLastBV),GM_GEN_LAST,0,line);
          }
          else
          {
            SET_ORD_GCL(BVNUMBER(theLastBV),GM_GEN_LAST,cycle,line);
          }
          BVFIRSTVECTOR(theLastBV)        = LAST_line_end;                                      /* end because order is reverted below */

          line++;
          LAST_line_start = PREDVC(LAST_next_out);
          LAST_line_end   = LAST_last_in;
        }
      }
      LAST_next_out = LAST_last_in;

      /* check if there is a !used vec in the strong fifo */
      if ((theVector=FindOptimalStrong(&LFifo))!=NULL)
        LAST_last_in = CutStrongLine(theVector,LAST_last_in);
    }
    while (theVector!=NULL);

    init = 0;

    /****************************************************************************/
    /*	get CUT (or Feedback Vertex)-set and do what needs to be done			*/
    /****************************************************************************/

    FIRST_last_in = (*FindCutSet)(theGrid,FIRST_last_in,&nCUT);
    if (FIRST_last_in==NULL) nCUT = 0;
    else PREDVC(FIRST_last_in) = NULL;

    nCutTot += nCUT;

  } while (nCUT>0);

  /* TODO: use dlmgr macros for lists after streamwise ordering */

  /* insert FIRST list one-by-one to LASTVECTOR list of grid */
  LASTVECTOR(theGrid) = NULL;
  for (theVector=PREDVC(&FIRST_handle); theVector!=NULL; theVector=predVector)
  {
    predVector                      = PREDVC(theVector);
    PREDVC(theVector)   = LASTVECTOR(theGrid);
    LASTVECTOR(theGrid) = theVector;
  }

  /* insert LAST list as a whole to LASTVECTOR list of grid */
  PREDVC(LAST_last_in) = LASTVECTOR(theGrid);
  LASTVECTOR(theGrid)  = PREDVC(&LAST_handle);

  /* construct SUCCVC list */
  succVector = NULL;
  for (theVector=LASTVECTOR(theGrid); theVector!=NULL; theVector=PREDVC(theVector))
  {
    SUCCVC(theVector) = succVector;
    succVector = theVector;
  }
  SFIRSTVECTOR(theGrid) = succVector;
  PREDVC(succVector)   = NULL;

  /* set pointers in BLOCKVECTORs */
  BVENDVECTOR(GLASTBV(theGrid)) = NULL;
  for (theBV=GLASTBV(theGrid); theBV!=NULL; theBV=BVPRED(theBV))
  {
    if (BVSUCC(theBV)!=NULL && BVENDVECTOR(theBV)==NULL)
      BVENDVECTOR(theBV)=BVFIRSTVECTOR(BVSUCC(theBV));
    if (BVFIRSTVECTOR(theBV)==NULL)
      BVFIRSTVECTOR(theBV)=BVENDVECTOR(theBV);
  }

  /* set VCCUT */
  for (theBV=GFIRSTBV(theGrid); theBV!=NULL; theBV=BVSUCC(theBV))
  {
    i = 0;
    if (BV_GEN(theBV)==BV_GEN_C)
      for (theVector=BVFIRSTVECTOR(theBV); theVector!=BVENDVECTOR(theBV); theVector=SUCCVC(theVector))
      {
        SETVCCUT(theVector,1);
        i++;
      }
    else
      for (theVector=BVFIRSTVECTOR(theBV); theVector!=BVENDVECTOR(theBV); theVector=SUCCVC(theVector))
      {
        SETVCCUT(theVector,0);
        i++;
      }
    bvn = BVNUMBER(theBV);
    if (verboselevel>1)
      UserWriteF("# %d members in %c%d,%d\n",i,gen_label[ORD_GEN(bvn)],ORD_CYC(bvn),ORD_LIN(bvn));
  }

  if (verboselevel>0)
  {
    if (verboselevel>1)
      UserWrite("#\n# summary:\n");
    UserWriteF("# %d cycles: %d cut from %d\n",(int)cycle,(int)nCutTot,(int)NVEC(theGrid));
    a = POW((DOUBLE)NVEC(theGrid),(DOUBLE)(DIM-1)/(DOUBLE)DIM);     a = (DOUBLE)nCutTot/a;
    UserWriteF("# corr. to %6.2f hyp. planes\n",(float)a);
  }

  /* check # members of succ list */
  i=0;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector)) i++;
  if (NVEC(theGrid) != i)
  {
    UserWrite("vectorstructure corrupted\n");
    return (__LINE__);
  }

  /* check # members of pred list */
  i=0;
  for (theVector=LASTVECTOR(theGrid); theVector!= NULL; theVector=PREDVC(theVector)) i++;
  if (NVEC(theGrid) != i)
  {
    UserWrite("vectorstructure corrupted\n");
    return (__LINE__);
  }

  /* set index w.r.t. new order, beginning with 1 */
  i = 1;
  for (theVector=FIRSTVECTOR(theGrid); theVector!= NULL; theVector=SUCCVC(theVector))
    VINDEX(theVector) = i++;

  ReleaseTmpMem(MGHEAP(MYMG(theGrid)),MarkKey);

  return (0);
}

/****************************************************************************/
/*D
   LineOrderVectors	-  Driver for general vector ordering

   SYNOPSIS:
   INT LineOrderVectors (MULTIGRID *theMG, INT levels,
                        const char *dependency, const char *dep_options, const char *findcut);

   PARAMETERS:
   .  theMG -  multigrid to order
   .  levels -  GM_ALL_LEVELS or GM_CURRENT_LEVEL
   .  mode - GM_FCFCLL or GM_FFCCLL (see 'orderv' command)
   .  dependency - name of user defined dependency item
   .  dep_options - options for user dependency function

   DESCRIPTION:
   This function orders 'VECTOR's in a multigrid according to the dependency
   function provided.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR if error occured.
   D*/
/****************************************************************************/

INT LineOrderVectors (MULTIGRID *theMG, INT levels, const char *dependency, const char *dep_options, const char *findcut, INT verboselevel)
{
  INT i, currlevel, baselevel;
  ALG_DEP *theAlgDep;
  FIND_CUT *theFindCut;
  GRID *theGrid;
  DependencyProcPtr DependencyProc;
  INT err;

  /* current level */
  currlevel = theMG->currentLevel;

  /* get algebraic dependency */
  theAlgDep = (ALG_DEP *) SearchEnv(dependency,"/Alg Dep",theAlgDepVarID,theAlgDepDirID);
  if (theAlgDep==NULL)
  {
    UserWrite("algebraic dependency not found\n");
    return(GM_ERROR);
  }
  DependencyProc = theAlgDep->DependencyProc;
  if (DependencyProc==NULL)
  {
    UserWrite("don't be stupid: implement a dependency!\n");
    return(GM_ERROR);
  }

  /* get find cut dependency */
  if (findcut==NULL)
  {
    FindCutSet = FeedbackVertexVectors;
    UserWrite("default cut set proc:\n    leaving order of cyclic dependencies unchanged\n");
  }
  else
  {
    theFindCut = (FIND_CUT *) SearchEnv(findcut,"/FindCut",theFindCutVarID,theFindCutDirID);
    if (theFindCut==NULL)
    {
      UserWrite("find cut proc not found\n");
      return(GM_ERROR);
    }
    FindCutSet = theFindCut->FindCutProc;
    if (FindCutSet==NULL)
    {
      UserWrite("don't be stupid: implement a find cut proc!\n");
      return(GM_ERROR);
    }
  }
  if (AllocateControlEntry(VECTOR_CW,VCSTRONG_LEN,&ce_VCSTRONG) != GM_OK)
    return (GM_ERROR);

  /* go */
  if (levels==GM_ALL_LEVELS)
    baselevel = 0;
  else
    baselevel = currlevel;
  for (i=baselevel; i<=currlevel; i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if ((*DependencyProc)(theGrid,dep_options))
    {
      PrintErrorMessage('E',"LineOrderVectors","DependencyProc failed");
      return (GM_ERROR);
    }
    if ((err=LineOrderVectorsAlgebraic(theGrid,verboselevel))!=0)
    {
      PrintErrorMessage('E',"LineOrderVectors","LineOrderVectorsAlgebraic failed");
      return (GM_ERROR);
    }
  }

  FreeControlEntry(ce_VCSTRONG);

  return (GM_OK);
}

/****************************************************************************/
/*D
   PrepareForLineorderVectors - prepare for lineorderv

   SYNOPSIS:
   INT PrepareForLineorderVectors (GRID *theGrid)

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function resets flags in vectors for the line ordering algorithm.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT PrepareForLineorderVectors (GRID *theGrid)
{
  VECTOR *vec;

  for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=SUCCVC(vec))
  {
    SETVCUSED(vec,FALSE);
    SETVCFLAG(vec,FALSE);
  }

  return (0);
}

/****************************************************************************/
/*D
   MarkBeginEndForLineorderVectors - flag begin or end vectors for lineorderv

   SYNOPSIS:
   INT MarkBeginEndForLineorderVectors (ELEMENT *elem, INT dt, INT ot, const INT *mark)

   PARAMETERS:
   .  elem - set flag in marked vectors of the element
   .  dt - data types of vectors needed
   .  ot - object types of vectors needed
   .  mark - set GM_LOV_BEGIN or GM_LOV_END or no flag

   DESCRIPTION:
   This function sets flags in vectors for the line ordering algorithm.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT MarkBeginEndForLineorderVectors (ELEMENT *elem, INT dt, INT ot, const INT *mark)
{
  VECTOR *vList[MAX_ELEM_VECTORS];
  INT i,cnt;

  if (GetVectorsOfDataTypesInObjects(elem,dt,ot,&cnt,vList)!=GM_OK)
    REP_ERR_RETURN (1);

  for (i=0; i<cnt; i++)
    switch (mark[i])
    {
    case GM_LOV_BEGIN :      SETVCUSED(vList[i],TRUE);
    case GM_LOV_END :        SETVCFLAG(vList[i],TRUE);
    }

  return (0);
}

/****************************************************************************/
/*D
   RevertVecOrder - revert order in vector list

   SYNOPSIS:
   INT RevertVecOrder (GRID *theGrid);

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   This function revertes the order in the vector list.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT RevertVecOrder (GRID *theGrid)
{
  VECTOR *vec,*tmp;
  BLOCKVECTOR     *bv;

  for (vec=FIRSTVECTOR(theGrid); vec!=NULL; vec=PREDVC(vec))
  {
    SWAP(PREDVC(vec),SUCCVC(vec),tmp);
  }
  SWAP(SFIRSTVECTOR(theGrid),LASTVECTOR(theGrid),tmp);

  /* also change the blockvectors */
  for (bv=GFIRSTBV(theGrid); bv!=NULL; bv=BVSUCC(bv))
  {
    tmp = BVFIRSTVECTOR(bv);
    BVFIRSTVECTOR(bv) = (BVENDVECTOR(bv)==NULL) ? FIRSTVECTOR(theGrid) : SUCCVC(BVENDVECTOR(bv));
    BVENDVECTOR(bv)   = SUCCVC(tmp);
  }

  return (0);
}

/****************************************************************************/
/*D
   LexAlgDep - Dependency function for lexicographic ordering

   SYNOPSIS:
   static INT LexAlgDep (GRID *theGrid, char *data);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  data - option string from 'orderv' command

   DESCRIPTION:
   This function defines a dependency function for lexicographic ordering.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT LexAlgDep (GRID *theGrid, const char *data)
{
  MULTIGRID *theMG;
  VECTOR *theVector,*NBVector;
  MATRIX *theMatrix;
  DOUBLE_VECTOR pos,nbpos;
  DOUBLE diff[DIM];
  INT i,order,res;
  INT Sign[DIM],Order[DIM],xused,yused,zused,error,SpecialTreatSkipVecs;
  char ord[3];
  BVP *theBVP;
  BVP_DESC *theBVPDesc;

  /* read ordering directions */
        #ifdef __TWODIM__
  res = sscanf(data,expandfmt("%2[rlud]"),ord);
        #else
  res = sscanf(data,expandfmt("%3[rlbfud]"),ord);
        #endif
  if (res!=1)
  {
    PrintErrorMessage('E',"LexAlgDep","could not read order type");
    return(1);
  }
  if (strlen(ord)!=DIM)
  {
                #ifdef __TWODIM__
    PrintErrorMessage('E',"LexAlgDep","specify 2 chars out of 'rlud'");
                #else
    PrintErrorMessage('E',"LexAlgDep","specify 3 chars out of 'rlbfud'");
                #endif
    return(1);
  }
  error = xused = yused = zused = FALSE;
  for (i=0; i<DIM; i++)
    switch (ord[i])
    {
    case 'r' :
      if (xused) error = TRUE;
      xused = TRUE;
      Order[i] = _X_; Sign[i] =  1; break;
    case 'l' :
      if (xused) error = TRUE;
      xused = TRUE;
      Order[i] = _X_; Sign[i] = -1; break;

                        #ifdef __TWODIM__
    case 'u' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] =  1; break;
    case 'd' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] = -1; break;
                        #else
    case 'b' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] =  1; break;
    case 'f' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] = -1; break;

    case 'u' :
      if (zused) error = TRUE;
      zused = TRUE;
      Order[i] = _Z_; Sign[i] =  1; break;
    case 'd' :
      if (zused) error = TRUE;
      zused = TRUE;
      Order[i] = _Z_; Sign[i] = -1; break;
                        #endif
    }
  if (error)
  {
    PrintErrorMessage('E',"LexAlgDep","bad combination of 'rludr' or 'rlbfud' resp.");
    return(1);
  }

  /* treat vectors with skipflag set specially? */
  SpecialTreatSkipVecs = FALSE;
  if              (strchr(data,'<')!=NULL)
    SpecialTreatSkipVecs = GM_PUT_AT_BEGIN;
  else if (strchr(data,'>')!=NULL)
    SpecialTreatSkipVecs = GM_PUT_AT_END;

  theMG   = MYMG(theGrid);

  /* find an approximate measure for the mesh size */
  theBVP = MG_BVP(theMG);
  theBVPDesc = MG_BVPD(theMG);
  InvMeshSize = POW2(GLEVEL(theGrid)) * pow(NN(GRID_ON_LEVEL(theMG,0)),1.0/DIM) / BVPD_RADIUS(theBVPDesc);

  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    VectorPosition(theVector,pos);

    for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      NBVector = MDEST(theMatrix);

      SETMUP(theMatrix,0);
      SETMDOWN(theMatrix,0);

      if (SpecialTreatSkipVecs)
      {
        if (VECSKIP(theVector) && !VECSKIP(NBVector))
          if (SpecialTreatSkipVecs==GM_PUT_AT_BEGIN)
            order = -1;
          else
            order =  1;

        if (VECSKIP(NBVector) && !VECSKIP(theVector))
          if (SpecialTreatSkipVecs==GM_PUT_AT_BEGIN)
            order =  1;
          else
            order = -1;
      }
      else
      {
        VectorPosition(NBVector,nbpos);

        V_DIM_SUBTRACT(nbpos,pos,diff);
        V_DIM_SCALE(InvMeshSize,diff);

        if (fabs(diff[Order[DIM-1]])<ORDERRES)
        {
                                        #ifdef __THREEDIM__
          if (fabs(diff[Order[DIM-2]])<ORDERRES)
          {
            if (diff[Order[DIM-3]]>0.0) order = -Sign[DIM-3];
            else order =  Sign[DIM-3];
          }
          else
                                        #endif
          if (diff[Order[DIM-2]]>0.0) order = -Sign[DIM-2];
          else order =  Sign[DIM-2];
        }
        else
        {
          if (diff[Order[DIM-1]]>0.0) order = -Sign[DIM-1];
          else order =  Sign[DIM-1];
        }
      }
      if (order==1) SETMUP(theMatrix,1);
      else SETMDOWN(theMatrix,1);
    }
  }

  return (0);
}

/****************************************************************************/
/*D
   StrongLexAlgDep - Dependency function for lexicographic ordering

   SYNOPSIS:
   static INT StrongLexAlgDep (GRID *theGrid, char *data);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  data - option string from 'orderv' command

   DESCRIPTION:
   This function defines a dependency function for lexicographic ordering.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT StrongLexAlgDep (GRID *theGrid, const char *data)
{
  MULTIGRID *theMG;
  VECTOR *theVector,*NBVector;
  MATRIX *theMatrix;
  DOUBLE_VECTOR pos,nbpos;
  DOUBLE diff[DIM];
  INT i,order,res;
  INT Sign[DIM],Order[DIM],xused,yused,zused,error,SpecialTreatSkipVecs;
  char ord[3];
  BVP *theBVP;
  BVP_DESC *theBVPDesc;

  /* read ordering directions */
        #ifdef __TWODIM__
  res = sscanf(data,expandfmt("%2[rlud]"),ord);
        #else
  res = sscanf(data,expandfmt("%3[rlbfud]"),ord);
        #endif
  if (res!=1)
  {
    PrintErrorMessage('E',"LexAlgDep","could not read order type");
    return(1);
  }
  if (strlen(ord)!=DIM)
  {
                #ifdef __TWODIM__
    PrintErrorMessage('E',"LexAlgDep","specify 2 chars out of 'rlud'");
                #else
    PrintErrorMessage('E',"LexAlgDep","specify 3 chars out of 'rlbfud'");
                #endif
    return(1);
  }
  error = xused = yused = zused = FALSE;
  for (i=0; i<DIM; i++)
    switch (ord[i])
    {
    case 'r' :
      if (xused) error = TRUE;
      xused = TRUE;
      Order[i] = _X_; Sign[i] =  1; break;
    case 'l' :
      if (xused) error = TRUE;
      xused = TRUE;
      Order[i] = _X_; Sign[i] = -1; break;

                        #ifdef __TWODIM__
    case 'u' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] =  1; break;
    case 'd' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] = -1; break;
                        #else
    case 'b' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] =  1; break;
    case 'f' :
      if (yused) error = TRUE;
      yused = TRUE;
      Order[i] = _Y_; Sign[i] = -1; break;

    case 'u' :
      if (zused) error = TRUE;
      zused = TRUE;
      Order[i] = _Z_; Sign[i] =  1; break;
    case 'd' :
      if (zused) error = TRUE;
      zused = TRUE;
      Order[i] = _Z_; Sign[i] = -1; break;
                        #endif
    }
  if (error)
  {
    PrintErrorMessage('E',"LexAlgDep","bad combination of 'rludr' or 'rlbfud' resp.");
    return(1);
  }

  /* treat vectors with skipflag set specially? */
  SpecialTreatSkipVecs = FALSE;
  if              (strchr(data,'<')!=NULL)
    SpecialTreatSkipVecs = GM_PUT_AT_BEGIN;
  else if (strchr(data,'>')!=NULL)
    SpecialTreatSkipVecs = GM_PUT_AT_END;

  theMG   = MYMG(theGrid);

  /* find an approximate measure for the mesh size */
  theBVP = MG_BVP(theMG);
  theBVPDesc = MG_BVPD(theMG);
  InvMeshSize = POW2(GLEVEL(theGrid)) * pow(NN(GRID_ON_LEVEL(theMG,0)),1.0/DIM) / BVPD_RADIUS(theBVPDesc);

  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    VectorPosition(theVector,pos);

    for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
    {
      NBVector = MDEST(theMatrix);

      SETMUP(theMatrix,0);
      SETMDOWN(theMatrix,0);
      SETMUSED(theMatrix,0);

      VectorPosition(NBVector,nbpos);

      V_DIM_SUBTRACT(nbpos,pos,diff);
      V_DIM_SCALE(InvMeshSize,diff);

      if (fabs(diff[Order[DIM-1]])<ORDERRES)
      {
                                #ifdef __THREEDIM__
        if (fabs(diff[Order[DIM-2]])<ORDERRES)
        {
          if (diff[Order[DIM-3]]>0.0) order = -Sign[DIM-3];
          else order =  Sign[DIM-3];
        }
        else
                                #endif
        if (diff[Order[DIM-2]]>0.0) order = -Sign[DIM-2];
        else order =  Sign[DIM-2];
        SETMUSED(theMatrix,1);
      }
      else
      {
        if (diff[Order[DIM-1]]>0.0) order = -Sign[DIM-1];
        else order =  Sign[DIM-1];
      }
      switch (order)
      {
      case  1 : SETMUP(theMatrix,1); break;
      case -1 : SETMDOWN(theMatrix,1); break;
      case  0 : SETMUP(theMatrix,1); SETMDOWN(theMatrix,1); break;
      }
    }
  }
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    SETVCUSED(theVector,0);
    SETVCFLAG(theVector,0);

    for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
      if (MUP(theMatrix) && !MUSED(theMatrix))
        break;
    if (theMatrix==NULL)
      SETVCUSED(theVector,1);
  }
  for (theVector=FIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
    for (theMatrix=MNEXT(VSTART(theVector)); theMatrix!=NULL; theMatrix=MNEXT(theMatrix))
      if (MUSED(theMatrix) && MUSED(MADJ(theMatrix)))
      {
        SETMUP(theMatrix,1);
        SETMDOWN(theMatrix,1);
      }

  return (0);
}

#ifdef __BLOCK_VECTOR_DESC__
/****************************************************************************/
/*D
   SetLevelnumberBV - set the level-field in a blockvector-tree

   SYNOPSIS:
   void SetLevelnumberBV( BLOCKVECTOR *bv, INT level );

   PARAMETERS:
   .  bv - first of the blockvectors whose levelnumber should be set
   .  level - levelnumber of this blockvector

   DESCRIPTION:
   The blockvectors form a hierarchical structure; thus each blockvector can
   be assigned a 'level'. This function traverses the blockvector-tree
   starting at the given 'bv' and sets in each blockvector the level field
   starting with the given 'level'. In this way it is possible to number also
   parts of the global blockvector-tree (to skip for example meaningless
   administrative, covering blocks).

   SEE ALSO:
   BLOCKVECTOR

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void SetLevelnumberBV( BLOCKVECTOR *bv, INT level )
{
  ASSERT( level < (1<<BVLEVEL_LEN) );           /* increase BVLEVEL_LEN in gm.h */

  if ( bv == NULL )
    return;

  if ( level >= 0 )
    SETBVLEVEL(bv,level);
  else
    SETBVLEVEL(bv,0);

  if ( !BV_IS_LEAF_BV( bv ) )
    for ( bv = BVDOWNBV(bv); bv != NULL; bv = BVSUCC(bv) )
      SetLevelnumberBV( bv, level + 1 );
}


/****************************************************************************/
/* (without *D, since only internal function)

   CreateBVPlane - Creates a stripewise decomposition of a 2D rectangular domain with a regular mesh

   SYNOPSIS:
        static INT CreateBVPlane( BLOCKVECTOR **bv_plane, const BV_DESC *bvd_plane, const BV_DESC_FORMAT *bvdf, VECTOR **v, INT stripes, INT vectors_per_stripe, GRID *grid );

   PARAMETERS:
   .  bv_plane - handle for the blockvector covering the plane
   .  bvd_plane - blockvector description for the block 'bv_plane'
   .  bvdf - blockvector description format for 'bvd_plane'
   .  v - input: handle to the first vector; output: handle to the 'last' vector in this plane
   .  stripes  - number of stripes to be constructed
   .  vectors_per_stripe - number of vectors a stripe should contain
   .  grid - the grid containing the vectors to be structured

   DESCRIPTION:
   From the list of vectors the blockvector-list is generated in
   the following way:
   First a blockvector 'bv_plane' corresponding to 'bvd_plane' is constructed
   which will cover all the blockvectors for the single stripes and
   thus all the vectors of the plane. Beginning
   at the start of the vector-list ('v'), blocks are constructed containing
   'vectors_per_stripe' consecutive vectors. The start of the new
   blockvector-list is anchored in the blockvector 'bv_plane'.
   This blockvectors are
   numbered beginning from 0 along the construction. If the vector-list
   is not long enough, the function will not be troubled but construct
   as far as possible the structure and GM_INCONSISTANCY is returned.

   REQUIREMENTS:
   The macro '__BLOCK_VECTOR_DESC__' must be defined in 'gm.h' to enable
   blockvector descriptions.

   REMARK:
   This is only an auxiliary function. In the returned 'bv_plane' the
   components BVSUCC, BVPRED, BVNUMBER are not yet set.

   APPLICATION:
   If the vector-list is ordered lexicographic, each line resp. column
   containing 'vectors_per_stripe' vectors, then this function constructs
   the blockvector structure belonging to a
   `linewise domain decomposition`.

   If 'vectors_per_stripe' is a natural multiple of the actual linewidth
   (resp. columnlength), this function constructs the blockvector
   structure belonging to a `stripewise domain decomposition`, each stripe
   consisting of 'vectors_per_stripe'/linewidth lines (resp. columns).

   RETURN VALUE:

   .n GM_OK if ok
   .n GM_OUT_OF_MEM if there was not enough memory to allocate all blockvectors
   .n GM_INCONSISTANCY if the vector-list was too short

   SEE ALSO:
   BLOCKVECTOR, CreateBVStripe2D, CreateBVStripe3D, CreateBVDomainHalfening

   D*/
/****************************************************************************/

static INT CreateBVPlane( BLOCKVECTOR **bv_plane, const BV_DESC *bvd_plane, const BV_DESC_FORMAT *bvdf, VECTOR **v, INT stripes, INT vectors_per_stripe, GRID *grid )
{
  BLOCKVECTOR *bv;
  register BLOCKVECTOR *prev;
  register INT i, j;
  register VECTOR *pred_v;

  (void)CreateBlockvector( grid, bv_plane );
  if ( *bv_plane == NULL )
    return GM_OUT_OF_MEM;

  SETBVDOWNTYPE( *bv_plane, BVDOWNTYPEBV );             /* further blocks follow */
  SETBVTVTYPE( *bv_plane, BV1DTV );
  SETBVORIENTATION( *bv_plane, BVNOORIENTATION );
  BVFIRSTVECTOR( *bv_plane ) = *v;

  /* construct for each stripe a blockvector */
  for ( i = 0; (i < stripes) && (*v != NULL); i++ )
  {
    (void)CreateBlockvector( grid, &bv );

    if ( i == 0 )
    {
      BVDOWNBV( *bv_plane ) = bv;                   /* anchor stripe-blocks in the block of the plane */
      if ( bv == NULL )
        return GM_OUT_OF_MEM;
      BVPRED( bv ) = NULL;
    }
    else
    {
      BVSUCC( prev ) = bv;                              /* continue blockvector list */
      if ( bv == NULL )
        return GM_OUT_OF_MEM;
      BVPRED( bv ) = prev;
    }

    prev = bv;                                                          /* prepare for the next loop */

    SETBVDOWNTYPE( bv, BVDOWNTYPEVECTOR );              /* block is at last block level */
    SETBVTVTYPE( bv, BV1DTV );
    SETBVORIENTATION( bv, BVHORIZONTAL );
    BVNUMBER( bv ) = i;

    /* let the blockvector point to the vector */
    BVDOWNVECTOR( bv ) = *v;
    BVFIRSTVECTOR( bv ) = *v;

    /* find successor of the last vector of this blockvector and store it;
       update the blockvector description for all vectors of this
       blockvector
     */
    for ( j = vectors_per_stripe; (j > 0) && (*v != NULL); j-- )
    {
      VBVD( *v ) = *bvd_plane;
      BVD_PUSH_ENTRY( &VBVD( *v ), i, bvdf );
      pred_v = *v;
      *v = SUCCVC( *v );
    }
    BVLASTVECTOR( bv ) = pred_v;
    BVNUMBEROFVECTORS( bv ) = vectors_per_stripe - j;
  }
  BVSUCC( bv ) = NULL;                          /* end of the blockvector list */

  *v = pred_v;
  BVLASTVECTOR( *bv_plane ) = *v;
  BVDOWNBVLAST( *bv_plane ) = bv;
  BVNUMBEROFVECTORS( *bv_plane ) = stripes * vectors_per_stripe - j;

  /* the for loop exited premature because of v == NULL */
  if ( (i < stripes) || (j != 0) )
    return GM_INCONSISTANCY;

  return GM_OK;
}
#endif /* __BLOCK_VECTOR_DESC__ */


/****************************************************************************/
/*D
   CreateBVStripe2D - Creates a stripewise domain decomposition of a 2D rectangular domain with a regular mesh

   SYNOPSIS:
   INT CreateBVStripe2D( GRID *grid, INT vectors, INT vectors_per_stripe );

   PARAMETERS:
   .  grid - the grid containing the vectors to be structured
   .  vectors  - number of vectors, i.e. gridpoints, in the inner
   .  vectors_per_stripe - number of vectors a stripe should contain

   DESCRIPTION:
   From the list of vectors the blockvector-list is generated in
   the following way:
   The first blockvector level consists of 2 blockvectors: the first
   covering all inner vectors (number 0), the second all boundary vectors
   (number 1). Beginning
   at the start of the vector-list, blocks are constructed containing
   'vectors_per_stripe' consecutive vectors. The start of the new
   blockvector-list is anchored in the blockvector for the inner vectors.
   This blockvectors are
   numbered beginning from 0 along the construction. If the vector-list
   is shorter than 'vectors', the function will not be troubled
   and GM_INCONSISTANCY is returned. If the vector-list
   is longer than 'vectors', these vectors are covered in a seperate
   blockvector with number 1 (this overlapping vectors might be for
   example dirichlet boundary vectors).

   For the blockvectors in the tree with root 'all inner' the level-field is
   set beginning with 0.

   REQUIREMENTS:
   The macro '__BLOCK_VECTOR_DESC__' must be defined in 'gm.h' to enable
   blockvector descriptions.

   APPLICATION:
   If the vector-list is ordered lexicographic, each line resp. column
   containing 'vectors_per_stripe' vectors, then this function constructs
   the blockvector structure belonging to a
   `linewise domain decomposition`.

   If 'vectors_per_stripe' is a natural multiple of the actual linewidth
   (resp. columnlength), this function constructs the blockvector
   structure belonging to a `stripewise domain decomposition`, each stripe
   consisting of 'vectors_per_stripe'/linewidth lines (resp. columns).

   RETURN VALUE:

   .n GM_OK if ok
   .n GM_OUT_OF_MEM if there was not enough memory to allocate all blockvectors
   .n GM_INCONSISTANCY if the vector-list was too short

   SEE ALSO:
   BLOCKVECTOR, CreateBVStripe3D, CreateBVDomainHalfening

   D*/
/****************************************************************************/

INT CreateBVStripe2D( GRID *grid, INT vectors, INT vectors_per_stripe )
{
#ifdef __BLOCK_VECTOR_DESC__
  BLOCKVECTOR *bv_inner, *bv_boundary;
  VECTOR *v;
  INT nr_blocks, ret;
  BV_DESC bvd;

  if ( GFIRSTBV( grid ) != NULL )
    FreeAllBV( grid );

  /* number of blockvectors to be constructed */
  nr_blocks = ( vectors + vectors_per_stripe - 1) / vectors_per_stripe;

  BVD_INIT( &bvd );
  BVD_PUSH_ENTRY( &bvd, 0, &two_level_bvdf );
  v = FIRSTVECTOR( grid );
  ret = CreateBVPlane( &bv_inner, &bvd, &two_level_bvdf, &v, nr_blocks, vectors_per_stripe, grid );
  if ( ret != GM_OK )
  {
    FreeBVList( grid, bv_inner );
    return ret;
  }
  v = SUCCVC( v );

  (void)CreateBlockvector( grid, &bv_boundary );
  if ( bv_boundary == NULL )
  {
    FreeBVList( grid, bv_inner );
    return GM_OUT_OF_MEM;
  }

  GFIRSTBV( grid ) = bv_inner;          /* anchor blockvector list in the grid */
  GLASTBV( grid ) = bv_boundary;

  BVPRED( bv_inner ) = NULL;
  BVSUCC( bv_inner ) = bv_boundary;
  BVNUMBER( bv_inner ) = 0;

  BVPRED( bv_boundary ) = bv_inner;
  BVSUCC( bv_boundary ) = NULL;
  BVNUMBER( bv_boundary ) = 1;
  SETBVDOWNTYPE( bv_boundary, BVDOWNTYPEVECTOR );       /* block is at last block level */
  SETBVTVTYPE( bv_boundary, BV1DTV );
  SETBVORIENTATION( bv_boundary, BVNOORIENTATION );
  BVNUMBEROFVECTORS( bv_boundary ) = NVEC( grid ) - BVNUMBEROFVECTORS( bv_inner );
  BVFIRSTVECTOR( bv_boundary ) = v;
  BVLASTVECTOR( bv_boundary ) = LASTVECTOR( grid );

  /* complete the boundary vector list */
  BVD_INC_LAST_ENTRY( &bvd, 1, &two_level_bvdf );
  for ( ; v != NULL; v = SUCCVC( v ) )
    VBVD( v ) = bvd;

  SetLevelnumberBV( bv_inner, 0 );
  return GM_OK;
#else
  return(1);
#endif /* __BLOCK_VECTOR_DESC__ */
}


/****************************************************************************/
/*D
   CreateBVStripe3D - Creates a plane- and stripewise domain decomposition of a 3D rectangular domain with a regular mesh

   SYNOPSIS:
   INT CreateBVStripe3D( GRID *grid, INT inner_vectors, INT stripes_per_plane, INT vectors_per_stripe );

   PARAMETERS:
   .  grid - the grid containing the vectors to be structured
   .  inner_vectors  - number of vectors, i.e. gridpoints, in the inner
   .  stripes_per_plane - number of stripes forming a plane
   .  vectors_per_stripe - number of vectors a stripe should contain

   DESCRIPTION:
   From the list of vectors the blockvector-list is generated in
   the following way:
   The first blockvector level consists of 2 blockvectors: the first
   covering all inner vectors (number 0), the second all other vectors
   (number 1). Beginning
   at the start of the vector-list, plane-blockvectors are constructed: each
   plane-blockvector consists of a blockvector list of 'stripes_per_plane'
   stripe-blockvectors and each
   stripe-blockvector contains 'vectors_per_stripe' vectors. Plane-blocks are
   constructed as long as not 'inner_vectors' vectors are processed
   or all vectors in the grid are processed. The start of the new
   plane-blockvector-list is anchored in the blockvector for the
   inner vectors. The blockvectors in each list are
   numbered beginning from 0 along the construction. If the vector-list
   is shorter than 'inner_vectors', the function will not be troubled
   and GM_INCONSISTANCY is returned. If the vector-list
   is longer than 'vectors', these vectors are covered in a seperate
   blockvector with number 1 (this overlapping vectors might be for
   example dirichlet boundary vectors).

   REQUIREMENTS:
   The macro '__BLOCK_VECTOR_DESC__' must be defined in 'gm.h' to enable
   blockvector descriptions.

   APPLICATION:
   If the vector-list is ordered lexicographic and each line contains
   'vectors_per_stripe' vectors and each plane contains 'stripes_per_plane'
   lines, then this function constructs the blockvector structure belonging
   to a `plane- and linewise domain decomposition`.

   RETURN VALUE:

   .n GM_OK if ok
   .n GM_OUT_OF_MEM if there was not enough memory to allocate all blockvectors
   .n GM_INCONSISTANCY if the vector-list was too short

   SEE ALSO:
   BLOCKVECTOR, CreateBVStripe2D, CreateBVDomainHalfening

   D*/
/****************************************************************************/

INT CreateBVStripe3D( GRID *grid, INT inner_vectors, INT stripes_per_plane, INT vectors_per_stripe )
{
#ifdef __BLOCK_VECTOR_DESC__
  BLOCKVECTOR *bv_inner, *bv_boundary, *bv_plane, *prev;
  VECTOR *v;
  INT i, nr_planes, nr_vectors, ret;
  BV_DESC bvd;

  if ( GFIRSTBV( grid ) != NULL )
    FreeAllBV( grid );

  /* number of planes to be constructed */
  nr_planes = ( inner_vectors + stripes_per_plane*vectors_per_stripe - 1) / (stripes_per_plane*vectors_per_stripe);

  v = FIRSTVECTOR( grid );

  (void)CreateBlockvector( grid, &bv_inner );
  if ( bv_inner == NULL )
    return GM_OUT_OF_MEM;

  (void)CreateBlockvector( grid, &bv_boundary );
  if ( bv_boundary == NULL )
  {
    DisposeBlockvector( grid, bv_inner );
    return GM_OUT_OF_MEM;
  }

  GFIRSTBV( grid ) = bv_inner;          /* anchor blockvector list in the grid */
  GLASTBV( grid ) = bv_boundary;

  BVPRED( bv_inner ) = NULL;
  BVSUCC( bv_inner ) = bv_boundary;
  BVNUMBER( bv_inner ) = 0;
  SETBVDOWNTYPE( bv_inner, BVDOWNTYPEBV );              /* block is parent for all planes */
  SETBVTVTYPE( bv_inner, BV1DTV );
  SETBVORIENTATION( bv_inner, BVNOORIENTATION );
  BVFIRSTVECTOR( bv_inner ) = v;

  BVPRED( bv_boundary ) = bv_inner;
  BVSUCC( bv_boundary ) = NULL;
  BVNUMBER( bv_boundary ) = 1;
  SETBVDOWNTYPE( bv_boundary, BVDOWNTYPEVECTOR );       /* block is at last block level */
  SETBVTVTYPE( bv_boundary, BV1DTV );
  SETBVORIENTATION( bv_boundary, BVNOORIENTATION );

  /* construct for each plane a blockvector-tree */
  BVD_INIT( &bvd );
  BVD_PUSH_ENTRY( &bvd, 0, &three_level_bvdf );         /* inner vectors */
  BVD_PUSH_ENTRY( &bvd, 0, &three_level_bvdf );         /* 1. plane */
  nr_vectors = 0;
  for ( i = 0; (i < nr_planes) && (v != NULL); i++ )
  {
    if ( i != 0 )
      v = SUCCVC( v );                   /* set v from the last vector of the previous block to the first of the current block */

    ret = CreateBVPlane( &bv_plane, &bvd, &three_level_bvdf, &v, stripes_per_plane, vectors_per_stripe, grid );
    SETBVTVTYPE( bv_plane, BV2DTV );                    /* use tensor product testvector (full 2D testvector) */
    if ( ret == GM_OUT_OF_MEM )
    {
      FreeAllBV( grid );
      return ret;
    }

    if ( i == 0 )
    {
      BVDOWNBV( bv_inner ) = bv_plane;                   /* anchor stripe-blocks in the block of the plane */
      BVPRED( bv_plane ) = NULL;
    }
    else
    {
      BVSUCC( prev ) = bv_plane;                                /* continue blockvector list */
      BVPRED( bv_plane ) = prev;
    }

    prev = bv_plane;                                                    /* prepare for the next loop */

    BVNUMBER( bv_plane ) = i;
    nr_vectors += BVNUMBEROFVECTORS( bv_plane );

    BVD_INC_LAST_ENTRY( &bvd, 1, &three_level_bvdf );
  }
  BVSUCC( bv_plane ) = NULL;                            /* end of the blockvector list */

  BVLASTVECTOR( bv_inner ) = v;
  BVDOWNBVLAST( bv_inner ) = bv_plane;
  BVNUMBEROFVECTORS( bv_inner ) = nr_vectors;

  v = SUCCVC( v );
  BVFIRSTVECTOR( bv_boundary ) = v;
  BVLASTVECTOR( bv_boundary ) = LASTVECTOR( grid );
  BVNUMBEROFVECTORS( bv_boundary ) = NVEC( grid ) - BVNUMBEROFVECTORS( bv_inner );

  /* complete the boundary vector list */
  BVD_INIT( &bvd );
  BVD_PUSH_ENTRY( &bvd, 1, &three_level_bvdf );         /* boundary vectors */
  for ( ; v != NULL; v = SUCCVC( v ) )
    VBVD( v ) = bvd;

  if ( ret != GM_OK )
    return ret;

  if ( BVNUMBEROFVECTORS(bv_inner) != nr_planes*stripes_per_plane*vectors_per_stripe )
    return GM_INCONSISTANCY;

  SetLevelnumberBV( bv_inner, 0 );

  return GM_OK;

#else
  return(1);
#endif /* __BLOCK_VECTOR_DESC__ */
}


/****************************************************************************/
/*D
   CreateBVDomainHalfening - Creates a recursive domain halfening decomposition

   SYNOPSIS:
   INT CreateBVDomainHalfening( GRID *grid, INT side )

   PARAMETERS:
   .  grid - the grid containing the vectors to be structured
   .  side - number of points on the side of the quadratic mesh

   DESCRIPTION:
   The grid must be similar to a quadratic mesh. The vectors must be numbered
   linewise lexicographic. Then this function constructs a hierarchy
   of blockvectors describing a recursive domain halfening.

   A quadratic domain is halfened by cutting along a vertical line. The
   points on this line are gathered into subdomain with number 2. The points
   left to this interface-line are collected in subdomain 0, the right ones
   into subdomain 1. The subdomains are listed according to their numbers.
   Each of the 2 resulting subdomains 0 and 1 are now halfened horizontally
   and processed analogical. The halfening is recursively proceeded until
   a subdomain consists of less than 10 points. The initial grid should have
   a side length allowing all halfenings without a remainder.

   The first blockvector in the hierarchy is covering all the vectors.

   NOTE:
   If the grid contains Dirichlet vectors, they must be hidden from this
   routine; e.g. order the whole vector list such that the Dirichlet vectors
   are at the end of the list; then let the LASTVECTOR of the grid be the
   vector in front of the first Dirichlet vector and correct this temporarly
   setting after returning from this routine.

   RETURN VALUE:
   INT
   .n GM_OK if ok
   .n GM_OUT_OF_MEM if there is not enough memory to allocate the blockvectors

   SEE ALSO:
   BLOCKVECTOR, CreateBVStripe

   D*/
/****************************************************************************/

#ifndef ModelP
INT CreateBVDomainHalfening( GRID *grid, INT side, INT leaf_size )
{
#ifdef __BLOCK_VECTOR_DESC__
  BLOCKVECTOR *bv;
  INT ret;
  register VECTOR *v, *end_v;

  if ( GFIRSTBV( grid ) != NULL )
    FreeAllBV( grid );

  /* create first block of the hierarchy */
  if ( CreateBlockvector( grid, &bv ) != GM_OK )
    return GM_OUT_OF_MEM;

  GFIRSTBV( grid ) = bv;
  GLASTBV( grid ) = bv;
  SETBVDOWNTYPE( bv, BVDOWNTYPEVECTOR );
  SETBVTVTYPE( bv, BV1DTV );
  BVDOWNVECTOR( bv ) = FIRSTVECTOR( grid );
  BVPRED( bv ) = NULL;
  BVSUCC( bv ) = NULL;
  BVNUMBER( bv ) = 0;
  BVLASTVECTOR( bv ) = LASTVECTOR( grid );
  BVFIRSTVECTOR( bv ) = FIRSTVECTOR( grid );

  /* set this global blockvector as the first stage of the hierarchy */
  end_v = SUCCVC(LASTVECTOR( grid ));
    #ifdef __BLOCK_VECTOR_DESC__
  for ( v = FIRSTVECTOR( grid ); v != end_v; v = SUCCVC( v ) )
    BVD_PUSH_ENTRY( &VBVD( v ), 0, &DH_bvdf );
    #endif

  if ( (ret = BlockHalfening( grid, bv, 0, 0, side, side, side, BV_VERTICAL, leaf_size )) != GM_OK )
  {
    FreeAllBV( grid );
    return ret;
  }

  /* set the number of covered VECTORS */
  BVNUMBEROFVECTORS( bv ) = BVNUMBEROFVECTORS( BVDOWNBV(bv) ) +
                            BVNUMBEROFVECTORS( BVSUCC(BVDOWNBV(bv)) ) +
                            BVNUMBEROFVECTORS( BVDOWNBVLAST(bv) );

  /* set FIRST- and LASTVECTOR */
  FIRSTVECTOR( grid ) = BVFIRSTVECTOR( bv );
  LASTVECTOR( grid ) = BVLASTVECTOR( bv );

  return GM_OK;
#else
  return(1);
#endif /* __BLOCK_VECTOR_DESC__ */
}
#endif /* ModelP */

#ifdef __BLOCK_VECTOR_DESC__
/****************************************************************************/
/* (without *D, since only internal function)

   BlockHalfening - internal function to perform the work of CreateBVDomainHalfening

   SYNOPSIS:
   static INT BlockHalfening( GRID *grid, BLOCKVECTOR *bv, INT left, INT bottom, INT width, INT height, INT side, INT orientation )

   PARAMETERS:
   .  grid - the grid containing the vectors to be structured
   .  bv - blockvector to be subdivided
   .  left - left side of the rectangle representing the 'bv' in the virtual mesh
   .  bottom - bottom side of the rectangle representing the 'bv' in the virtual mesh
   .  width -  sidewidth of the rectangle representing the 'bv' in the virtual mesh
   .  height - sideheight of the rectangle representing the 'bv' in the virtual mesh
   .  side - sidelength of the original square
   .  orientation - gives the orientation of the cutting line (horizontal or vertical)

   DESCRIPTION:
   Performs the work described for the function CreateBVDomainHalfening
   above.

   Consider a virtual equidistant mesh with meshwidth 1 in both directions.
   Together with the linewise lexicographic numbering of the vectors you can
   easily calculate the geometric position of each vector according to his
   number in this mesh. In this way it is easy to determine to which
   subdomain a given vector belongs. The "geometric" information to this
   function are related to this virtual mesh.

   RETURN VALUE:
   INT
   .n GM_OK if ok
   .n GM_OUT_OF_MEM if there is not enough memory to allocate the blockvectors

 */
/****************************************************************************/

static INT BlockHalfening( GRID *grid, BLOCKVECTOR *bv, INT left, INT bottom, INT width, INT height, INT side, INT orientation, INT leaf_size )
{
  VECTOR *v, *vpred, *end_v, **v0, **v1, **vi;
  BLOCKVECTOR *bv0, *bv1, *bvi;
  register INT index, interface, nr_0, nr_1, nr_i;
  VECTOR *pred_first_vec;

    #ifndef __BLOCK_VECTOR_DESC__
  return(1);
    #endif

  /* create a double linked list of 3 blockvectors */
  v = BVFIRSTVECTOR( bv );
  end_v = BVENDVECTOR( bv );
  pred_first_vec = PREDVC( v );

  if ( CreateBlockvector( grid, &bv0 ) != GM_OK )
    return GM_OUT_OF_MEM;

  BVPRED( bv0 ) = NULL;
  SETBVDOWNTYPE( bv0, BVDOWNTYPEVECTOR );
  SETBVTVTYPE( bv0, BV1DTV );
  BVNUMBER( bv0 ) = 0;
  v0 = &BVDOWNVECTOR( bv0 );

  if ( CreateBlockvector( grid, &bv1 ) != GM_OK )
  {
    DisposeBlockvector( grid, bv0 );
    return GM_OUT_OF_MEM;
  }

  BVSUCC( bv0 ) = bv1;
  BVPRED( bv1 ) = bv0;
  SETBVDOWNTYPE( bv1, BVDOWNTYPEVECTOR );
  SETBVTVTYPE( bv1, BV1DTV );
  BVNUMBER( bv1 ) = 1;
  v1 = &BVDOWNVECTOR( bv1 );

  if ( CreateBlockvector( grid, &bvi ) != GM_OK )
  {
    DisposeBlockvector( grid, bv1 );
    DisposeBlockvector( grid, bv0 );
    return GM_OUT_OF_MEM;
  }

  BVSUCC( bv1 ) = bvi;
  BVPRED( bvi ) = bv1;
  SETBVDOWNTYPE( bvi, BVDOWNTYPEVECTOR );
  SETBVTVTYPE( bvi, BV1DTV );
  BVNUMBER( bvi ) = 2;
  BVSUCC( bvi ) = NULL;
  vi = &BVDOWNVECTOR( bvi );

  /* all necessary memory allocated; insert new block list in father block */
  SETBVDOWNTYPE( bv, BVDOWNTYPEBV );
  SETBVTVTYPE( bv, BV1DTV );
  BVDOWNBV( bv ) = bv0;
  BVDOWNBVLAST( bv ) = bvi;

  /* coordinate of the interface line */
  interface = ( (orientation == BV_VERTICAL) ? (width-1)/2 + left : bottom + (height-1)/2 );

  /* sort each vector in one of the 3 subdomains */
  nr_0 = nr_1 = nr_i = 0;
  for ( vpred = NULL; v != end_v; vpred = v, v = SUCCVC( v ) )
  {
    /* calculate row/column number from index */
    index = ( (orientation == BV_VERTICAL) ? VINDEX(v) % side : VINDEX(v) / side );

    if ( index < interface )
    {
      /* vector belongs to subdomain 0 */
      *v0 = v;
      v0 = &SUCCVC( v );
                        #ifdef __BLOCK_VECTOR_DESC__
      BVD_PUSH_ENTRY( &VBVD( v ), 0, &DH_bvdf );
            #endif
      nr_0++;
    }
    else if ( index > interface )
    {
      /* vector belongs to subdomain 1 */
      *v1 = v;
      v1 = &SUCCVC( v );
                        #ifdef __BLOCK_VECTOR_DESC__
      BVD_PUSH_ENTRY( &VBVD( v ), 1, &DH_bvdf );
            #endif
      nr_1++;
    }
    else
    {
      /* vector belongs to subdomain 2 (interface) */
      *vi = v;
      vi = &SUCCVC( v );
                        #ifdef __BLOCK_VECTOR_DESC__
      BVD_PUSH_ENTRY( &VBVD( v ), 2, &DH_bvdf );
            #endif
      nr_i++;
    }

  }

  /* relink the 3 subdomain lists to the global vector list */
  *v0 = BVDOWNVECTOR( bv1 );
  *v1 = BVDOWNVECTOR( bvi );
  *vi = end_v;
  v = BVDOWNVECTOR( bv0 );
  vpred = pred_first_vec;
  if ( pred_first_vec != NULL )
    SUCCVC( pred_first_vec ) = v;
  /* correct the PRED-pointer chain */
  for ( ; v != end_v; vpred = v, v = SUCCVC( v ) )
    PREDVC( v ) = vpred;
  if ( v != NULL )              /* at the end of the vectors of this blockvector */
    PREDVC( v ) = vpred;

  BVNUMBEROFVECTORS( bv0 ) = nr_0;
  BVNUMBEROFVECTORS( bv1 ) = nr_1;
  BVNUMBEROFVECTORS( bvi ) = nr_i;

  /* complete the relation between blockvectors and vectors */
  BVLASTVECTOR( bv0 ) = BVPRED(BVDOWNVECTOR( bv1 ));
  BVLASTVECTOR( bv1 ) = BVPRED(BVDOWNVECTOR( bvi ));
  BVLASTVECTOR( bvi ) = vpred;
  BVLASTVECTOR( bv ) = vpred;

  BVFIRSTVECTOR( bv0 ) = BVDOWNVECTOR( bv0 );
  BVFIRSTVECTOR( bv1 ) = BVDOWNVECTOR( bv1 );
  BVFIRSTVECTOR( bvi ) = BVDOWNVECTOR( bvi );
  BVFIRSTVECTOR( bv ) = BVFIRSTVECTOR( bv0 );

  /* next loop */

  if ( orientation == BV_VERTICAL )
  {
    /* new dimensions */
    width = ( width - 1 ) / 2;

    /* if there are too many points per subdomain halfen them */
    if ( ( width * height ) > leaf_size )
    {
      if ( BlockHalfening( grid, bv0, left, bottom, width, height, side, BV_HORIZONTAL, leaf_size ) == GM_OUT_OF_MEM )
        return GM_OUT_OF_MEM;
      if ( BlockHalfening( grid, bv1, left+width+1, bottom, width, height, side, BV_HORIZONTAL, leaf_size ) == GM_OUT_OF_MEM )
        return GM_OUT_OF_MEM;
    }
  }
  else       /* orientation == BV_HORIZONTAL */
  {
    /* new dimensions */
    height = ( height - 1 ) / 2;

    /* if there are too many points per subdomain halfen them */
    if ( ( width * height ) > leaf_size )
    {
      if ( BlockHalfening( grid, bv0, left, bottom, width, height, side, BV_VERTICAL, leaf_size ) == GM_OUT_OF_MEM )
        return GM_OUT_OF_MEM;
      if ( BlockHalfening( grid, bv1, left, bottom +height+1, width, height, side, BV_VERTICAL, leaf_size ) == GM_OUT_OF_MEM )
        return GM_OUT_OF_MEM;
    }
  }

  return GM_OK;
}
#endif /* __BLOCK_VECTOR_DESC__ */


/****************************************************************************/
/*D
   MoveVector - Move vector within the vector list of the grid

   SYNOPSIS:
   INT MoveVector (GRID *theGrid, VECTOR *moveVector,
   VECTOR *destVector, INT after);

   PARAMETERS:
   .  theGrid - pointer to grid
   .  moveVector - vector to move
   .  destVector - vector before or after which 'moveVector' will be inserted
   .  after - true (1) or false (0)

   DESCRIPTION:
   This function moves a 'VECTOR' within the double linked list. If 'after' is
   true then 'moveVector' will be inserted immediately after 'destVector', if
   'after' is false then it will be inserted immediately before. If 'destVector'
   is 'NULL' then 'moveVector' is inserted at the beginning when 'after' is true
   and at the end of the list when 'after' is false.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT MoveVector (GRID *theGrid, VECTOR *moveVector, VECTOR *destVector, INT after)
{
  if (theGrid==NULL || moveVector==NULL) return (1);
  if (moveVector==destVector) return (0);

  /* take vector out of list */
  if (PREDVC(moveVector)!=NULL) SUCCVC(PREDVC(moveVector))      = SUCCVC(moveVector);
  else SFIRSTVECTOR(theGrid)           = SUCCVC(moveVector);
  if (SUCCVC(moveVector)!=NULL) PREDVC(SUCCVC(moveVector))      = PREDVC(moveVector);
  else LASTVECTOR(theGrid)             = PREDVC(moveVector);

  /* put it in list */
  if (destVector!=NULL)
  {
    if (after)
    {
      if (SUCCVC(destVector)!=NULL) PREDVC(SUCCVC(destVector))      = moveVector;
      else LASTVECTOR(theGrid)             = moveVector;
      SUCCVC(moveVector) = SUCCVC(destVector);
      PREDVC(moveVector) = destVector;
      SUCCVC(destVector) = moveVector;
    }
    else
    {
      if (PREDVC(destVector)!=NULL) SUCCVC(PREDVC(destVector))      = moveVector;
      else SFIRSTVECTOR(theGrid)           = moveVector;
      PREDVC(moveVector) = PREDVC(destVector);
      SUCCVC(moveVector) = destVector;
      PREDVC(destVector) = moveVector;
    }
  }
  else
  {
    if (after)
    {
      SUCCVC(moveVector) = SFIRSTVECTOR(theGrid);
      PREDVC(moveVector) = NULL;
      SFIRSTVECTOR(theGrid) = moveVector;
      if (SUCCVC(moveVector)!=NULL) PREDVC(SUCCVC(moveVector)) = moveVector;
    }
    else
    {
      SUCCVC(moveVector) = NULL;
      PREDVC(moveVector) = LASTVECTOR(theGrid);
      LASTVECTOR(theGrid) = moveVector;
      if (PREDVC(moveVector)!=NULL) SUCCVC(PREDVC(moveVector)) = moveVector;
    }
  }

  return (0);
}

#ifdef __INTERPOLATION_MATRIX__

/****************************************************************************/
/*D
   CreateIMatrix -  Return pointer to a new interpolation matrix structure

   SYNOPSIS:
   MATRIX *CreateIMatrix (GRID *theGrid, VECTOR *fvec, VECTOR *cvec);

   PARAMETERS:
   .  theGrid - fine grid
   .  fvec - fine grid vector
   .  cvec - coarse grid vector

   DESCRIPTION:
   This function allocates a new 'MATRIX' structures in the
   'imatrix' list of 'fvec'.

   RETURN VALUE:
   MATRIX *
   .n    pointer to the new matrix
   .n    NULL if error occured.
   D*/
/****************************************************************************/

MATRIX *CreateIMatrix (GRID *theGrid, VECTOR *fvec, VECTOR *cvec)
{
  MULTIGRID *theMG;
  HEAP *theHeap;
  MATRIX *pm;
  INT RootType, DestType, MType, ds, Size;

  pm = GetIMatrix(fvec,cvec);
  if (pm != NULL)
    return(pm);

  ASSERT(fvec != NULL);
  ASSERT(cvec != NULL);

  RootType = VTYPE(fvec);
  DestType = VTYPE(cvec);
  MType = MatrixType[RootType][DestType];

  /* check expected size */
  theMG = MYMG(theGrid);
  theHeap = theMG->theHeap;
  ds = FMT_S_IMAT_TP(MGFORMAT(theMG),MType);
  if (ds == 0)
    return (NULL);
  Size = sizeof(MATRIX)-sizeof(DOUBLE)+ds;
  if (MSIZEMAX<Size) return (NULL);
  assert (Size % ALIGNMENT == 0);
  pm = (MATRIX *)GetMemoryForObject (theMG,Size,MAOBJ);
  if (pm==NULL)
    return (NULL);

  SETOBJT(pm,MAOBJ);
  SETMDIAG(pm,0);
  SETMROOTTYPE(pm,RootType);
  SETMDESTTYPE(pm,DestType);
  SETMSIZE(pm,Size);
  MDEST(pm) = cvec;
  MNEXT(pm) = VISTART(fvec);
  VISTART(fvec) = pm;

  /* counters */
  theGrid->nIMat++;
  return(pm);
}

/****************************************************************************/
/*D
   DisposeIMatrices - Remove interpolation matrix from the data structure

   SYNOPSIS:
   static INT DisposeIMatrices (GRID *theGrid, MATRIX *theMatrix);

   PARAMETERS:
   .  theGrid - the grid to remove from
   .  theMatrix - start of matrix list to dispose

   DESCRIPTION:
   This function removes an interpolation  matrix list from the data
   structure.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT DisposeIMatrices (GRID *theGrid, MATRIX *theMatrix)
{
  MATRIX *Matrix, *NextMatrix;

  for (Matrix=theMatrix; Matrix!=NULL; )
  {
    NextMatrix = NEXT(Matrix);
    PutFreeObject(theGrid->mg,Matrix,MSIZE(Matrix),MAOBJ);
    theGrid->nIMat--;
    Matrix = NextMatrix;
  }

  return(0);
}

INT DisposeIMatrixList (GRID *theGrid, VECTOR *theVector)
{
  if (DisposeIMatrices(theGrid,VISTART(theVector)))
    RETURN (1);
  VISTART(theVector) = NULL;
  return (0);
}

INT DisposeIMatricesInGrid (GRID *theGrid)
{
  VECTOR *theV;

  for (theV=PFIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    if (DisposeIMatrices(theGrid,VISTART(theV))) RETURN (1);
    VISTART(theV) = NULL;
  }

  return (0);
}

/****************************************************************************/
/*D
   GetIMatrix - Return pointer to interpolation matrix if it exists

   SYNOPSIS:
   MATRIX *GetIMatrix (VECTOR *FineVector, VECTOR *CoarseVector);

   PARAMETERS:
   .  FineVector - fine grid vector
   .  CoarseVector - coarse grid vector

   DESCRIPTION:
   This function returns pointer to interpolation matrix.
   If it does not exist already, it returns NULL.

   RETURN VALUE:
   MATRIX *
   .n       pointer to Matrix,
   .n       NULL if no Matrix exists.
   D*/
/****************************************************************************/

MATRIX *GetIMatrix (VECTOR *FineVector, VECTOR *CoarseVector)
{
  MATRIX *theMatrix;

  for (theMatrix=VISTART(FineVector); theMatrix!=NULL;
       theMatrix = MNEXT(theMatrix))
    if (MDEST(theMatrix)==CoarseVector)
      return (theMatrix);

  return (NULL);
}

#endif
/* __INTERPOLATION_MATRIX__ */

/****************************************************************************/
/*
   InitAlgebra - Init algebra

   SYNOPSIS:
   INT InitAlgebra (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function inits algebra.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

INT InitAlgebra (void)
{
  INT i, j, n;

  /* install the /Alg Dep directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitAlgebra","could not changedir to root");
    return(__LINE__);
  }
  theAlgDepDirID = GetNewEnvDirID();
  if (MakeEnvItem("Alg Dep",theAlgDepDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitAlgebra","could not install '/Alg Dep' dir");
    return(__LINE__);
  }
  theAlgDepVarID = GetNewEnvVarID();

  /* install the /FindCut directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitAlgebra","could not changedir to root");
    return(__LINE__);
  }
  theFindCutDirID = GetNewEnvDirID();
  if (MakeEnvItem("FindCut",theFindCutDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitAlgebra","could not install '/FindCut' dir");
    return(__LINE__);
  }
  theFindCutVarID = GetNewEnvVarID();

  /* set MatrixType-field */
  /* example for MAXVECTORS=4
     0 4 7 9
     4 1 5 8
     7 5 2 6
     9 8 6 3
   */
  n=0;
  for (i=0; i<MAXVECTORS; i++)
    for (j=0; j<MAXVECTORS-i; j++)
    {
      MatrixType[j][j+i] = n;
      MatrixType[j+i][j] = n;
      n++;
    }
  if (n != MAXMATRICES)
    return (__LINE__);

  if (n+MAXVECTORS != MAXCONNECTIONS)
    return (__LINE__);

  /* init standard algebraic dependencies */
  if (CreateAlgebraicDependency ("lex",LexAlgDep)==NULL) return(__LINE__);
  if (CreateAlgebraicDependency ("stronglex",StrongLexAlgDep)==NULL) return(__LINE__);

  /* init default find cut proc */
  if (CreateFindCutProc ("lex",FeedbackVertexVectors)==NULL) return(__LINE__);

  for (i=0; i<MAXVOBJECTS; i++)
    switch (i)
    {
    case NODEVEC : ObjTypeName[i] = "nd"; break;
    case EDGEVEC : ObjTypeName[i] = "ed"; break;
    case ELEMVEC : ObjTypeName[i] = "el"; break;
    case SIDEVEC : ObjTypeName[i] = "si"; break;
    default : ObjTypeName[i] = "";
    }

  return (0);
}
