// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ff_gen.h														*/
/*																			*/
/* Purpose:   general frequency filtering decompostion routines             */
/*                numproc mechanism												*/
/*																			*/
/* Author:	  Christian Wrobel                                                                              */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
/*																			*/
/* History:   16.09.96 begin, ug version 3.3								*/
/*																			*/
/* Remarks:   FF is used as the abbreviation for "frequency filtering"		*/
/*                        TFF is used as the abbreviation for the "tangential frequency */
/*			       filtering" method due to Christian Wagner, 1995			*/
/*																			*/
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __FF_GEN__
#define __FF_GEN__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* if defined, special BV communication routines are used */
#ifdef FF_ModelP
#define FFCOMM
#endif

#ifndef M_LN2
#define M_LN2         0.69314718055994530942
#endif

#define DUMMY_COMP -1

/* number of FF decompositions to be stored as decomposed matrixes */
#define NUMBER_FF_DECOMPS 1

/* max. depth of matrix hierarchy */
#define FF_MAX_MATS 10

/* max. number of auxiliary vectors */
#define FF_MAX_VECS 20


#define STIFFMAT_ON_LEVEL(bv)                           (FF_Mats[BVLEVEL(bv)])
#define DECOMPMAT_ON_LEVEL(bv)                          (FF_Mats[BVLEVEL(bv)+1])
#define DECOMP_MATDATA_DESC_ON_LEVEL(bv)                        (FF_MATDATA_DESC_ARRAY[BVLEVEL(bv)+1])
/* if you are already on the level of the single blocks */
#define STIFFMAT_ON_LEVEL_BLOCKWISE(bv)         (FF_Mats[BVLEVEL(bv)-1])
#define DECOMPMAT_ON_LEVEL_BLOCKWISE(bv)        (FF_Mats[BVLEVEL(bv)])

/* realizes a simple stack of aux vectors in the array FF_Vecs
   CAUTION: the sequence of free's must be exctly the reverse of the get's! */
#ifdef Debug
#define GET_AUX_VEC                                     ( (TOS_FF_Vecs<FF_MAX_VECS) ? FF_Vecs[TOS_FF_Vecs++] : -1 )
/* tricky: enforce, that TOS_FF_Vecs-- and ASSERT(FALSE) have the same type */
#define FREE_AUX_VEC(vec)                       ( (vec==FF_Vecs[TOS_FF_Vecs-1]) ? (TOS_FF_Vecs--) : 0,ASSERT(FALSE) );
#else
#define GET_AUX_VEC                                     (FF_Vecs[TOS_FF_Vecs++])
#define FREE_AUX_VEC(vec)                       (TOS_FF_Vecs--);
#endif

/* solve subproblems in MultWithMInv exactly by recursive solving */
#define MINV_EXACTQQQ

#ifdef MINV_EXACT
#define SOLVE_ITERQQQ
#endif

/* stop iteration if residuum is damped by eps */
/* usually one iterates until the res reaches a given accuracy */
#define EPS_RELATIVEQQQ

/* solve in calc. theta A^-1 * B * tv exactly */
#define TV_INV_EXACTQQQ

/* use the harmonic extended topmost testvector as  the global tv */
#define TV_HARMONIC_EXTENSIONQQQQ

/* use the exact solution as testvector */
#define TV_EXACT_SOLUTIONQQQ

/* extends global tv to the subdomains to fulfil the global filtering condition */
#define EXTEND_GLOBAL_TVQQQ

/* shift all frequ (but not the first ) 1 up */
#define SHIFT_FREQQQQ

#define BACKUP_MEM(bv) (DOUBLE*)BVUSERDATA(bv)

#define CHECK_CALCULATIONQQQ

/* invert for Theta calc exactly */
#define THETA_EXACTQQQ

/* invert 2D subproblem in MultWithMInv exactly */
#define MINV_2D_EXACTQQQ

/* determine theta analytically as eigenvalue */
#define THETA_ANAQQQ

/* macros for communication buffer */
#define FFMAX_TRIES     500000  /* max. number of tries til timeout in communication */

#define FF_LINES_NR -100
#define FF_CROSS_NR -101

#ifdef ModelP

#define FFCommBufferBV(bv)      ((FFCommBuffer*)BVUSERDATA(bv))
#define FFVChannelIn(bv)    (FFCommBufferBV(bv)->vcin)
#define FFVChannelOut(bv)   (FFCommBufferBV(bv)->vcout)
#define FFMsgIn(bv)             (FFCommBufferBV(bv)->msgIn)
#define FFMsgOut(bv)            (FFCommBufferBV(bv)->msgOut)
#define FFBufferIn(bv)      (FFCommBufferBV(bv)->inbuffer)
#define FFBufferOut(bv)     (FFCommBufferBV(bv)->outbuffer)
#define FFSteps(bv)             (FFCommBufferBV(bv)->steps)
#define FFStepCounter(bv)   (FFCommBufferBV(bv)->stepcounter)

#define FFHasCommBuffer(bv) ((bv)!=NULL && FFCommBufferBV(bv)!=NULL && FFBufferIn(bv)!=NULL && FFBufferOut(bv)!=NULL && FFVChannelIn(bv)!=NULL && FFVChannelOut(bv)!=NULL)

#define FFCommTridiagMatSize(bv)        ((3 * BVNUMBEROFVECTORS(bv) - 2) * sizeof(DOUBLE))
#define FFCommVecSize(bv)                       (BVNUMBEROFVECTORS(bv) * sizeof(DOUBLE))
#define FFCommSize(bv,ffcommobjt)       ((ffcommobjt)==FFCommVec ? FFCommVecSize(bv) : (ffcommobjt)==FFCommTridiagMat ? FFCommTridiagMatSize(bv) : -1)
#define FFCommMaxSize(bv)                       FFCommTridiagMatSize(bv)

#define FF_MAX_CROSS_PER_PE             4
#define FF_CROSS_MAX_NEIGHBOURS 3
#define FF_MAX_CROSS_MATS               (FF_MAX_CROSS_PER_PE*FF_CROSS_MAX_NEIGHBOURS)

#define FFCROSSMAT(r,c)                         FFCrossMatMem[2*FFCrossBw*(r)+(c)]

#define FFCrossVC(pe)                           (PEInfoArray[pe].vc)
#define FFCrossRow(pe)                          (PEInfoArray[pe].row)
#define FFCrossCol(pe)                          (PEInfoArray[pe].col)
#define FFCrossFirstPos(pe)                     (PEInfoArray[pe].fpos)
#define FFCrossOffsetPtr(pe)            (PEInfoArray[pe].offsets)
#define FFCrossVecs(pe)                         (FFCrossOffsetPtr(pe)->numberofvectors)
#define FFCrossNbs(pe)                          (FFCrossOffsetPtr(pe)->numberofneighbours)
#define FFCrossVecSize(pe)                      (FFCrossOffsetPtr(pe)->sizeofvector)
#define FFCrossMatSize(pe)                      (FFCrossOffsetPtr(pe)->sizeofmatrix)
#define FFCrossVecPosOffset(pe,vnr)             (FFCrossOffsetPtr(pe)->v_rel_pos[vnr])
#define FFCrossNbPosOffset(pe,vnr,nbnr) (FFCrossOffsetPtr(pe)->nb_rel_pos[vnr][nbnr])

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct ffcommbuffer
{
  VChannelPtr vcin, vcout;
  msgid msgIn, msgOut;
  INT steps, stepcounter;
  char *inbuffer, *outbuffer;
  DOUBLE bufferspace[1];                /* will be extended at allocation time! */
  /* type DOUBLE neccessary for maximal alignment */
};

typedef struct ffcommbuffer FFCommBuffer;
typedef enum ffcommobject {FFCommNone,FFCommVec,FFCommTridiagMat} FFCommObjectType;
typedef INT (*FFBufferFunc)(BLOCKVECTOR *bv, FFCommObjectType ffcommobjt);

struct ffcrossoffsets
{
  INT numberofvectors;          /* number of cross vectors in the pe */
  INT numberofneighbours;       /* number of neighbours in the matrix graph of each cross point in the pe */
  INT sizeofvector;                     /* buffer size to store all cross vectors */
  INT sizeofmatrix;                     /* buffer size to store all cross point matrix */
  INT v_rel_pos[FF_MAX_CROSS_PER_PE];           /* relative position of the cross vectors (in their lex. ordering) resp. the lex. first cross vector in this pe */
  INT nb_rel_pos[FF_MAX_CROSS_PER_PE][FF_CROSS_MAX_NEIGHBOURS];         /* relative position of the neighbours of the cross vectors (in pred/succ ordering) resp. the according cross vector */
};
typedef struct ffcrossoffsets FFCrossOffsets;

struct peinfo
{
  VChannelPtr vc;                               /* channel to the slave pe */
  INT row;                                              /* row of this pe in the domain decomposition grid */
  INT col;                                              /* column of this pe in the domain decomposition grid */
  INT fpos;                                             /* global number of the lex. first cross point vector in this pe */
  FFCrossOffsets *offsets;              /* pointer to the offsets of the further cross points within this pe */
};
typedef struct peinfo PEInfo;

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

extern PEInfo *PEInfoArray;             /* dyn. allocacted array for each slave pe */
extern int FFFarmer;                    /* pe number of the master for cross system calculation */
extern VChannelPtr FFFarmerVC;  /* for slave: channel to the farmer */
extern DOUBLE *FFCrossVecMem;   /* mem for global cross point vector */
extern INT FFCrossVecSize;              /* size of the FFCrossVecMem vector */
extern DOUBLE *FFCrossMatMem;   /* mem for band matrix for cross point system */
extern INT FFCrossMatSize;              /* size of the FFCrossMatMem matrix */
extern INT FFCrossBw;                   /* bandwidth of the FFCrossMatMem matrix */
extern DOUBLE FFFarmerBuffer[FF_MAX_CROSS_MATS];        /* buffer for communication with FFFarmer */
extern INT FFCrossOffdiagMats;  /* number of offdiagonal matrix entries for cross system */
#endif

extern INT mute_level;

/* value below them a division is refused calculating a testvector */
extern DOUBLE FFsmallTV;

/* ratio for a jump to be detected */
extern DOUBLE FFmuchBigger;

/* value below them a number is considered as 0.0 */
extern DOUBLE FFEPS;

/* value below them an approximation error is considered as ok */
extern DOUBLE FFaccuracy;

/* global array to hold the matrix hierarchy */
extern INT FF_Mats[FF_MAX_MATS];
extern MATDATA_DESC *FF_MATDATA_DESC_ARRAY[FF_MAX_MATS];

/* global array to hold the auxiliary vectors */
extern INT FF_Vecs[FF_MAX_VECS];
extern INT TOS_FF_Vecs;


/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
extern INT aux2_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
extern INT aux3_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
extern INT aux4_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
extern INT aux5_COMP;

/* auxiliary component; only for checking the results (if CHECK_CALCULATION is on) */
extern INT aux6_COMP;


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

#ifdef ModelP
INT FFStartDeallocBuffer( BLOCKVECTOR *bv, FFCommObjectType ffcommobjt );
INT FFFinishDeallocBuffer( BLOCKVECTOR *bv, FFCommObjectType ffcommobjt );

INT FFStartComm( BLOCKVECTOR *bv, FFBufferFunc b_func, FFCommObjectType ffcommobjt );
void FFFinishComm( BLOCKVECTOR *bv, FFBufferFunc b_func, FFCommObjectType ffcommobjt, INT calls );
void FFMakeConsistent( BLOCKVECTOR *bv, FFCommObjectType ffcommobjt, INT comp );
#define FFVectorConsistent(bv,v_comp)           (FFMakeConsistent(bv,FFCommVec,v_comp))
#define FFTridiagMatConsistent(bv,m_comp)       (FFMakeConsistent(bv,FFCommTridiagMat,m_comp))

void FFInitCrossComm( BLOCKVECTOR *bv );
void FFFinishCrossComm( void );
void FFGatherCrossMat( BLOCKVECTOR *bv, INT m_comp );
void FFSolveCrossVec( BLOCKVECTOR *bv, INT v_comp );

#endif



INT     storeVectorBS( BLOCKVECTOR *bv, INT x_comp, GRID *grid );
INT     restoreVectorBS( BLOCKVECTOR *bv, INT x_comp );

void printv( INT x_nr );
void printvgrid( GRID *g, INT x_nr );
void printvBS( const BLOCKVECTOR *bv, INT x_nr );
void printm( INT m_nr );
void printmgrid( GRID *g, INT m_nr );
void printmMG( MULTIGRID *theMG, INT m_nr );
void printmBS( const BLOCKVECTOR *bv_row, const BLOCKVECTOR *bv_col, INT m_nr );
void printPatternBS( const BLOCKVECTOR *bv_row, const BLOCKVECTOR *bv_col, INT m_nr );

void printBV( const BV_DESC_FORMAT *bvdf );
void printBVgrid( GRID *grid, const BV_DESC_FORMAT *bvdf );

DOUBLE FFMeshwidthOfGrid( GRID *grid );
INT FF_PrepareGrid( GRID *grid, DOUBLE *meshwidth, INT init, INT K_comp, INT x_comp, INT b_comp, const BV_DESC_FORMAT *bvdf );

void FFConstructTestvector( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr, DOUBLE wavenr_3D );
void FFConstructTestvector_loc( const BLOCKVECTOR *bv, INT tv_comp, DOUBLE wavenr, DOUBLE wavenr_3D );

INT FFMultWithM( const BLOCKVECTOR *bv, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf, INT y_comp, INT x_comp );

INT FFMultWithMInv(
  const BLOCKVECTOR *bv,
  const BV_DESC *bvd,
  const BV_DESC_FORMAT *bvdf,
  INT v_comp,
  INT b_comp
#ifdef ModelP
  ,const VECDATA_DESC *v,
  GRID *grid
#endif
  );

INT InitFF (void);

#endif
