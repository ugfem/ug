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

#ifndef M_LN2
#define M_LN2         0.69314718055994530942
#endif

#define DUMMY_COMP -1

/* number of FF decompositions to be stored as decomposed matrixes */
#define NUMBER_FF_DECOMPS 1

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

extern INT mute_level;

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

INT InitFF (void);

#endif
