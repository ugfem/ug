// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gridt.h														*/
/*																			*/
/* Purpose:   header file for gridt.c										*/
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
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

#ifndef __NGIN__
#define __NGIN__

#define __USE_IN_UG__

#ifdef __USE_IN_UG__
#include "general.h"
#include "heaps.h"
#include "lgm_transfer.h"
#include "ugdevices.h"
#endif

#include "namespace.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define NG_LPMAX                40
#define ERROR_PREFIX    "\nngin: parse-error: "

typedef struct {
  INT line_id;
  FLOAT local;
} LINE_POSITION;

typedef struct {
  INT n_lp;
  LINE_POSITION lp[NG_LPMAX];
  DOUBLE global[2];
} BND_NODE;

typedef struct {
  DOUBLE global[2];
} INNER_NODE;

typedef struct {
  INT c_id[2];
} ELEM_SIDE;

typedef struct {
  INT subdom;
  INT n_c;
  INT c_id[4];
  INT n_s;
  ELEM_SIDE side[4];
} NG_ELEMENT;

#ifdef __USE_IN_UG__

#define NG_MALLOC(h,s,m)       GetTmpMem(h,s,m)
#define NG_Print               UserWriteF
#define NG_FOPEN(d,n)           {if (lgmdomainpathes_set) \
                                   (d)=FileOpenUsingSearchPaths(n,"r","lgmdomainpathes"); \
                                 else (d)=fopen(n,"r");}

#else

#define NG_MALLOC(h,s,m)       malloc(s)
#define NG_Print               printf
#define NG_FOPEN(d,n)          {(d)=fopen(n,"r");}
typedef INT HEAP;
typedef struct {

  INT nBndP;                         /* nb. of boundary points                  */
  INT *BndP_nLine;                   /* nb. of lines per bound. point                   */
  INT **BndP_LineID;                 /* id of each line                                 */
  FLOAT **BndP_lcoord;               /* local coord of BndP on each line                */
  FLOAT **BndPosition;               /* list of boundary points                 */
  INT nInnP;                         /* nb. of inner nodes                      */
  DOUBLE **InnPosition;               /* positions of inner nodes               */
  INT nSubDomains;                   /* nb. of subdomains                       */
  INT *nSides;                       /* nb. of boundary sides per subdomain     */
  INT ***Side_corner_ids;            /* corner ids                              */
  INT *nElements;                    /* nb. of element corners                  */
  INT **Element_corners;             /* nb. of element corners                  */
  INT **Element_SideOnBnd;           /* used bitwise: sides on bnd for elem     */
  INT ***Element_corner_ids;         /* nb. of side corners                     */
  INT ***nbElements;                 /* nb. of side corners                     */

} LGM_MESH_INFO;

#error MIST!!

#endif

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

int PutBndNode (BND_NODE *BndNode);
int PutInnerNode (INNER_NODE *InnNode);
int PutElement (NG_ELEMENT *Elem);
int NG_ReadMesh (const char *name, HEAP *Heap, LGM_MESH_INFO *theMesh, INT MarkKey);

#ifdef __USE_IN_UG__
int NG_Init (int domainpaths_set);
#endif

void ngbreak (void);

END_UGDIM_NAMESPACE

#endif
