// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/***************************************************************************/
/*                                                                          */
/* File:      lb4.h                                                         */
/*                                                                          */
/* Purpose:   clustered partitioning with several strategies header file    */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr1.iwr.uni-heidelberg.de                  */
/*                                                                          */
/*            Stefan Lang                                                   */
/*            Institut fuer Mathematische Maschinen und                     */
/*            Datenverarbeitung III                                         */
/*            Universitaet Erlangen-Nuernberg                               */
/*            Martensstrasse 3                                              */
/*            91058 Erlangen                                                */
/*                                                                          */
/* History:   10 Jan 1994, begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __LB4__
#define __LB4__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGTYPES2__
/* #include "ugtypes2.h" */
#endif
#ifndef __GM_H__
#include "gm.h"
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAXDEPTH            8           /* max depth of a cluster           */


/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/


typedef struct cluster {
  SHORT source;                             /* where cluster resides currently  */
  SHORT destination;                        /* where cluster is mapped to       */
  SHORT minlevel;                           /* where cluster starts             */
  SHORT depth;                              /* maxlevel = minlevel + depth      */
  unsigned SHORT level_size[MAXDEPTH];      /* weight on each level (from minl.)*/
  COORD sx,sy;                              /* center of mass                   */
                #ifdef __THREEDIM__
  COORD sz;                                 /* center of mass                   */
                #endif
  ELEMENT  *root_element;                   /* root element for dual graph      */
  /* following part of struct is similar to struct vtx_data of Chaco      */
  int size;                                 /* total weight of cluster          */
  int nedges;                               /* number of neighbors              */
                                            /* above field includes self-edge   */
  int edges[MAX_SIDES_OF_ELEM+1];                       /* dual graph of grid on minlevel   */
  /* float ewgts[SIDES+1]; */              /* weights of all the edges         */
  /*above 2fields have self-edge first*/
} CLUSTER ;


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

int Balance_CCPTM (MULTIGRID *mg,                                 /* data    */
                   int minlevel, int cluster_depth, int threshold, /* clusters*/
                   int Const,                                                    /* balance */
                   int element_limit, int channel_limit,                         /*transfer*/
                   int strategy,                                                 /* strategy*/
                   int eigen,                                                    /* eigensol*/
                   int loc,                                                      /* KL?     */
                   int dims,                                                     /*dimension*/
                   int weights,                                                  /* weights?*/
                   int coarse,                                                   /* #coarse */
                   int mode,                                                     /*RSBtype */
                   int iter);                                                    /*Transfer*/


#endif
