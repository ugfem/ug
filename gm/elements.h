// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file elements.h
 * \ingroup gm
 */

/****************************************************************************/
/*                                                                          */
/* File:      elements.h                                                    */
/*                                                                          */
/* Purpose:   general element concept (header file)                         */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   24.03.95 begin, ug version 3.0                                */
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

#ifndef __ELEMENTS__
#define __ELEMENTS__

#include "namespace.h"

START_UGDIM_NAMESPACE

/** \brief Number of different element types    */
#define TAGS 8

extern INT n_offset[TAGS];
extern INT father_offset[TAGS];
extern INT sons_offset[TAGS];
extern INT nb_offset[TAGS];
extern INT evector_offset[TAGS];
extern INT svector_offset[TAGS];
extern INT side_offset[TAGS];
extern INT data_offset[TAGS];

/* the element descriptions are also globally available, these are pointers ! */
extern GENERAL_ELEMENT *element_descriptors[TAGS];
extern GENERAL_ELEMENT *reference_descriptors[MAX_CORNERS_OF_ELEM+1];
extern INT reference2tag[MAX_CORNERS_OF_ELEM+1];


/****************************************************************************/
/*                                                                          */
/* function definitions                                                     */
/*                                                                          */
/****************************************************************************/

INT PreInitElementTypes         (void);
INT InitElementTypes            (MULTIGRID *theMG);


END_NAMESPACE


#endif
