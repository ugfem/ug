// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      dddresources.h                                                */
/*                                                                          */
/* Purpose:   basic ddd resource manager                                    */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 3a                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   951124 kb  begin                                              */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD_RESOURCE_H__
#define __DDD_RESOURCE_H__

#ifdef __cplusplus
extern "C" {
#endif



/****************************************************************************/
/*                                                                          */
/* data structures and new types                                            */
/*                                                                          */
/****************************************************************************/


/*
        basic resource descriptor
 */
typedef struct
{
  int maxObjs;
  int maxCpls;
  int nObjs;
  int nCpls;
} DDD_RESOURCES;



/****************************************************************************/
/*                                                                          */
/* declaration of DDD functional interface                                  */
/*                                                                          */
/****************************************************************************/


/*
        Resource Manager Module
 */
DDD_RESOURCES *DDD_InfoResources (void);


#ifdef __cplusplus
}
#endif
#endif
