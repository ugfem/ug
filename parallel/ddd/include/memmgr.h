// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      memmgr.h                                                      */
/*                                                                          */
/* Purpose:   header file for ddd memory management interface               */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/*                                                                          */
/* History:   94/04/27 kb  begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __MEMMGR__
#define __MEMMGR__
#ifdef __cplusplus
extern "C" {
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



/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/




/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/




/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/




/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


void *memmgr_AllocOMEM (size_t size, int typeid , int prio, int attr);
void  memmgr_FreeOMEM (void *mem, size_t size, int typeid );

void *memmgr_AllocPMEM (size_t size);
void  memmgr_FreePMEM (void *mem);

void *memmgr_AllocAMEM (size_t size);
void  memmgr_FreeAMEM (void *mem);

void *memmgr_AllocTMEM (size_t size);
void  memmgr_FreeTMEM (void *mem);

void *memmgr_AllocHMEM (size_t size);
void  memmgr_FreeHMEM (void *mem);

void  memmgr_MarkHMEM (void);
void  memmgr_ReleaseHMEM (void);


#ifdef __cplusplus
}
#endif
#endif
