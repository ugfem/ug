// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file GenerateRules.h
 * \ingroup gm
 */

/****************************************************************************/
/*                                                                                                                                                      */
/* File:        GenerateRules.h                                                                                                 */
/*                                                                                                                                                      */
/* Purpose: header file for GenerateRules.c                                                             */
/*                                                                                                                                                      */
/* Author:      Henrik Reichert                                                                                                 */
/*                      Institut fuer Angewandte Mathematik                                                     */
/*                      Universitaet Heidelberg                                                                                 */
/*                      Im Neuenheimer Feld 294                                                                                 */
/*                      6900 Heidelberg                                                                                                 */
/*                                                                                                                                                      */
/* History: 10.9.1993 begin                                                                                             */
/*                                                                                                                                                      */
/****************************************************************************/


/****************************************************************************/
/*                                                                                                                                                      */
/* auto include mechanism and other include files                                                       */
/*                                                                                                                                                      */
/****************************************************************************/

#ifndef __RULEGEN__
#define __RULEGEN__

#include <cstdio>

/****************************************************************************/
/*                                                                                                                                                      */
/* defines in the following order                                                                                       */
/*                                                                                                                                                      */
/*                compile time constants defining static data size (i.e. arrays)        */
/*                other constants                                                                                                       */
/*                macros                                                                                                                        */
/*                                                                                                                                                      */
/****************************************************************************/

/* macros for referencing of sons */
/* 4 high bits for no of neighbours to be passed */
#define PATHDEPTHMASK 0xF0000000
#define PATHDEPTHSHIFT 28
#define PATHDEPTH(i)                            (((i) & PATHDEPTHMASK)>>PATHDEPTHSHIFT)
#define SETPATHDEPTH(i,val)             (i) = ((i)&(~PATHDEPTHMASK))|(((val)<<PATHDEPTHSHIFT)&PATHDEPTHMASK)

/* 2 bits at position n for element side */
#define NEXTSIDEMASK 0x00000003
#define NEXTSIDE(i,n)                           (((i) & (NEXTSIDEMASK<<(2*(n))))>>(2*(n)))
#define SETNEXTSIDE(i,n,val)            (i) = ((i)&(~(NEXTSIDEMASK<<(2*(n)))))|(((val)&NEXTSIDEMASK)<<(2*(n)))


/****************************************************************************/
/*                                                                                                                                                      */
/* data structures exported by the corresponding source file                            */
/*                                                                                                                                                      */
/****************************************************************************/

/****************************************************************************/
/*                                                                                                                                                      */
/* definition of exported global variables                                                                      */
/*                                                                                                                                                      */
/****************************************************************************/

/****************************************************************************/
/*                                                                                                                                                      */
/* function declarations                                                                                                        */
/*                                                                                                                                                      */
/****************************************************************************/


#endif
