// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                                                                                                      */
/* File:          cw.h                                                                                                                  */
/*                                                                                                                                                      */
/* Purpose:   control word definitions header file                                                      */
/*                                                                                                                                                      */
/* Author:        Peter Bastian                                                                                                 */
/*                        Institut fuer Computeranwendungen III                                                 */
/*                        Universitaet Stuttgart                                                                                */
/*                        Pfaffenwaldring 27                                                                                    */
/*                        70569 Stuttgart                                                                                               */
/*                        email: ug@ica3.uni-stuttgart.de                                                           */
/*                                                                                                                                                      */
/* History:   11.01.95 begin, ug version 3.0                                                            */
/*                                                                                                                                                      */
/* Remarks:                                                                                                                             */
/*                                                                                                                                                      */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                                                                                                      */
/* auto include mechanism and other include files                                                       */
/*                                                                                                                                                      */
/****************************************************************************/

#ifndef __CW__
#define __CW__

/**************************************************/
/* A namespace for the c++ version                */
/**************************************************/
#ifdef __cplusplus
#ifdef __TWODIM__
namespace UG2d {
#else
namespace UG3d {
#endif
#endif

/** \brief Description of a control word */
typedef struct {

  /** \brief this struct is used */
  INT used;

  /** \brief name string */
  char *name;

  /** \brief where in object is it ? */
  unsigned INT offset_in_object;

  /** \brief bitwise object ID */
  INT objt_used;

  /** \brief used bits */
  unsigned INT used_mask;

} CONTROL_WORD;

/** \brief Manage part of a control word */
typedef struct {

  /** \brief this struct is used                          */
  INT used;

  /** \brief name string */
  char *name;

  /** \brief pointer to corresponding control word */
  INT control_word;

  /** \brief shift in control word */
  INT offset_in_word;

  /** \brief number of bits used */
  INT length;

  /** \brief bitwise object ID  */
  INT objt_used;

  /** \brief copy from control word (faster)      */
  unsigned INT offset_in_object;

  /** \brief 1 where bits are used                        */
  unsigned INT mask;

  /** \brief 0 where bits are used                        */
  unsigned INT xor_mask;

} CONTROL_ENTRY;

/** \brief Global array with descriptions */
extern CONTROL_WORD
  control_words[];

/** \brief Predefined control words */
extern CONTROL_ENTRY
  control_entries[MAX_CONTROL_ENTRIES];

/****************************************************************************/
/*                                                                                                                                                      */
/* function definitions                                                                                                         */
/*                                                                                                                                                      */
/****************************************************************************/

INT InitCW      (void);

#ifdef __cplusplus
}  /* namespace UG{2|3}d */
#endif

#endif
