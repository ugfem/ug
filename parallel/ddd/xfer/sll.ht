/****************************************************************************/
/*                                                                          */
/* File:      sll.ht                                                        */
/*                                                                          */
/* Purpose:   single linked list templates, header file                     */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   960826 kb  created                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* macros for SLL data structure definitions                                */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/* linked list of items */
extern T *list(T);

/* number of items in list */
extern int n(T);


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* from sll.ct */
T *New(T) (SLLNewArgProtos);
T **SortedArray(T) (int (*) (const void *, const void *));
int Unify(T) (T **, int (*) (T **, T **));
#ifdef SLL_WithOrigOrder
void OrigOrder(T) (T **, int);
#endif
void Init(T) (void);
void FreeAll(T) (void);
void GetSizes(T) (int *, int *, size_t *, size_t *);



/****************************************************************************/

/* undefine switches and helper macros */

#ifdef SLL_WithOrigOrder
#undef SLL_WithOrigOrder
#endif

/****************************************************************************/

