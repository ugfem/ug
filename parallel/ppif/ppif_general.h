// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ppif_general.h												*/
/*																			*/
/* Purpose:   ppif general header file                                                                  */
/*																			*/
/* Author:	  Stefan Lang                                                                           */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   24.06.96 start: PPIF_VERSION string, RCSID macro              */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __PPIF_GENERAL__
#define __PPIF_GENERAL__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* defintions of RCS macros */
#define PPIF_VERSION "PPIF_1_0"

#define PPIF_RCS_STRING PPIFRCSSTRING(PPIF_VERSION,ARCH_VERSION)
#define PPIFRCSSTRING(ppif_version,arch_version)\
  PPIFRCSSTRINGAUX(ppif_version,arch_version)
#define PPIFRCSSTRINGAUX(ppif_version,arch_version)\
  "$" "State: PPIF_VERSION=" # ppif_version " ARCH_VERSION=" # arch_version " $"

#define RCSID(header,module_rcs_string) RCSIDAUX(header,module_rcs_string)
#define RCSIDAUX(header,module_rcs_string) static char rcsid[] = header ## module_rcs_string;


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

#endif
