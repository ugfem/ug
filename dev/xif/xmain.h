// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  xmain.h														*/
/*																			*/
/* Purpose:   global variables exported by xmain							*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					*/
/*																			*/
/* History:   17.02.94 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __XMAIN__
#define __XMAIN__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/


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

extern Display *display;                                        /* the display					*/
extern int screen_num;                                          /* screen on display			*/
extern char *prog_name;                                         /* my own name					*/
extern Screen *screen_ptr;                                      /* dont know for what			*/
extern unsigned int display_width;                      /* size of screen if needed     */
extern unsigned int display_height;
extern int if_argc;                                             /* command line args			*/
extern char **if_argv;


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

#endif
