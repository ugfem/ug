// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  uginterface.h                                                                                                 */
/*																			*/
/* Purpose:   defines data structure for ug interface						*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*																			*/
/* History:   14.06.93 begin, ug version ug21Xmas3d                                             */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UGINTERFACE__
#define __UGINTERFACE__


#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __WPM__
#include "wpm.h"
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#ifdef __TWODIM__
#define ALLOWED_TOOL(p,t)               (1)
#define REFRESH_TOOL(p,t)               (0)
#endif

#ifdef __THREEDIM__
#define ALLOWED_TOOL(p,t)               ((t)==arrowTool || ((t)==gnoedelTool && (p)->PlotCommandPtr!=NULL))
#define REFRESH_TOOL(p,t)               ((t)==gnoedelTool && (p)->PlotCommandPtr!=NULL)
#endif

#define MAXLEN_INFOSTRING               20

/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

INT                     InitUgInterface                                 ();

INT                     DelCmdKey                                               (char c);
INT                     SetCmdKey                                               (char c, const char *String);
INT                     DelAllCmdKeys                                   (void);
INT                     ListCmdKeys                                     ();

INT                     SetCurrentPicture                               (PICTURE *thePicture);
PICTURE                 *GetCurrentPicture                              (void);
INT                     SetCurrentUgWindow                              (UGWINDOW *theUgWindow);
UGWINDOW                *GetCurrentUgWindow                     (void);

INT                     UserInterrupt                                   (const char *text);
INT                     UserIn                                                  (char *String);
INT                     UserRead                                                (char *String);

INT                     SetRefreshState                                 (INT status);

#endif
