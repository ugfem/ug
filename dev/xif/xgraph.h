// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  xgraph.h														*/
/*																			*/
/* Purpose:   header file for graph window functionality					*/
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

#ifndef __XGRAPH__
#define __XGRAPH__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define CONTROLSIZE     24                              /* total size is 32 including line	*/

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef struct graphwindow {

  /* this is a linked list */
  struct graphwindow *next;                                     /* pointer to next window		*/

  /* drawing parameters */
  short marker_size;                                                    /* size of markers in pixels	*/
  short marker_id;                                                      /* number of marker                     */

  /* font metrics for fast access */
  int font_ascent;                                                      /* font ascent					*/
  int font_height;                                                      /* ascent+descent				*/
  int font_width;                                                       /* widest character                     */

  /* window size and position to filter resize & drag events */
  int window_x;
  int window_y;
  int window_width;
  int window_height;

  /* windows current point */
  int x;
  int y;

  /* X things */
  Window win;                                                           /* window id					*/
  GC gc;                                                                        /* a graphics context			*/
  XFontStruct *font_info;                                       /* the font structure			*/
  Region region;                                                        /* accumulate clipping region	*/
  char font_name[128];                                          /* font name from resource		*/
  Pixmap icon_pixmap;                                           /* icon to use					*/
  XTextProperty icon_name;                                      /* icons name					*/
  XTextProperty window_name;                                    /* windows name                                 */
} GraphWindow ;


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

OUTPUTDEVICE    *InitXOutputDevice      (void);
unsigned long   UGBlack                         (void);
unsigned long   UGWhite                         (void);
int                     DrawControls            (GraphWindow *gw);
int                     InitControls            (Window win);
int                     WhichTool                       (GraphWindow *gwin, int x, int y, int *tool);
int                     DrawRegion                      (GraphWindow *gwin, int x, int y);
GraphWindow     *WhichGW                        (Window win);
void                    SetCurrentGW            (GraphWindow *g);

#endif
