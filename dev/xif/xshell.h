// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  xshell.h														*/
/*																			*/
/* Purpose:   header file for shell window functionality					*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					        */
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

#ifndef __XSHELL__
#define __XSHELL__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* maximum size of text buffer */
#define MAXLINES                100
#define MAXCOLS                 130

#define SHELLICONNAME   "ug3.2 shell"
#define SHELLWINNAME    "ug3.2 shell"
#define RESOURCENAME    "ug3"

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef struct {
  /* the text buffer */
  char *lineStart[MAXLINES];                                    /* pointers to start of lines	*/
  int topLine;                                                          /* index of top line on screen	*/
  int numLines;                                                         /* current size of shell in x	*/
  int numCols;                                                          /* current size of shell in y	*/
  int line;                                                                     /* current line position		*/
  int col;                                                                      /* current col position                 */

  /* font metrics for fast access */
  int font_ascent;                                                      /* font ascent					*/
  int font_height;                                                      /* ascent+descent				*/
  int font_width;                                                       /* widest character                     */

  /* window size to filter resize event */
  int window_width;
  int window_height;

  /* X things */
        #ifdef USE_XAW
  Widget wid;                                                                   /* widget id					*/
        #endif /* USE_XAW */
  Window win;                                                           /* window id					*/
  GC gc;                                                                        /* a graphics context			*/
  XFontStruct *font_info;                                       /* the font structure			*/
  Region region;                                                        /* accumulate clipping region	*/
  char font_name[128];                                          /* font name from resource		*/
  Pixmap icon_pixmap;                                           /* icon to use					*/
  XTextProperty icon_name;                                       /* icons name					 */
  XTextProperty window_name;                                    /* windows name                                 */
} ShellWindow ;


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

int   ShellOpen                           (ShellWindow *sh);
void  ShellInsertChar             (ShellWindow *sh, char c);
void  ShellInsertString           (ShellWindow *sh, const char *s);
void  ShellHandleExposeEvent  (ShellWindow *sh, XEvent *report);
void  ShellHandleResizeEvent  (ShellWindow *sh, XEvent *report);
char *ShellHandleKeybordEvent (ShellWindow *sh, XEvent *report,int *cmdKey, int onlyCmdKey);
void AppendOrInsertCutbuffer  (ShellWindow *sh, char *cutbuf, int cnt);

#endif
