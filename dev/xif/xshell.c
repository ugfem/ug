// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  xshell.c														*/
/*																			*/
/* Purpose:   shell window for ug3 based on Xlib or X Athena Widget set     */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de		                                */
/*																			*/
/* History:   16.02.94 begin, ug version 3.0								*/
/*            06.06.95 implementation of Athena Text widget as Shellwindow      */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* X11 includes */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/Xresource.h>
#include <X11/keysym.h>

#define USE_XAW     /* enables use of Athena Text Widget  (libXaw, libXt) */
/* undefine USE_XAW, if you want pure X (only libX11) */

#ifdef USE_XAW
/* Xt & Xaw includes */
#include <X11/Intrinsic.h>
#include <X11/Vendor.h>
#include <X11/StringDefs.h>
#include <X11/Xaw/XawInit.h>
#include <X11/Xaw/AsciiText.h>
#include <X11/Xaw/Text.h>
#endif /* USE_XAW */

/* standard C includes */
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

/* interface includes */
#include "compiler.h"
#include "devices.h"
#include "initdev.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"

/* Xif includes */
#include "xmain.h"
#include "xshell.h"
#include "xgraph.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* hardcoded default values */
#define DEFAULTXPOS     50
#define DEFAULTYPOS     100
#define DEFAULTLINES    24
#define DEFAULTCOLS     80
#define DEFAULTFONT     "9x15"
#define DEFAULTMINX     300
#define DEFAULTMINY     200
#define DEFAULTBORDER   4

/* bitmaps for icons */
#include "shell-icon"

/* macros for access to text buffer */
#define MIN(x,y)                (((x)<(y)) ? (x) : (y))
#define MAX(x,y)                (((x)>(y)) ? (x) : (y))
#define TEXTLINE(sh,i)          (sh)->lineStart[((sh)->topLine+(i))%((sh)->numLines)]
#define XPOS(sh,i,j)    XTextWidth((sh)->font_info,TEXTLINE(sh,j),MIN(strlen(TEXTLINE(sh,j)),i))
#define YPOS(sh,i,j)    (sh)->font_ascent+(j)*((sh)->font_height)

/* change color */
#define BLACK(sh)               XSetForeground(display,(sh)->gc,UGBlack())
#define WHITE(sh)               XSetForeground(display,(sh)->gc,UGWhite())

#ifdef USE_XAW
#define HISTLEN         30
#define HISTNEXT(pos)   (((pos)+INPUTBUFFERLEN)%(HISTLEN*INPUTBUFFERLEN))
#define HISTPREV(pos)   ((pos==0) ? ((HISTLEN-1)*INPUTBUFFERLEN) : ((pos)-INPUTBUFFERLEN)%(HISTLEN*INPUTBUFFERLEN))
#endif

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

#ifdef USE_XAW
extern Widget toplevel;
extern Widget ugshell;
extern XtAppContext context;
#endif /* USE_XAW */

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static char buf[INPUTBUFFERLEN];                        /* input line buffer			*/
static int nchars=0;                                            /* number of charcters in input */
static char cmdKeyBuf;                                          /* buffer for a command key     */
static int AltFlag;                                             /* 1 if Alt has been pressed	*/
static int InsertPos=0;                                         /* only used with Xaw			*/

#ifdef USE_XAW
static int HistLastPos=0;
static int HistInsertPos=0;
static int HistReadPos=0;
static char HistBuf[HISTLEN*INPUTBUFFERLEN];
static int MaxLines=400;
static int s=0;
static int      *LineBuffer;
static int LineFirst=0;
static int LineInsert=0;
#endif

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

#ifdef USE_XAW
void UpAction (Widget w, XEvent *event, String *params, Cardinal *num_params);
void DownAction (Widget w, XEvent *event, String *params, Cardinal *num_params);
void LeftAction (Widget w, XEvent *event, String *params, Cardinal *num_params);
void RightAction (Widget w, XEvent *event, String *params, Cardinal *num_params);
#endif


/****************************************************************************/
/*																			*/
/* Function:   ShellShowCursor, ShellHideCursor                                                         */
/*																			*/
/* Purpose:    show/hide cursor at current point							*/
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellShowCursor (ShellWindow *sh)
{
  int x,y;

  /* get current position */
  x = XPOS(sh,sh->col,sh->line);
  y = YPOS(sh,sh->col,sh->line);

  BLACK(sh);
  XFillRectangle(display,sh->win,sh->gc,x,y-sh->font_ascent,sh->font_width,sh->font_height);
}

void ShellHideCursor (ShellWindow *sh)
{
  int x,y;

  /* get current position */
  x = XPOS(sh,sh->col,sh->line);
  y = YPOS(sh,sh->col,sh->line);

  WHITE(sh);
  XFillRectangle(display,sh->win,sh->gc,x,y-sh->font_ascent,sh->font_width,sh->font_height);
}


/****************************************************************************/
/*																			*/
/* Function:   ShellOpen													*/
/*																			*/
/* Purpose:    Open a shell window and initialize a ShellWindow data str.	*/
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*																			*/
/* Output:	   0: OK														*/
/*			   1: error, could not complete                                                                 */
/*																			*/
/****************************************************************************/

int ShellOpen (ShellWindow *sh)
{
  unsigned int width, height;
  int x, y, i;
  unsigned int border_width = DEFAULTBORDER;
  char *window_name = SHELLWINNAME;
  char *icon_name = SHELLICONNAME;
  char *Geometry, *Fontname;
  char *ug_name = RESOURCENAME;
  int flags;
  unsigned long valuemask = 0;
  static char dash_list[] = {12,24};
  XGCValues values;
  XSizeHints size_hints;
  XWMHints wm_hints;
  XClassHint class_hints;
        #ifdef USE_XAW
  char buffer[128];
  int n;
  Arg args[20];
  XtTranslations NewTranslations;
  XtActionsRec actions[] = {
    {"UpAction", UpAction},
    {"DownAction", DownAction},
    {"LeftAction", LeftAction},
    {"RightAction",RightAction},
  };
        #endif /* USE_XAW */

  /* allocate text buffer see, to know that memory is sufficient */
  sh->topLine = 0;
  sh->line = 0;
  sh->col = 0;
  for (i=0; i<MAXLINES; i++)
  {
    /* malloc a buffer line */
    sh->lineStart[i] = (char *)calloc(MAXCOLS+1,sizeof(char));
    if (sh->lineStart[i] == NULL)
    {
      fprintf(stderr,"%s could not allocate line %d\n",prog_name,i);
      return(1);
    }
  }

  /* find out geometry */
  width = DEFAULTCOLS;
  height = DEFAULTLINES;
  x = DEFAULTXPOS;
  y = DEFAULTYPOS;
  Geometry = XGetDefault(display,ug_name,"shellgeom");
  if (Geometry!=NULL)
  {
    flags = XParseGeometry(Geometry,&x,&y,&width,&height);
    if (!(WidthValue & flags)) width = DEFAULTCOLS;
    if (!(HeightValue & flags)) height = DEFAULTLINES;
    if (XValue & flags)
      if (XNegative&flags) x = display_width+x;
    if (YValue & flags)
      if (YNegative&flags) y = display_height+y;
  }
  width = MIN(width,MAXCOLS);       /* width & height in characters! */
  height = MIN(height,MAXLINES);
  sh->numLines = height;
  sh->numCols = width;

  /* first load font to know size of window */
  Fontname = XGetDefault(display,ug_name,"shellfont");
  if (Fontname==NULL)
  {
    printf("font 'shellfont' not found in Xdefaults\n"
           "using default "DEFAULTFONT);
    strcpy(sh->font_name,DEFAULTFONT);
  }
  else
    strcpy(sh->font_name,Fontname);
  if ( (sh->font_info=XLoadQueryFont(display,sh->font_name))==NULL)
  {
    /* font could not be loaded, try default font */
    strcpy(sh->font_name,DEFAULTFONT);
    if ( (sh->font_info=XLoadQueryFont(display,sh->font_name))==NULL)
    {
      fprintf(stderr,"%s could not load font %s\n",prog_name,sh->font_name);
      return(1);
    }
  }
  sh->font_ascent = sh->font_info->ascent;
  sh->font_height = sh->font_info->ascent+sh->font_info->descent;
  sh->font_width = sh->font_info->max_bounds.width;

  width = width*sh->font_width;
  height = height*sh->font_height;
  sh->window_width = width;
  sh->window_height = height;

        #ifdef USE_XAW
  n=0;
  XtSetArg(args[n], XtNx, x); n++;
  XtSetArg(args[n], XtNy, y); n++;
  XtSetArg(args[n], XtNwidth, width+sh->font_width); n++;
  XtSetArg(args[n], XtNheight, height+sh->font_height+sh->font_info->descent); n++;
  XtSetArg(args[n], XtNborderWidth, border_width); n++;
  XtSetArg(args[n], XtNscrollVertical, XawtextScrollWhenNeeded); n++;
  XtSetArg(args[n], XtNscrollHorizontal, XawtextScrollWhenNeeded); n++;
  XtSetArg(args[n], XtNstring, sh->lineStart[0]); n++;
  XtSetArg(args[n], XtNlength, MAXLINES*(MAXCOLS+1)); n++;
  XtSetArg(args[n], XtNpieceSize, MAXLINES*(MAXCOLS+1)); n++;
  XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
  /* TODO: delete this, line wrapping does not work correctly
     XtSetArg(args[n], XtNwrap,XawtextWrapLine); n++; */
  XtSetArg(args[n], XtNfont,sh->font_info); n++;
  XtSetArg(args[n], XtNbackground,WhitePixel(display,screen_num)); n++;
  XtSetArg(args[n], XtNforeground,BlackPixel(display,screen_num)); n++;

  ugshell = XtCreateManagedWidget ("ugshell", asciiTextWidgetClass,
                                   toplevel, args, n);

  /* modify Translation table */
  NewTranslations = XtParseTranslationTable
                      ("<Key>BackSpace: no-op(RingBell)\n\
						  <Key>Delete:    no-op(RingBell)\n\
						  <Key>Up:        no-op(RingBell)\n\
						  <Key>Down:      no-op(RingBell)\n\
						  <Key>Left:      no-op(RingBell)\n\
						  <Key>Right:     no-op(RingBell)\n");
  /*
     Ctrl<Key>/:     search(forward,"")\n\
     Ctrl Shift<Key>/: search(backward,"")\n");
   */
  XtOverrideTranslations(ugshell,NewTranslations);

  /* add special actions */
  XtAppAddActions(context,actions,XtNumber(actions));

  /* realize widget tree */
  XtRealizeWidget (toplevel);

  sh->wid = ugshell;
  sh->win = XtWindow(ugshell);

        #else /* USE_XAW */

  /* open a window */
  sh->win = XCreateSimpleWindow(display,RootWindow(display,screen_num),
                                x,y,width,height,border_width,BlackPixel(display,screen_num),
                                WhitePixel(display,screen_num));
        #endif /* USE_XAW */

  /* generate icon, needed for properties */
  sh->icon_pixmap = XCreateBitmapFromData(display,sh->win,shell_icon_bits,shell_icon_width,
                                          shell_icon_height);

  /* set standard properties */
  size_hints.flags = USPosition | USSize | PMinSize;
  size_hints.min_width = DEFAULTMINX;
  size_hints.min_height = DEFAULTMINY;
  if (XStringListToTextProperty(&window_name,1,&(sh->window_name))==0)
  {
    fprintf(stderr,"%s: structure alloc for window_name failed.\n",prog_name);
    exit(-1);
  }
  if (XStringListToTextProperty(&icon_name,1,&(sh->icon_name))==0)
  {
    fprintf(stderr,"%s: structure alloc for icon_name failed.\n",prog_name);
    exit(-1);
  }
  wm_hints.initial_state = NormalState;
  wm_hints.input = True;
  wm_hints.icon_pixmap = sh->icon_pixmap;
  wm_hints.flags = StateHint | IconPixmapHint | InputHint;
  class_hints.res_name = prog_name;
  class_hints.res_class = ug_name;
  XSetWMProperties(display,sh->win,&sh->window_name,&sh->icon_name,if_argv,if_argc,
                   &size_hints,&wm_hints,&class_hints);

        #ifndef USE_XAW
  /* select event types that will be received */
  XSelectInput(display,sh->win,ExposureMask|KeyPressMask|ButtonPressMask|StructureNotifyMask);
        #endif

  /* prepare graphics context */
  sh->gc = XCreateGC(display,sh->win,valuemask,&values);
  XSetFont(display,sh->gc,sh->font_info->fid);
  XSetForeground(display,sh->gc,BlackPixel(display,screen_num));
  XSetLineAttributes(display,sh->gc,0,LineSolid,CapRound,JoinRound);
  XSetDashes(display,sh->gc,0,dash_list,2);

  /* now map the window */
  XMapWindow(display,sh->win);

  /* create a region to accumulate update region */
  sh->region = XCreateRegion();

        #ifndef USE_XAW
  /* show cursor for first time */
  ShellShowCursor(sh);
        #endif

        #ifdef USE_XAW
  if (GetDefaultValue(DEFAULTSFILENAME,"shellbuffersize",buffer)==0)
  {
    sscanf(buffer," %d ",&n);
    if (n>0) MaxLines = n;
    else printf("shellbuffersize variable in defaults file = %d  - should be a positve integer (>0)\n",n);
  }
  LineBuffer = (int *) calloc(MaxLines,sizeof(int));
        #endif

  /* und tschuess */
  return(0);
}


/****************************************************************************/
/*																			*/
/* Function:   ShellScrollLine												*/
/*																			*/
/* Purpose:    scroll forward and backward by a number of lines                         */
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellScrollLine (ShellWindow *sh)
{
  /* copy area up one line */
  XCopyArea(display,sh->win,sh->win,sh->gc,0,sh->font_height,sh->window_width,
            sh->window_height-sh->font_height,0,0);

  /* delete last line */
  WHITE(sh);
  XFillRectangle(display,sh->win,sh->gc,0,sh->window_height-sh->font_height,
                 sh->window_width,sh->font_height);

  /* wrap text buffer */
  sh->topLine = (sh->topLine+1)%(sh->numLines);

  return;
}


/****************************************************************************/
/*																			*/
/* Function:   ShellInsertChar												*/
/*																			*/
/* Purpose:    Insert a character into shell window                                             */
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*			   char c														*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellInsertChar (ShellWindow *sh, char c)
{
  int x,y;
  int error;
        #ifdef USE_XAW
  Cardinal n;
  Arg args[10];
  XawTextPosition appendPos;
  XawTextBlock text;
  char t = c;
  char NewString[INPUTBUFFERLEN];

  text.firstPos = 0;
  text.length = 1;
  text.ptr = &t;
  text.format = FMT8BIT;

  appendPos = XawTextGetInsertionPoint(sh->wid);
  if (nchars > InsertPos)
  {
    n=0;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetValues(ugshell,args,n);
  }
  error = XawTextReplace(sh->wid,appendPos,appendPos,&text);
  if (nchars > InsertPos)
  {
    n=0;
    XtSetArg(args[n], XtNeditType, XawtextAppend); n++;
    XtSetValues(ugshell,args,n);
  }
  XawTextSetInsertionPoint(sh->wid,appendPos+1);
        #else /* USE_XAW */

  /* hide cursor */
  ShellHideCursor(sh);

  switch (c)
  {
  /* backspace */
  case 8 :
    /* decrement buffer position */
    if (sh->col>0)
      (sh->col)--;
    else
    {
      if (sh->line>0)
      {
        (sh->line)--;
        sh->col = strlen(TEXTLINE(sh,sh->line))-1;
      }
    }

    /* put 0 character at new position */
    TEXTLINE(sh,sh->line)[sh->col] = (char) 0;
    break;

  /* return */
  case '\n' :
    if (sh->line==sh->numLines-1)
      ShellScrollLine(sh);
    else
      (sh->line)++;
    sh->col = 0;
    TEXTLINE(sh,sh->line)[sh->col] = (char) 0;

    break;

  /* printable */
  default :
    if ((c<32)||(c>=127)) break;

    /* print character on screen */
    x = XPOS(sh,sh->col,sh->line);
    y = YPOS(sh,sh->col,sh->line);
    BLACK(sh);
                        #ifndef USE_XAW
    XDrawString(display,sh->win,sh->gc,x,y,&c,1);
                        #endif /* USE_XAW */
    /* insert character in text buffer */
    TEXTLINE(sh,sh->line)[sh->col] = c;

    /* place end of string after character */
    TEXTLINE(sh,sh->line)[sh->col+1] = (char) 0;

    /* advance cursor */
    (sh->col)++;
    if (sh->col>=sh->numCols) {
      if (sh->line==sh->numLines-1)
        ShellScrollLine(sh);
      else
        (sh->line)++;
      sh->col = 0;
    }
    IFDEBUG(dev,1)
    printf("char = %d x=%d y=%d line=%d col=%d bufferline=%d\n",(int)c,x,y,sh->line,sh->col,
           (sh->topLine+sh->line)%sh->numLines);
    ENDDEBUG
    break;
  }

  /* show cursor again */
  ShellShowCursor(sh);
        #endif /* USE_XAW */
}

/************************************************************************/
/*									*/
/* Function:	ShellUpdateTextBuffer					*/
/*									*/
/* Purpose:     Update the buffer of a shell window			*/
/*									*/
/* Input:	ShellWindow *sh						*/
/*		number of characters to remove   int i                  */
/*									*/
/* Output:	   none                                                 */
/*									*/
/************************************************************************/


#ifdef USE_XAW
void ShellUpdateTextBuffer (ShellWindow *sh,int i)
{
  XawTextPosition appendPos;
  XawTextBlock textInsert;

  textInsert.firstPos = 0;
  textInsert.length = 0;
  textInsert.ptr = NULL;
  textInsert.format = FMT8BIT;

  LineInsert=LineInsert%MaxLines;
  LineBuffer[LineInsert]=i;
  LineInsert=(LineInsert+1)%MaxLines;

  if(LineInsert==LineFirst)
  {
    appendPos=XawTextGetInsertionPoint(sh->wid);
    LineFirst=(LineFirst+1)%MaxLines;
    XawTextReplace(sh->wid,0,LineBuffer[LineInsert]-s,&textInsert);
    s=LineBuffer[LineInsert];
    XawTextSetInsertionPoint(sh->wid,appendPos);
  }
}
#endif

/****************************************************************************/
/*																			*/
/* Function:   ShellInsertString											*/
/*																			*/
/* Purpose:    Insert a string into shell window							*/
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*			   char *s														*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellInsertString (ShellWindow *sh, const char *s)
{
    #ifdef USE_XAW
  XawTextPosition appendPos;
  XawTextBlock text;

  text.firstPos = 0;
  text.length = strlen(s);
  text.ptr = (char *)s;
  text.format = FMT8BIT;

  appendPos=XawTextGetInsertionPoint(sh->wid);
  if( text.length!=0 && s[text.length-1]=='\n' )
    ShellUpdateTextBuffer(sh,appendPos+text.length);
  XawTextSetInsertionPoint(sh->wid,appendPos);
  XawTextReplace(sh->wid,appendPos,appendPos+text.length,&text);
  XawTextSetInsertionPoint(sh->wid,appendPos+text.length);
    #else
  int i;
  for (i=0; i<strlen(s); i++) ShellInsertChar(sh,s[i]);
    #endif
}


/****************************************************************************/
/*																			*/
/* Function:   ShellRefresh                                                                                             */
/*																			*/
/* Purpose:    refresh contents of shell window                                                         */
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellRefresh (ShellWindow *sh)
{
  int i,x,y,l;
  XRectangle rectangle;

  /* set clipping region to accumulated region in shell */
  XSetRegion(display,sh->gc,sh->region);

  /* draw characters */
  BLACK(sh);
  for (i=0; i<sh->numLines; i++)
    if ((l=strlen(TEXTLINE(sh,i)))>0)
    {
      x = 0;
      y = YPOS(sh,0,i);
      XDrawString(display,sh->win,sh->gc,x,y,TEXTLINE(sh,i),l);
    }

  /* show cursor */
  ShellShowCursor(sh);

  /* set clipping region to whole window for subsequent drawings */
  rectangle.x = 0;
  rectangle.y = 0;
  rectangle.width = (short) sh->window_width;
  rectangle.height = (short) sh->window_height;
  XUnionRectWithRegion(&rectangle,sh->region,sh->region);
  XSetRegion(display,sh->gc,sh->region);

  /* clear clipping region */
  XDestroyRegion(sh->region);
  sh->region = XCreateRegion();

}


/****************************************************************************/
/*																			*/
/* Function:   ShellResize													*/
/*																			*/
/* Purpose:    resize the shell window										*/
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*			   int width		 new size of window in pixels !                         */
/*			   int height													*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellResize (ShellWindow *sh, int width, int height)
{
  int i;
  int c_width,c_height;
  int obl,ntl;

  /* compute new size in characters */
  c_width = MIN(width/sh->font_width,MAXCOLS);
  c_height = MIN(height/sh->font_height,MAXLINES);
  IFDEBUG(dev,1)
  printf("RESIZE: OLD: top=%2d nl=%2d line=%2d    ",sh->topLine,sh->numLines,sh->line);
  ENDDEBUG
  /* clear lines not in use */
  for (i=sh->numLines; i<MAXLINES; i++) sh->lineStart[i][0] = (char) 0;

  /* adjust buffer vertically */
  if (c_height!=sh->numLines)
  {
    obl = (sh->topLine+sh->numLines-1)%sh->numLines;
    if (c_height>sh->numLines)
    {
      ntl = obl+1;
      for (i=0; i<sh->numLines-ntl; i++)
        strcpy(sh->lineStart[c_height-1-i],sh->lineStart[sh->numLines-1-i]);
      for (i=0; i<sh->numLines-ntl; i++)
        sh->lineStart[ntl+i][0] = (char) 0;
      sh->line += c_height-sh->numLines;
    }
    else
    {
      if (obl<=c_height-1)
      {
        ntl = (obl+1)%c_height;
        for (i=0; i<c_height-1-obl; i++)
          strcpy(sh->lineStart[(obl+1)+i],sh->lineStart[(sh->numLines-c_height)+(obl+1)+i]);

      }
      else
      {
        ntl = 0;
        for (i=0; i<=c_height-1; i++)
          strcpy(sh->lineStart[i],sh->lineStart[obl-c_height+1+i]);
      }
      sh->line = MAX(0,sh->line-(sh->numLines-c_height));
    }
    sh->topLine = ntl;
    sh->numLines = c_height;
  }

  /* adjust buffer horizontally */
  for (i=0; i<MAXLINES; i++) sh->lineStart[i][c_width] = (char) 0;
  sh->numCols = c_width;
  sh->col = MIN(sh->col,c_width-1);
  IFDEBUG(dev,1)
  printf("NEW: top=%2d nl=%2d line=%2d\n",sh->topLine,sh->numLines,sh->line);
  ENDDEBUG
}


/****************************************************************************/
/*																			*/
/* Function:   ShellHandleExposeEvent										*/
/*																			*/
/* Purpose:    handle expose event for shell window                                             */
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*			   XEvent *report												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellHandleExposeEvent (ShellWindow *sh, XEvent *report)
{
  XRectangle rectangle;

  rectangle.x = (short) report->xexpose.x;
  rectangle.y = (short) report->xexpose.y;
  rectangle.width = (short) report->xexpose.width;
  rectangle.height = (short) report->xexpose.height;
  XUnionRectWithRegion(&rectangle,sh->region,sh->region);
  IFDEBUG(dev,1)
  printf("EXPOSE EVENT: rect= %d %d %d %d cnt= %d\n",rectangle.x,rectangle.y,
         rectangle.width,rectangle.height,report->xexpose.count);
  ENDDEBUG
  if (report->xexpose.count!=0) return;
  ShellRefresh(sh);
}


/****************************************************************************/
/*																			*/
/* Function:   ShellHandleResizeEvent										*/
/*																			*/
/* Purpose:    handle resize event for shell window                                             */
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*			   XEvent *report												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellHandleResizeEvent (ShellWindow *sh, XEvent *report)
{
  XRectangle rectangle;

  if (  (report->xconfigure.width==sh->window_width)
        &&(report->xconfigure.height==sh->window_height)) return;
  rectangle.x = 0;
  rectangle.y = 0;
  rectangle.width = (short) report->xconfigure.width;
  rectangle.height = (short) report->xconfigure.height;
  XUnionRectWithRegion(&rectangle,sh->region,sh->region);
  WHITE(sh);
  XFillRectangle(display,sh->win,sh->gc,0,0,report->xconfigure.width,
                 report->xconfigure.height);
  sh->window_width = report->xconfigure.width;
  sh->window_height = report->xconfigure.height;

  ShellResize(sh,report->xconfigure.width,report->xconfigure.height);
}


/****************************************************************************/
/*																			*/
/* Function:   ShellHandleKeybordEvent										*/
/*																			*/
/* Purpose:    handle keypress event for shell window						*/
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*			   XEvent *report												*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

char *ShellHandleKeybordEvent (ShellWindow *sh, XEvent *report,int *cmdKey, int onlyCmdKey)
{
  char buffer[32];
  int bufsize = 32;
  KeySym keysym;
  XKeyEvent keyevent;
  XComposeStatus compose;
  int count;
  char c;
        #ifdef USE_XAW
  Cardinal n;
  Arg args[10];
  XtTranslations NewTranslations;
        #endif

  /* no command key is default */
  *cmdKey = 0;
  keyevent = report->xkey;
  count = XLookupString(&keyevent,buffer,bufsize,&keysym,&compose);
  /* vorher XLookupString(report,buffer,bufsize,&keysym,&compose) */

  if((keysym==XK_Alt_L)||(keysym==XK_Alt_R)||
     (keysym==XK_Meta_L)||(keysym==XK_Meta_R))
  {
    AltFlag = 1;
    return(NULL);
  }
  else if (((keysym==XK_Return)||(keysym==XK_KP_Enter)||(keysym==XK_Linefeed)) && !onlyCmdKey)
  {
    c = '\n';
                #ifndef USE_XAW
    ShellInsertChar(sh,c);
                #endif
    buf[nchars] = (char) 0;
                #ifdef USE_XAW
    NewTranslations = XtParseTranslationTable
                        ("<Key>BackSpace: no-op(RingBell)\n\
							  <Key>Delete:    no-op(RingBell)\n\
							  <Key>Left:	  no-op(RingBell)\n");
    XtOverrideTranslations(ugshell,NewTranslations);
    if (nchars>0)
    {
      memcpy(&HistBuf[HistInsertPos],buf,nchars);
      NewTranslations = XtParseTranslationTable
                          ("<Key>Up:	UpAction()\n\
								  <Key>Down:no-op(RingBell)\n");
      XtOverrideTranslations(ugshell,NewTranslations);
      HistInsertPos = HISTNEXT(HistInsertPos);
      HistBuf[HistInsertPos] = (char) 0;
      if (HistInsertPos == HistLastPos)
        HistLastPos = HISTNEXT(HistLastPos);
      HistReadPos = HistInsertPos;
    }

    XawTextSetInsertionPoint(ugshell,XawTextGetInsertionPoint(ugshell)-InsertPos+nchars);
                #endif /* USE_XAW */

    nchars = 0;
    InsertPos = nchars;
    AltFlag = 0;
    return(buf);
  }
  else if ((keysym==XK_BackSpace)||(keysym==XK_Delete))
  {
    if (nchars>0)
    {
                        #ifdef USE_XAW
      if (nchars==InsertPos)
      {
                        #endif
      c = (char)8;
                                #ifndef USE_XAW
      ShellInsertChar(sh,c);
                                #endif
      nchars--;
      InsertPos = nchars;
                        #ifdef USE_XAW
    }
    else
    {
      memmove(&buf[InsertPos-1],&buf[InsertPos],nchars-InsertPos);
      nchars--;
      InsertPos--;
    }
                        #endif
    }
                #ifdef USE_XAW
    else
    {
      NewTranslations = XtParseTranslationTable
                          ("<Key>BackSpace: no-op(RingBell)\n\
								  <Key>Delete:    no-op(RingBell)\n\
								  <Key>Left:	  no-op(RingBell)\n");
      XtOverrideTranslations(ugshell,NewTranslations);
    }
                #endif


    AltFlag = 0;
    return(NULL);
  }

  else if (  ((keysym>=XK_KP_Space)&&(keysym<=XK_KP_9))
             ||((keysym>=XK_space)&&(keysym<=XK_asciitilde)))
  {

    /* vorher		if (report->xkey.state&Mod1Mask) */
    if(AltFlag)
    {
      cmdKeyBuf = buffer[0];
      *cmdKey = 1;
      AltFlag = 0;
      return(&cmdKeyBuf);
    }
    c = buffer[0];
                #ifndef USE_XAW
    ShellInsertChar(sh,c);
                #else
    if (nchars==0)
    {
      NewTranslations = XtParseTranslationTable
                          ("<Key>BackSpace: delete-previous-character() \n\
							  <Key>Delete:    delete-previous-character() \n\
							  <Key>Left:	  LeftAction() \n");
      XtOverrideTranslations(ugshell,NewTranslations);
    }
                #endif
                #ifdef USE_XAW
    if (InsertPos==nchars)
    {
                #endif
    buf[nchars] = c;
    nchars++;
    InsertPos = nchars;
                #ifdef USE_XAW
  }
  else
  {
    memmove(&buf[InsertPos+1],&buf[InsertPos],nchars-InsertPos+1);
    buf[InsertPos] = c;
    nchars++;
    InsertPos++;
  }
                #endif

    if (nchars==INPUTBUFFERLEN-1)
    {
      buf[nchars] = (char) 0;
                        #ifdef USE_XAW
      memcpy(&HistBuf[HistReadPos],buf,nchars);
      HistInsertPos = HISTNEXT(HistInsertPos);
      HistReadPos = HistInsertPos;
                        #endif
      nchars = 0;
      InsertPos = nchars;
      AltFlag = 0;
      return(buf);
    }
    AltFlag = 0;
    return(NULL);
  }

  AltFlag = 0;
  return(NULL);
}

#ifdef USE_XAW
/****************************************************************************/
/*																			*/
/* Function:   ShellReplaceString											*/
/*																			*/
/* Purpose:    Replace textposition 'from' to 'to' in shell window by           */
/*			   string s. The textposition is counted from the end.          */
/*																			*/
/* Input:	   ShellWindow *sh												*/
/*			   char *s														*/
/*			   int from			                                                                                */
/*			   int to                                                                                                       */
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ShellReplaceString (ShellWindow *sh, const char *s, int from, int to)
{
  XawTextPosition appendPos;
  XawTextBlock text;

  text.firstPos = 0;
  text.length = strlen(s);
  text.ptr = (char *)s;
  text.format = FMT8BIT;

  appendPos = XawTextGetInsertionPoint(ugshell);
  XawTextReplace(ugshell,appendPos-from,appendPos-to,&text);
  XawTextSetInsertionPoint(ugshell,appendPos+text.length);
}

/****************************************************************************/
/*																			*/
/* Function:   AppendOrInsertCutbuffer                                                                          */
/*																			*/
/* Purpose:    insert X cutbuffer to the shell cmdline buffer and			*/
/*			   display modified cmdline in ugshell window, if insertion		*/
/*																			*/
/* Input:	   char *cutbuf                                                                                         */
/*			   int  cnt                                                                                             */
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void AppendOrInsertCutbuffer (ShellWindow *sh, char *cutbuf, int cnt)
{
  int i;
  Cardinal n;
  XtTranslations NewTranslations;
  Arg args[10];

  if (nchars != InsertPos)
  {
    n=0;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetValues(ugshell,args,n);
  }

  if (cnt>0)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>BackSpace: delete-previous-character() \n\
							  <Key>Delete:    delete-previous-character() \n\
							  <Key>Left:      LeftAction() \n");
    XtOverrideTranslations(sh->wid,NewTranslations);
  }

  for (i=0; (i<cnt) && (nchars<INPUTBUFFERLEN-1); i++)
  {
    if (isgraph(cutbuf[i])||isspace(cutbuf[i]))
    {
      if (InsertPos==nchars)
      {
        buf[nchars] = cutbuf[i];
        nchars++;
        InsertPos = nchars;
      }
      else
      {
        memmove(&buf[InsertPos+1],&buf[InsertPos],nchars-InsertPos+1);
        buf[InsertPos] = cutbuf[i];
        nchars++;
        InsertPos++;
      }
    }
  }
}

/****************************************************************************/
/*																			*/
/* Function:   ReplaceCmd                                                                                       */
/*																			*/
/* Purpose:    replace the commandbuffer by a command from history			*/
/*																			*/
/* Input:	   Widget w: widget, which has gotten the action				*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void ReplaceCmd(Widget w)
{
  char *ActCommandPtr;
  int ActCommandLen;
  XawTextPosition ReplacePos;
  XawTextBlock text;
  XtTranslations NewTranslations;

  ActCommandPtr = &HistBuf[HistReadPos];
  ActCommandLen = strlen(ActCommandPtr);
  memcpy(buf,ActCommandPtr,ActCommandLen);

  text.firstPos = 0;
  text.length = ActCommandLen;
  text.ptr = ActCommandPtr;
  text.format = FMT8BIT;

  ReplacePos = XawTextGetInsertionPoint(w);
  XawTextReplace(w,ReplacePos-nchars,ReplacePos,&text);
  XawTextSetInsertionPoint(w,ReplacePos-nchars+ActCommandLen);

  nchars = ActCommandLen;
  InsertPos = nchars;
  if (nchars <=0)
  {
    NewTranslations = XtParseTranslationTable
                        ("<Key>BackSpace: no-op(RingBell) \n\
							  <Key>Delete:    no-op(RingBell) \n\
							  <Key>Left:      no-op(RingBell) \n\
							  <Key>Right:     no-op(RingBell) \n");
    XtOverrideTranslations(ugshell,NewTranslations);
  }
  else
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>BackSpace: delete-previous-character() \n\
							  <Key>Delete:    delete-previous-character() \n\
							  <Key>Left:      LeftAction() \n\
							  <Key>Right:     no-op(RingBell) \n");
    XtOverrideTranslations(ugshell,NewTranslations);
  }
}

/****************************************************************************/
/*																			*/
/* Function:   UpAction                                                                                         */
/*																			*/
/* Purpose:    display previous command from history, if exists                         */
/*																			*/
/* Input:	   Widget w: widget, which has gotten the action				*/
/*			   XEvent event, String *params, Cardinal *num_params:          */
/*			                dummy parameter to match funktionprototype				*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void UpAction (Widget w, XEvent *event, String *params, Cardinal *num_params)
{
  XtTranslations NewTranslations;

  if (HistReadPos == HistInsertPos)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Down:      DownAction()\n");
    XtOverrideTranslations(ugshell,NewTranslations);
  }

  HistReadPos = HISTPREV(HistReadPos);
  if (HistReadPos == HistLastPos)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Up:        no-op(RingBell)\n");
    XtOverrideTranslations(ugshell,NewTranslations);
  }

  ReplaceCmd(w);
}

/****************************************************************************/
/*																			*/
/* Function:   DownAction                                                                                       */
/*																			*/
/* Purpose:    display next command from history, if exists                             */
/*																			*/
/* Input:	   Widget w: widget, which has gotten the action				*/
/*			   XEvent event, String *params, Cardinal *num_params:          */
/*			                dummy parameter to match funktionprototype				*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void DownAction (Widget w, XEvent *event, String *params, Cardinal *num_params)
{
  XtTranslations NewTranslations;

  if (HistReadPos == HistLastPos)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Up:        UpAction()\n");
    XtOverrideTranslations(ugshell,NewTranslations);
  }
  HistReadPos = HISTNEXT(HistReadPos);
  if (HistReadPos == HistInsertPos)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Down:      no-op(RingBell)\n");

    XtOverrideTranslations(ugshell,NewTranslations);
  }
  ReplaceCmd(w);
}

/****************************************************************************/
/*																			*/
/* Function:   LeftAction                                                                                       */
/*																			*/
/* Purpose:    move cursor one position left in cmdline                     */
/*																			*/
/* Input:	   Widget w: widget, which has gotten the action				*/
/*			   XEvent event, String *params, Cardinal *num_params:          */
/*			                dummy parameter to match funktionprototype				*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void LeftAction (Widget w, XEvent *event, String *params, Cardinal *num_params)
{
  int NewPos;
  Cardinal n;
  Arg args[10];
  XtTranslations NewTranslations;

  if (InsertPos == nchars)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Right:      RightAction()\n");

    XtOverrideTranslations(ugshell,NewTranslations);
  }

  NewPos = XawTextGetInsertionPoint(ugshell);
  NewPos--;
  XawTextSetInsertionPoint(ugshell,NewPos);
  InsertPos--;

  if (InsertPos == 0)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Left:      no-op(RingBell)\n");

    XtOverrideTranslations(ugshell,NewTranslations);
  }
}

/****************************************************************************/
/*																			*/
/* Function:   RightAction                                                                                      */
/*																			*/
/* Purpose:    move cursor one position right in cmdline                    */
/*																			*/
/* Input:	   Widget w: widget, which has gotten the action				*/
/*			   XEvent event, String *params, Cardinal *num_params:          */
/*			                dummy parameter to match funktionprototype				*/
/*																			*/
/* Output:	   none                                                                                                                 */
/*																			*/
/****************************************************************************/

void RightAction (Widget w, XEvent *event, String *params, Cardinal *num_params)
{
  int NewPos;
  Cardinal n;
  Arg args[10];
  XtTranslations NewTranslations;

  if (InsertPos == 0)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Left:      LeftAction()\n");

    XtOverrideTranslations(ugshell,NewTranslations);
  }

  NewPos = XawTextGetInsertionPoint(ugshell);
  NewPos++;
  XawTextSetInsertionPoint(ugshell,NewPos);
  InsertPos++;

  if (InsertPos == nchars)
  {
    /* modify Translation table */
    NewTranslations = XtParseTranslationTable
                        ("<Key>Right:      no-op(RingBell)\n");

    XtOverrideTranslations(ugshell,NewTranslations);
  }
}

#endif /* USE_XAW */
