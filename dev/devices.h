// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  devices.h                                                                                                     */
/*																			*/
/* Purpose:   implements a simple but portable graphical user interface         */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de					                */
/*																			*/
/* History:   14.06.93 begin, ug version ug21Xmas3d                                             */
/*			  16.12.94 restructured for ug version 3.0						*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __DEVICESH__
#define __DEVICESH__

#include <stdio.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

/****************************************************************************/
/*																			*/
/*	constants																*/
/*																			*/
/****************************************************************************/

/* markers																	*/
#define EMPTY_SQUARE_MARKER             0
#define GRAY_SQUARE_MARKER                      1
#define FILLED_SQUARE_MARKER            2
#define EMPTY_CIRCLE_MARKER             3
#define GRAY_CIRCLE_MARKER                      4
#define FILLED_CIRCLE_MARKER            5
#define EMPTY_RHOMBUS_MARKER            6
#define GRAY_RHOMBUS_MARKER             7
#define FILLED_RHOMBUS_MARKER           8
#define PLUS_MARKER                             9
#define CROSS_MARKER                            10

#define NMARKERS                                        11

/* tool numbers */
#define nboftools                               7
#define arrowTool                               0
#define crossTool                               1
#define choiceTool                              2
#define circleTool                              3
#define handTool                                4
#define heartTool                               5
#define gnoedelTool                     6

/* tool names */
#define arrowToolName                   "pointer tool"
#define crossToolName                   "insert bn tool"
#define choiceToolName                  "move nd tool"
#define circleToolName                  "insert in tool"
#define handToolName                    "select nd tool"
#define heartToolName                   "select el tool"
#define gnoedelToolName                 "mark red tool"

/* text position */
#define TEXT_NOT_CENTERED               0
#define TEXT_CENTERED                   1

/* text modes */
#define TEXT_REGULAR                    0
#define TEXT_INVERSE                    1
#define TEXT_INDEXED                    2

/* buffer space for command line input */
#define INPUTBUFFERLEN                  256

/* event types */
#define EVENT_ERROR                     0

#define EVERY_EVENT                     1

#define NO_EVENT                                2
#define TERM_GOAWAY                     3
#define TERM_CMDKEY                     4
#define TERM_STRING                     5
#define DOC_GOAWAY                              6
#define DOC_ACTIVATE                    7
#define DOC_DRAG                                8
#define DOC_GROW                                9
#define DOC_CHANGETOOL                  10
#define DOC_CONTENTCLICK                11
#define DOC_UPDATE                              12


/****************************************************************************/
/*																			*/
/*	macros																	*/
/*																			*/
/****************************************************************************/

#define EVENT_TYPE(p)                   ((p).Type)

/****************************************************************************/
/*																			*/
/*	data types																*/
/*																			*/
/****************************************************************************/

/* identification of windows */
typedef INT WINDOWID;

/* type for device coordinates */
typedef struct
{
  short x;
  short y;
} SHORT_POINT ;

/* function types exported by OUTPUTDEVICE */
typedef WINDOWID (*OpenOutputPtr)(const char *title, INT x, INT y, INT width, INT height, INT *Global_LL, INT *Global_UR, INT *Local_LL, INT *Local_UR, INT *error);
typedef INT (*CloseOutputPtr)(WINDOWID win);
typedef INT (*ActivateOutputPtr)(WINDOWID win);
typedef INT (*UpdateOutputPtr)(WINDOWID win, char *s, INT tool);

/* abstract graphical output device */
struct outputdevice {

  /* This is an environment variable */
  ENVVAR v;

  /* properties */
  long black;                                                   /* value for black										*/
  long white;                                                   /* value for white										*/
  long red;                                                             /* value for red										*/
  long green;                                                   /* value for green										*/
  long blue;                                                            /* value for blue										*/
  long cyan;                                                            /* value for cyan										*/
  long orange;                                                  /* value for orange                                                                     */
  long yellow;                                                  /* value for yellow                                                                     */
  long darkyellow;                                              /* value for yellow                                                                     */
  long magenta;                                                 /* value for magenta									*/
  short hasPalette;                                             /* 1 if device has a color lookup ta					*/
  long range;                                                   /* # of possible color indices							*/
  long spectrumStart;                                   /* usable range for a continuous						*/
  long spectrumEnd;                                             /* color spectrum										*/
  DOUBLE PixelRatio;                                            /* ratio of (physical) hight to width of a pixel		*/
  short signx;                                                  /* direction of increasing x-coordinates				*/
  short signy;                                                  /* direction of increasing y-coordinates				*/

  /* pointers to basic drawing functions */
  void (*Move)(SHORT_POINT);                                                            /* move in device coordinates		*/
  void (*Draw)(SHORT_POINT);                                                            /* draw from current point to given */
  void (*Polyline)(SHORT_POINT *, INT );                                        /* draw a polyline					*/
  void (*InversePolyline)(SHORT_POINT *, INT );                         /* draw an inverted polyline		*/
  void (*Polygon)(SHORT_POINT *, INT );                                         /* fill a polygon w. curr. col		*/
  void (*InversePolygon)(SHORT_POINT *, INT );                          /* invert a polygon w. curr. col	*/
  void (*ErasePolygon)(SHORT_POINT *, INT );                            /* erase a polygon w. curr. col         */
  void (*Polymark)(short, SHORT_POINT *);                                       /* place markers					*/
  void (*Text)(const char *, INT);                                                      /* draw text in current size/font	*/
  void (*CenteredText)(SHORT_POINT, const char *, INT);         /* draw text centered at x,y		*/
  void (*ClearViewPort)(void);                                                          /* clear a view port				*/

  /* pointers to set functions */
  void (*SetLineWidth)(short);                                                          /* line width in pixels (points)	*/
  void (*SetTextSize)(short);                                                           /* text size in pixels (points)         */
  void (*SetMarker)(short);                                                                     /* set marker id					*/
  void (*SetMarkerSize)(short);                                                         /* marker size in pixels (points)	*/
  void (*SetColor)(long);                                                                       /* arg is index or direct col value */
  void (*SetPaletteEntry)(long,short,short,short);                      /* set index to value				*/
  void (*SetNewPalette)(long,long,short*,short*,short*);        /* replace entrie					*/

  /* pointers to miscellaneous functions */
  void (*GetPaletteEntry)(long,short *,short *,short *);        /* read color tab					*/
  void (*Flush)(void);                                                                          /* flush graphics buffer			*/

  /* operations for managing windows */
  OpenOutputPtr OpenOutput;                       /* function to open a window								*/
  CloseOutputPtr CloseOutput;             /* function to close a window                                                         */
  ActivateOutputPtr ActivateOutput;       /* function to activate window							*/
  UpdateOutputPtr UpdateOutput;           /* function to draw outline with tool and info box		*/
};



/* event structure */
typedef struct {                                        /* no event                                                             */
  INT Type;                                                     /* event type								*/

  /* data */
  INT InterfaceEvent;                           /* 1 if the interface event was handled         */
} NO_UGEVENT;

typedef struct {                                        /* go away event for terminal window		*/
  INT Type;                                                     /* event type								*/
} TERM_GOAWAY_EVENT;

typedef struct {                                        /* cmd key event for terminal window		*/
  INT Type;                                                     /* event type								*/

  /* data */
  char CmdKey;                                          /* character from keyboard					*/
} TERM_CMDKEY_EVENT;

typedef struct {                                        /* string event for terminal window             */
  INT Type;                                                     /* event type								*/

  /* data */
  char String[INPUTBUFFERLEN];          /* string from keyboard                                         */
} TERM_STRING_EVENT;

typedef struct {                                        /* go away event for view					*/
  INT Type;                                                     /* event type								*/

  /* data */
  WINDOWID win;                                         /* the window								*/
} DOC_GOAWAY_EVENT;

typedef struct {                                        /* activate event for view					*/
  INT Type;                                                     /* event type								*/

  /* data */
  WINDOWID win;                                         /* the window								*/
} DOC_ACTIVATE_EVENT;

typedef struct {                                        /* drag event for view						*/
  INT Type;                                                     /* event type								*/

  /* data */
  WINDOWID win;                                         /* the window								*/
  INT Global_LL[2];                                     /* new absolute position of window on screen*/
  INT Global_UR[2];                                     /*											*/
} DOC_DRAG_EVENT;

typedef struct {                                        /* grow event for view						*/
  INT Type;                                                     /* event type								*/

  /* data */
  WINDOWID win;                                         /* the window								*/
  INT Global_LL[2];                                     /* new absolute position of window on screen*/
  INT Global_UR[2];                                     /*											*/
  INT Local_LL[2];                                      /* range of pixels used for plotting		*/
  INT Local_UR[2];                                      /*											*/
} DOC_GROW_EVENT;

typedef struct {                                        /* change tool event for view				*/
  INT Type;                                                     /* event type								*/

  /* data */
  WINDOWID win;                                         /* the window								*/
  INT Tool;                                                     /* change to that tool						*/
} DOC_CHANGETOOL_EVENT;

typedef struct {                                        /* content click event for view                         */
  INT Type;                                                     /* event type								*/

  /* data */
  WINDOWID win;                                         /* the window								*/
  INT MousePosition[2];                         /* mouse position							*/
} DOC_CONTENTCLICK_EVENT;

typedef struct {                                        /* update event for view					*/
  INT Type;                                                     /* event type								*/

  /* data */
  WINDOWID win;                                         /* the window								*/
} DOC_UPDATE_EVENT;

typedef union {
  INT Type;
  NO_UGEVENT NoEvent;
  TERM_GOAWAY_EVENT TermGoAway;
  TERM_CMDKEY_EVENT TermCmdKey;
  TERM_STRING_EVENT TermString;
  DOC_GOAWAY_EVENT DocGoAway;
  DOC_ACTIVATE_EVENT DocActivate;
  DOC_DRAG_EVENT DocDrag;
  DOC_GROW_EVENT DocGrow;
  DOC_CHANGETOOL_EVENT DocChangeTool;
  DOC_CONTENTCLICK_EVENT DocContentClick;
  DOC_UPDATE_EVENT DocUpdate;
} EVENT;


typedef struct outputdevice OUTPUTDEVICE;


/****************************************************************************/
/*																			*/
/* function exported by the devices module									*/
/*																			*/
/****************************************************************************/

/* initialization and clean up */
INT               InitDevices                           (int *argcp, char **argv);
void          ExitDevices               (void);

/* return the size of the monitor screen, return 0 if not available */
INT               GetScreenSize                         (INT size[2]);

/* set/get mute level for output control */
void              SetMuteLevel                          (INT mute);
INT               GetMuteLevel                          (void);

/* multiple output devices */
OUTPUTDEVICE *CreateOutputDevice                (char *name);
OUTPUTDEVICE *GetOutputDevice                   (const char *name);
OUTPUTDEVICE *GetDefaultOutputDevice    (void);

/* text output to shell with log file mechanism */
void              UserWrite                             (const char *s);
int               UserWriteF                            (const char *format, ...);
void              PrintErrorMessage             (char type, const char *procName, const char *text);
void              PrintErrorMessageF            (char type, const char *procName, const char *format, ...);
INT               OpenLogFile                           (const char *name);
INT               CloseLogFile                          (void);
INT                       SetLogFile                            (FILE *file);
INT               WriteLogFile                          (const char *text);

/* text output to shell window without log file */
void              WriteString                           (const char *s);

/* event input */
INT               GetNextUGEvent                        (EVENT *theEvent, INT Eventmask);

/* mouse functions */
void              MousePosition                         (INT *ScreenCoord);
INT               MouseStillDown                        (void);

/* tool name handling */
INT                       SetToolName                           (INT tool, const char *name);
const char   *GetToolName                               (INT tool);

#endif
