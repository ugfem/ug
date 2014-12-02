// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugdevices.h                                                   */
/*                                                                          */
/* Purpose:   implements a simple but portable graphical user interface     */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   14.06.93 begin, ug version ug21Xmas3d                         */
/*            16.12.94 restructured for ug version 3.0                      */
/*                                                                          */
/* Remarks:   was "devices.h" in earlier version of UG                      */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DEVICESH__
#define __DEVICESH__

#include <cstdio>

#include "ugtypes.h"

#include "ugenv.h"

#include "namespace.h"

START_UG_NAMESPACE

/****************************************************************************/
/*                                                                          */
/*      constants                                                           */
/*                                                                          */
/****************************************************************************/

/* markers */
enum {EMPTY_SQUARE_MARKER,
      GRAY_SQUARE_MARKER,
      FILLED_SQUARE_MARKER,
      EMPTY_CIRCLE_MARKER,
      GRAY_CIRCLE_MARKER,
      FILLED_CIRCLE_MARKER,
      EMPTY_RHOMBUS_MARKER,
      GRAY_RHOMBUS_MARKER,
      FILLED_RHOMBUS_MARKER,
      PLUS_MARKER,
      CROSS_MARKER};

#define NMARKERS                                        11

/* tool numbers */
enum {arrowTool,
      crossTool,
      choiceTool,
      circleTool,
      handTool,
      heartTool,
      gnoedelTool};

#define nboftools 7

/* toolbox text len */
#define INFO_SIZE                               128
#define INFO_LEN                                127

#define NO_TOOL_CHOSEN                  -1              /* possible return value of WhichTool */

/* text position */
enum {TEXT_NOT_CENTERED,
      TEXT_CENTERED};

/* text modes */
enum {TEXT_REGULAR,
      TEXT_INVERSE,
      TEXT_INDEXED};

/* buffer space for command line input */
#define INPUTBUFFERLEN                  4096

/* event types */
enum EventType {EVENT_ERROR,
                EVERY_EVENT,
                NO_EVENT,
                TERM_GOAWAY,
                TERM_CMDKEY,
                TERM_STRING,
                DOC_GOAWAY,
                DOC_ACTIVATE,
                DOC_DRAG,
                DOC_GROW,
                DOC_CHANGETOOL,
                DOC_CONTENTCLICK,
                DOC_UPDATE};

enum UG_PALETTE {

  COLOR_PALETTE,
  BLACK_WHITE_PALETTE,
  GRAY_PALETTE
};

/****************************************************************************/
/*                                                                                                                                                      */
/*      macros                                                                                                                                  */
/*                                                                                                                                                      */
/****************************************************************************/

#define EVENT_TYPE(p)                   ((p).Type)

/****************************************************************************/
/*                                                                                                                                                      */
/*      data types                                                                                                                              */
/*                                                                                                                                                      */
/****************************************************************************/

/** \brief Identification of windows, used to carry a pointer (xif, ppm, ps)
   or INT (rif) */
typedef void* WINDOWID;

/** \brief Type for device coordinates */
typedef struct
{
  short x;
  short y;
} SHORT_POINT ;

/* function types exported by OUTPUTDEVICE */
typedef WINDOWID (*OpenOutputPtr)(const char *title, INT rename, INT x, INT y, INT width, INT height, INT *Global_LL, INT *Global_UR, INT *Local_LL, INT *Local_UR, INT *error);
typedef INT (*CloseOutputPtr)(WINDOWID win);
typedef INT (*ActivateOutputPtr)(WINDOWID win);
typedef INT (*UpdateOutputPtr)(WINDOWID win, INT tool);

/** \brief Data structure to define an interface to output devices

    The struct 'OUTPUTDEVICE' defines an interface to an output device with
    graphics capabilites. ug uses a default output device which usually is your
    computer monitor. Additionally there can be defined several other output devices.

    Up to now there is implemented an interface to XWindows of UNIX and to the
    Macintosh OS with the possibilties of window handling and plotting.
    They serve as default output device.

    Another output is the 'meta' ouput device. This is a format to write graphics commands
    to file which later can be view with the 'xugv' tool or that can be translated to PostScript
    format using the 'm2ps' tool. It is a lean storage format
    suited quite well for producing and viewing "films" with many pictures of time dependent solutions.
    It is also a helpful tool for production runs on large problems which often will run
    in batch mode.

    In the near future probably also a PostScript output device will exist.

    The output device struct requires functions for opening, closing, activating and updating
    a window on the device. Then there is a collection of graphics functions to set color,
    line width move the cursor, draw lines and higher level functions to plot filled polygons
    and text.

    Additionally there is information specified on the color palette used and some standard
    colors are defined.
 */
struct outputdevice {

  /** \brief This is an environment variable */
  ENVVAR v;

  long black;                                                   /* value for black                                                                              */
  long gray;                                /* value for gray                                       */
  long white;                                                   /* value for white                                                                              */
  long red;                                                             /* value for red                                                                                */
  long green;                                                   /* value for green                                                                              */
  long blue;                                                            /* value for blue                                                                               */
  long cyan;                                                            /* value for cyan                                                                               */
  long orange;                                                  /* value for orange                                                                     */
  long yellow;                                                  /* value for yellow                                                                     */
  long darkyellow;                                              /* value for yellow                                                                     */
  long magenta;                                                 /* value for magenta                                                                    */
  short hasPalette;                                             /* 1 if device has a color lookup ta                                    */
  long range;                                                   /* # of possible color indices                                                  */
  long spectrumStart;                                   /* usable range for a continuous                                                */
  long spectrumEnd;                                             /* color spectrum                                                                               */
  DOUBLE PixelRatio;                                            /* ratio of (physical) hight to width of a pixel                */
  short signx;                                                  /* direction of increasing x-coordinates                                */
  short signy;                                                  /* direction of increasing y-coordinates                                */

  /* pointers to basic drawing functions */
  void (*Move)(SHORT_POINT);                                                            /* move in device coordinates           */
  void (*Draw)(SHORT_POINT);                                                            /* draw from current point to given */
  void (*Polyline)(SHORT_POINT *, INT );                                        /* draw a polyline                                      */
  void (*InversePolyline)(SHORT_POINT *, INT );                         /* draw an inverted polyline            */
  void (*Polygon)(SHORT_POINT *, INT );                                         /* fill a polygon w. curr. col          */
  void (*ShadedPolygon)(SHORT_POINT *, INT, DOUBLE );           /* shade a polygon w. curr. col     */
  void (*InversePolygon)(SHORT_POINT *, INT );                          /* invert a polygon w. curr. col        */
  void (*ErasePolygon)(SHORT_POINT *, INT );                            /* erase a polygon w. curr. col         */
  void (*Polymark)(short, SHORT_POINT *);                                       /* place markers                                        */
  void (*InvPolymark)(short, SHORT_POINT *);                                    /* place inverse markers                        */
  void (*DrawText)(const char *, INT);                                          /* draw text in current size/font       */
  void (*CenteredText)(SHORT_POINT, const char *, INT);         /* draw text centered at x,y            */
  void (*ClearViewPort)(void);                                                          /* clear a view port                            */

  /* pointers to set functions */
  void (*SetLineWidth)(short);                                                          /* line width in pixels (points)        */
  void (*SetTextSize)(short);                                                           /* text size in pixels (points)         */
  void (*SetMarker)(short);                                                                     /* set marker id                                        */
  void (*SetMarkerSize)(short);                                                         /* marker size in pixels (points)       */
  void (*SetColor)(long);                                                                       /* arg is index or direct col value */
  void (*SetPaletteEntry)(long,short,short,short);                      /* set index to value                           */
  void (*SetNewPalette)(long,long,short*,short*,short*);        /* replace entrie                                       */

  /* pointers to miscellaneous functions */
  void (*GetPaletteEntry)(long,short *,short *,short *);        /* read color tab                                       */
  void (*Flush)(void);                                                                          /* flush graphics buffer                        */
  void (*PlotPixelBuffer)(void *, void *, int, int, int, int);
  /* plot a pixel buffer              */

  /* operations for managing windows */
  OpenOutputPtr OpenOutput;                       /* function to open a window                                                          */
  CloseOutputPtr CloseOutput;             /* function to close a window                                                         */
  ActivateOutputPtr ActivateOutput;       /* function to activate window                                                        */
  UpdateOutputPtr UpdateOutput;           /* function to draw outline with tool and info box            */
};



/** \brief No event */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  INT InterfaceEvent;                           /* 1 if the interface event was handled         */
  WINDOWID GraphWinActive;                      /* WINDOWID if a uggraphwin is active,0 else*/
  INT Mouse[2];                                         /* current mouse coord (rel. to window)         */
} NO_UGEVENT;

/** \brief Go away event for terminal window            */
typedef struct {
  INT Type;                                                     /* event type                                                           */
} TERM_GOAWAY_EVENT;

/** \brief Cmd key event for terminal window            */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  char CmdKey;                                          /* character from keyboard                                      */
} TERM_CMDKEY_EVENT;

/** \brief String event for terminal window             */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  char String[INPUTBUFFERLEN];          /* string from keyboard                                         */
} TERM_STRING_EVENT;

/** \brief Go away event for view                                       */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  WINDOWID win;                                         /* the window                                                           */
} DOC_GOAWAY_EVENT;

/** \brief Activate event for view                                      */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  WINDOWID win;                                         /* the window                                                           */
} DOC_ACTIVATE_EVENT;

/** \brief Drag event for view                                          */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  WINDOWID win;                                         /* the window                                                           */
  INT Global_LL[2];                                     /* new absolute position of window on screen*/
  INT Global_UR[2];                                     /*                                                                                      */
} DOC_DRAG_EVENT;

/** \brief Grow event for view                                          */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  WINDOWID win;                                         /* the window                                                           */
  INT Global_LL[2];                                     /* new absolute position of window on screen*/
  INT Global_UR[2];                                     /*                                                                                      */
  INT Local_LL[2];                                      /* range of pixels used for plotting            */
  INT Local_UR[2];                                      /*                                                                                      */
} DOC_GROW_EVENT;

/** \brief Change tool event for view                           */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  WINDOWID win;                                         /* the window                                                           */
  INT Tool;                                                     /* change to that tool                                          */
  INT MousePosition[2];                         /* mouse position                                                       */
} DOC_CHANGETOOL_EVENT;

/** \brief Content click event for view                         */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  WINDOWID win;                                         /* the window                                                           */
  INT MousePosition[2];                         /* mouse position                                                       */
} DOC_CONTENTCLICK_EVENT;

/** \brief Update event for view                                        */
typedef struct {
  INT Type;                                                     /* event type                                                           */

  /* data */
  WINDOWID win;                                         /* the window                                                           */
} DOC_UPDATE_EVENT;

/** \brief Data structure for the ug event handling

    The event handling of ug defines several possible event types the interface
    function 'GetNextUGEvent' can return to 'ProcessEvent'. Depending on the
    event type data are transferred by the corresponding component in the union.
    The events are distinguished by the 'Type' component in the 'EVENT'
 */
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
/*                                                                                                                                                      */
/* function exported by the devices module                                                                      */
/*                                                                                                                                                      */
/****************************************************************************/

/* initialization and clean up */
INT               InitDevices                           (int *argcp, char **argv);
INT           ExitDevices               (void);

/* return the size of the monitor screen, return 0 if not available */
INT               GetScreenSize                         (INT size[2]);

/* set/get mute level for output control */
void              SetMuteLevel                          (INT mute);
INT               GetMuteLevel                          (void);

/* multiple output devices */
OUTPUTDEVICE *CreateOutputDevice                (const char *name);
OUTPUTDEVICE *GetOutputDevice                   (const char *name);
OUTPUTDEVICE *GetDefaultOutputDevice    (void);

/* changing the palette */
INT                       UgSetPalette                          (OUTPUTDEVICE *dev, INT palette);

/* text output to shell with log file mechanism */
void              UserWrite                             (const char *s);
int               UserWriteF                            (const char *format, ...);
void              PrintErrorMessage             (char type, const char *procName, const char *text);
void              PrintErrorMessageF            (char type, const char *procName, const char *format, ...);
INT               OpenLogFile                           (const char *name, int rename);
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
void              DrawInfoBox                           (WINDOWID win, const char *info);
INT                       WhichTool                                     (WINDOWID win, const INT mouse[2], INT *tool);

END_UG_NAMESPACE

#endif
