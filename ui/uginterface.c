// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  uginterface.c                                                                                                 */
/*																			*/
/* Purpose:   ug interface data structure manager							*/
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
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <assert.h>

/* low module */
#include "misc.h"
#include "evm.h"
#include "ugenv.h"
#include "devices.h"
#include "gm.h"
#include "wpm.h"
#include "wop.h"
#include "uginterface.h"
#include "cmdint.h"
#include "debug.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define INTERRUPT_CHAR          '.'

#define MAXCMDLEN                       256

#define POINTER                                 0               /* arrow tool can zoom pictures		*/
#define DRAG                                    1               /* arrow tool can drag pictures		*/
#define ZOOM                                    2               /* arrow tool can zoom pictures		*/

#define SEL_NODE                                0               /* hand tool can select nodes		*/
#define SEL_VECTOR                              1               /* hand tool can select vectors		*/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {

  /* fields for enironment list variable */
  ENVVAR v;

  /* cmd key specific stuff */
  char CommandName[MAXCMDLEN];                  /* command associated with the name */

} CMDKEY;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static UGWINDOW *currUgWindow;
static PICTURE *currPicture;

static INT autoRefresh;                                 /* ON or OFF						*/

static INT ArrowToolState=POINTER;              /* POINTER, ZOOM or DRAG			*/
static INT MarkToolState=RED;                   /* NO_REFINEMENT, RED, BLUE, COPY	*/
static INT HandToolState=SEL_NODE;              /* SEL_NODE, SEL_VECTOR				*/

static INT theCmdKeyDirID;                              /* env ID for the /Cmd Key dir		*/
static INT theCmdKeyVarID;                              /* env ID for the /Cmd Key dir		*/

static OUTPUTDEVICE *DefaultDevice;     /* our default ouput device             */

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*D
        SetCurrentPicture - make the specified picture the current picture

        SYNOPSIS:
        INT SetCurrentPicture (PICTURE *thePicture)

        PARAMETERS:
   .   thePicture - the picture that will be the current picture (may be 'NULL')

        DESCRIPTION:
        Make the specified picture the current picture (the frame of the current
        picture!='NULL' is highlighted and many actions refer by default to
        the current picture).

        RETURN VALUE:
        INT
   .n   0: ok
   D*/
/****************************************************************************/

INT SetCurrentPicture (PICTURE *thePicture)
{
  if (thePicture!=currPicture)
  {
    if (currPicture!=NULL)
    {
      DrawPictureFrame(currPicture,WOP_NOT_ACTIVE);
      InvalidateUgWindow(PIC_UGW(currPicture));
    }
    if (thePicture!=NULL)
    {
      DrawPictureFrame(thePicture,WOP_ACTIVE);
      InvalidateUgWindow(PIC_UGW(thePicture));
    }
  }
  currPicture = thePicture;

  return (0);
}

/****************************************************************************/
/*D
        GetCurrentPicture - return a pointer to the current picture

        SYNOPSIS:
        PICTURE *GetCurrentPicture (void)

        PARAMETERS:
        --

        DESCRIPTION:
        The function just returns a pointer to the current picture.

        RETURN VALUE:
        PICTURE *
   .n   a pointer to the current picture (may be 'NULL')
   D*/
/****************************************************************************/

PICTURE *GetCurrentPicture (void)
{
  return (currPicture);
}

/****************************************************************************/
/*D
        SetCurrentUgWindow - make the specified window the current window

        SYNOPSIS:
        INT SetCurrentUgWindow (UGWINDOW *theUgWindow)

        PARAMETERS:
   .   theUgWindow - the window that will be the current window (may be 'NULL')

        DESCRIPTION:
        Make the specified window the current window (many actions refer by default to
        the current window).

        RETURN VALUE:
        INT
   .n   0: ok
   D*/
/****************************************************************************/

INT SetCurrentUgWindow (UGWINDOW *theUgWindow)
{
  UGWINDOW *win;

  /* check validity */
  win = GetFirstUgWindow();
  if (win!=theUgWindow)
  {
    for (; win!=NULL; win=GetNextUgWindow(win))
      if (win==theUgWindow)
        break;
    if (win==NULL)
      return (1);
  }
  currUgWindow=theUgWindow;
  return (0);
}

/****************************************************************************/
/*D
        GetCurrentUgWindow - return a pointer to the current window

        SYNOPSIS:
        UGWINDOW *GetCurrentUgWindow (void)

        PARAMETERS:
        --

        DESCRIPTION:
        The function just returns a pointer to the current window.

        RETURN VALUE:
        UGWINDOW *
   .n   a pointer to the current window (may be 'NULL')
   D*/
/****************************************************************************/

UGWINDOW *GetCurrentUgWindow (void)
{
  return (currUgWindow);
}

/****************************************************************************/
/*																			*/
/* Function:  DoCmdKey														*/
/*																			*/
/* Purpose:   read out command string for command key						*/
/*																			*/
/* Input:	  char c: command key to read out								*/
/*			  char String: returned command string							*/
/*																			*/
/* Output:	  0: ok                                                                                                                 */
/*			  1: error														*/
/*																			*/
/****************************************************************************/

static INT DoCmdKey (char c, char *String)
{
  CMDKEY *theCmdKey;
  char theCmdKeyName[2],*s;

  /* find cmd key env item */
  theCmdKeyName[0] = c;
  theCmdKeyName[1] = '\0';
  theCmdKey = (CMDKEY*) SearchEnv(theCmdKeyName,"/Cmd Keys",theCmdKeyVarID,theCmdKeyDirID);
  if (theCmdKey != NULL)
  {
    strcpy(String, (const char *)theCmdKey->CommandName);
    for (s=String; *s!='\0'; s++)
      if (*s=='?')
        *s = '@';
    return (1);
  }
  return (0);
}

/****************************************************************************/
/*D
        SetCmdKey - create a command key

        SYNOPSIS:
        INT SetCmdKey (char c, const char *String)

        PARAMETERS:
   .   c - use this key in conjunction with the command key...
   .   String - ...to execute the command(s) specified here

        DESCRIPTION:
        Associate a command key with a command sequence. An environment
        item is created (or changed) in the '/Cmd Keys' directory that
        holds the specified commands.

        RETURN VALUE:
        INT
   .n   1: String too long, could not change to the '/Cmd Keys' dir,
   .n   could not create an command key environment item
   .n   0: ok

        SEE ALSO:
        DelCmdKey, ListCmdKeys, DelAllCmdKeys
   D*/
/****************************************************************************/

INT SetCmdKey (char c, const char *String)
{
  CMDKEY *theCmdKey;
  char theCmdKeyName[2];

  if (strlen(String)>=MAXCMDLEN)
    return (1);

  /* find cmd key env item */
  theCmdKeyName[0] = c;
  theCmdKeyName[1] = '\0';
  theCmdKey = (CMDKEY *) SearchEnv(theCmdKeyName,"/Cmd Keys",theCmdKeyVarID,theCmdKeyDirID);
  if (theCmdKey == NULL)
  {
    /* create cmd key */
    if (ChangeEnvDir("/Cmd Keys")==NULL)
      return(1);
    if ((theCmdKey=(CMDKEY *)MakeEnvItem(theCmdKeyName,theCmdKeyVarID,sizeof(CMDKEY)))==NULL)
      return(1);
  }
  strcpy(theCmdKey->CommandName, (const char *)String);
  return (0);
}

/****************************************************************************/
/*D
        DelCmdKey - remove a command key environment item

        SYNOPSIS:
        INT DelCmdKey (char c)

        PARAMETERS:
   .   c - the command character to delete

        DESCRIPTION:
        DelCmdKey removes the command key environment item associated with 'c'.

        RETURN VALUE:
        INT
   .n   1: could not remove the command key environment item
   .n   0: ok

        SEE ALSO:
        SetCmdKey, ListCmdKeys, DelAllCmdKeys
   D*/
/****************************************************************************/

INT DelCmdKey (char c)
{
  CMDKEY *theCmdKey;
  char theCmdKeyName[2];

  /* find cmd key env item */
  theCmdKeyName[0] = c;
  theCmdKeyName[1] = '\0';

  theCmdKey = (CMDKEY *) SearchEnv(theCmdKeyName,"/Cmd Keys",theCmdKeyVarID,theCmdKeyDirID);
  if (theCmdKey != NULL)
    if (RemoveEnvItem((ENVITEM *)theCmdKey))
      return (1);

  return (0);
}

/****************************************************************************/
/*D
        ListCmdKeys - list command keys with command sequences stored to shell

        SYNOPSIS:
        INT ListCmdKeys ()

        PARAMETERS:
        __

        DESCRIPTION:
        List command keys with command sequences stored to shell. That looks for
        example like':'

   .vb
   key command
   q quit
   ! <recent command>
   .ve


        RETURN VALUE:
        INT
   .n   1: could not change dir to '/Cmd Keys'
   .n   0: ok

        SEE ALSO:
        DelCmdKey, SetCmdKey, DelAllCmdKeys
   D*/
/****************************************************************************/

INT ListCmdKeys ()
{
  CMDKEY *theCmdKey;
  ENVDIR *theDir;
  char buffer[MAXCMDLEN+4];

  /* loop through and print all cmd keys */
  if ((theDir=ChangeEnvDir("/Cmd Keys"))==NULL)
    return (1);

  if (theDir->down == NULL)
    return(0);

  UserWrite("key command\n");
  for (theCmdKey = (CMDKEY *) theDir->down; theCmdKey!=NULL; theCmdKey=(CMDKEY *) NEXT_ENVITEM(theCmdKey))
    if (ENVITEM_TYPE(theCmdKey) == theCmdKeyVarID)
    {
      sprintf(buffer," %c  %s\n",ENVITEM_NAME(theCmdKey)[0],theCmdKey->CommandName);
      UserWrite(buffer);
    }

  return (0);
}

/****************************************************************************/
/*D
        DelAllCmdKeys - remove all command keys

        SYNOPSIS:
        INT DelAllCmdKeys (void)

        PARAMETERS:
        --

        DESCRIPTION:
        Remove all command key environment items

        RETURN VALUE:
        INT
   .n   1: could not change dir to '/Cmd Keys' or 'RemoveEnvItem' failed
   .n   0: ok

        SEE ALSO:
        DelCmdKey, SetCmdKey, ListCmdKeys
   D*/
/****************************************************************************/

INT DelAllCmdKeys (void)
{
  CMDKEY *theCmdKey;
  ENVDIR *theDir;

  /* loop through and delete all cmd keys */
  if ((theDir=ChangeEnvDir("/Cmd Keys"))==NULL)
    return (1);

  if (theDir->down == NULL)
    return(0);

  for (theCmdKey = (CMDKEY *) theDir->down; theCmdKey!=NULL; theCmdKey=(CMDKEY *) NEXT_ENVITEM(theCmdKey))
    if (ENVITEM_TYPE(theCmdKey) == theCmdKeyVarID)
      if (RemoveEnvItem((ENVITEM *)theCmdKey))
        return (1);

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  PrintEvent													*/
/*																			*/
/* Purpose:   print event													*/
/*																			*/
/* Input:	  VIEW *theView, INT Type										*/
/*																			*/
/* Output:	  none															*/
/*																			*/
/****************************************************************************/

static void PrintEvent (EVENT theEvent)
{
  switch (EVENT_TYPE(theEvent))
  {
  case EVENT_ERROR :
    UserWrite("EVENT_ERROR\n");
    break;
  case NO_EVENT :
    break;
  case TERM_GOAWAY :
    UserWrite("TERM_GOAWAY\n");
    break;
  case TERM_CMDKEY :
    UserWrite("TERM_CMDKEY\n");
    break;
  case TERM_STRING :
    UserWrite("TERM_STRING\n");
    break;
  case DOC_GOAWAY :
    UserWrite("DOC_GOAWAY\n");
    break;
  case DOC_ACTIVATE :
    UserWrite("DOC_ACTIVATE\n");
    break;
  case DOC_DRAG :
    UserWrite("DOC_DRAG\n");
    break;
  case DOC_GROW :
    UserWrite("DOC_GROW\n");
    break;
  case DOC_CHANGETOOL :
    UserWrite("DOC_CHANGETOOL\n");
    break;
  case DOC_CONTENTCLICK :
    UserWrite("DOC_CONTENTCLICK\n");
    break;
  case DOC_UPDATE :
    UserWrite("DOC_UPDATE\n");
    break;
  default :
    UserWrite("UNKNOWN\n");
    break;
  }
}

/****************************************************************************/
/*																			*/
/* Function:  ProcessEvent													*/
/*																			*/
/* Purpose:   process next EXT event or event from GUI						*/
/*																			*/
/* Input:	  none															*/
/*																			*/
/* Output:	  char String: string got from TERM_STRING - event				*/
/*																			*/
/* Return:	  INT : see below												*/
/*																			*/
/****************************************************************************/
/****************************************************************************/
/*D
        ProcessEvent - the event handler of ug

        SYNOPSIS:
        static INT ProcessEvent (char *String, INT EventMask)

        PARAMETERS:
   .   String - to store string got from TERM_STRING - event
   .   EventMask - specifies which events will be returned by 'GetNextUGEvent'
                possible values are EVERY_EVENT, TERM_STRING and TERM_CMDKEY

        DESCRIPTION:
        ProcessEvent is the event handler of ug. The event types it can handle
        and which are returned bye the interface function 'GetNextUGEvent' are':'

   .vb
   #define EVENT_ERROR                  0

   #define NO_EVENT				2
   #define TERM_GOAWAY                  3
   #define TERM_CMDKEY                  4
   #define TERM_STRING                  5
   #define DOC_GOAWAY				6
   #define DOC_ACTIVATE			7
   #define DOC_DRAG				8
   #define DOC_GROW				9
   #define DOC_CHANGETOOL			10
   #define DOC_CONTENTCLICK		11
   #define DOC_UPDATE				12
   .ve

   .   EVENT_ERROR              - an error occured in the interface function 'GetNextUGEvent'
   .   NO_EVENT			- no ug event (but may bay something to do for the output device
   .   TERM_GOAWAY              - close the shell window (and force ug to quit)
   .   TERM_CMDKEY              - a command key was typed
   .   TERM_STRING              - a command sequence typed into the shell was completed by typing <cr>
   .   DOC_GOAWAY			- close a ug graphics window
   .   DOC_ACTIVATE		- the first mouse click into a ug window produces an activate event
                                                        (and makes the corresponding window the current window)
   .   DOC_DRAG			- drag a ug graphics window
   .   DOC_GROW			- resize a ug graphics window (if the size changes in both directions
                                                        the pictures of that window will be resized accordingly; if it was
                                                        only one direction the scale pictures will not be changed)
   .   DOC_CHANGETOOL		- change the current mouse tool on the toolbar
   .   DOC_CONTENTCLICK	- a mouse click into a ug graphics window was encountered (the
                                                        action executed depends on the current mouse tool --> 'toolbar')
   .   DOC_UPDATE			- update event for a ug graphics window was encountered (if
                                                        the refresh state is on pictures will be reploted otherwise they will
                                                        just be invalidated)

        Always the last command string encountered from the shell is stored in the special command key "!"

        RETURN VALUE:
        INT
   .n   PE_STRING: command string encounterd from the shell
   .n   PE_INTERRUPT: user interrupt encountered (command key ".")
   .n   PE_NOTHING1: actually not of further interest
   .n   PE_NOTHING2: actually not of further interest
   .n   PE_OTHER: any other event occured (which is not of further interest)
   .n   PE_ERROR: an error occured

        SEE ALSO:
        UGWINDOW, PICTURE
   D*/
/****************************************************************************/

#define PE_STRING               0
#define PE_OTHER                1
#define PE_NOTHING1     2               /* Interface event: yes */
#define PE_NOTHING2     3               /* Interface event: no	*/
#define PE_INTERRUPT    4
#define PE_ERROR                5

static INT ProcessEvent (char *String, INT EventMask)
{
  EVENT theEvent;
  COORD qw, qh, scaling;
  UGWINDOW *theUgW;
  PICTURE *thePic;
  INT WinID, MousePosition[2], UGW_LLL_old[2], UGW_LUR_old[2], Offset[2];
  WORK theWork;

  if (GetNextUGEvent(&theEvent,EventMask))
    return (PE_ERROR);

  /* print event
     PrintEvent(theEvent); */

  switch (EVENT_TYPE(theEvent))
  {
  case NO_EVENT :
    if (!(EventMask&PE_INTERRUPT))
      /* do current work (not if UserInterrupt is calling) */
      for (theUgW=GetFirstUgWindow(); theUgW!=NULL; theUgW=GetNextUgWindow(theUgW))
      {
        if (UGW_VALID(theUgW)==NO) if (UpdateUgWindow(theUgW,currPicture)) return (PE_OTHER);
        if (autoRefresh)
          for (thePic=GetFirstPicture(theUgW); thePic!=NULL; thePic=GetNextPicture(thePic))
            if (PIC_VALID(thePic)==NO && VO_STATUS(PIC_VO(thePic))==ACTIVE)
            {
              if (DrawUgPicture(thePic))
              {
                autoRefresh = FALSE;
                PrintErrorMessage('W',"ProcessEvent","autorefresh is switched OFF");
                return (PE_OTHER);
              }
              if (thePic==currPicture) DrawPictureFrame(thePic,WOP_ACTIVE);
              else DrawPictureFrame(thePic,WOP_NOT_ACTIVE);
            }
      }

    if (theEvent.NoEvent.InterfaceEvent) return (PE_NOTHING1);
    return (PE_NOTHING2);
    break;
  case TERM_GOAWAY :
    /* tell interpreter to execute quit command */
    theEvent.Type = TERM_STRING;
    strcpy(String,"quit");
    break;
  case TERM_CMDKEY :
    if (DoCmdKey(theEvent.TermCmdKey.CmdKey, String))
    {
      theEvent.Type = TERM_STRING;
      strcpy(theEvent.TermString.String,(const char *)String);
      UserWrite(String);
      UserWrite("\n");
    }
    else if (theEvent.TermCmdKey.CmdKey==INTERRUPT_CHAR)
      return (PE_INTERRUPT);
  case TERM_STRING :
    assert (strlen(theEvent.TermString.String) < INPUTBUFFERLEN);
    strcpy(String,(const char *)theEvent.TermString.String);
    SetCmdKey('!',String);
    break;
  case DOC_GOAWAY :
    WinID = theEvent.DocGoAway.win;
    theUgW = WinID2UgWindow(WinID);
    if (theUgW == NULL) return (PE_OTHER);
    if ((currPicture!=NULL) && (PIC_UGW(currPicture)==theUgW))
      SetCurrentPicture(NULL);
    for (thePic=GetFirstPicture(theUgW); thePic!=NULL; thePic=GetFirstPicture(theUgW))
      if (DisposePicture(thePic))
        return (PE_OTHER);
    if (DisposeUgWindow(theUgW))
      return (PE_OTHER);
    if (theUgW == currUgWindow)
      if (SetCurrentUgWindow(GetFirstUgWindow())) return (PE_OTHER);
    break;
  case DOC_ACTIVATE :
    WinID = theEvent.DocActivate.win;
    SetCurrentUgWindow(WinID2UgWindow(WinID));
    if (currUgWindow == NULL) return (PE_OTHER);
    break;
  case DOC_DRAG :
    /* Update window position */
    WinID = theEvent.DocDrag.win;
    theUgW = WinID2UgWindow(WinID);
    if (theUgW == NULL) return (PE_OTHER);
    V2_COPY(theEvent.DocDrag.Global_LL,UGW_GLL(theUgW))
    V2_COPY(theEvent.DocDrag.Global_UR,UGW_GUR(theUgW))
    break;
  case DOC_GROW :
    WinID = theEvent.DocGrow.win;
    theUgW = WinID2UgWindow(WinID);
    if (theUgW == NULL) return (PE_OTHER);
    V2_COPY(UGW_LLL(theUgW),UGW_LLL_old)
    V2_COPY(UGW_LUR(theUgW),UGW_LUR_old)
    V2_COPY(theEvent.DocGrow.Global_LL,UGW_GLL(theUgW))
    V2_COPY(theEvent.DocGrow.Global_UR,UGW_GUR(theUgW))
    V2_COPY(theEvent.DocGrow.Local_LL,UGW_LLL(theUgW))
    V2_COPY(theEvent.DocGrow.Local_UR,UGW_LUR(theUgW))

    thePic=GetFirstPicture(theUgW);
    if (thePic == NULL) break;
    if (V2_ISEQUAL(UGW_LLL_old,PIC_GLL(thePic)) && V2_ISEQUAL(UGW_LUR_old,PIC_GUR(thePic)))
    {
      /* set new pixel range of picture */
      V2_COPY(UGW_LLL(theUgW),PIC_GLL(thePic))
      V2_COPY(UGW_LUR(theUgW),PIC_GUR(thePic))

      /* resize plane in physical space */
      if (ResizeViewPlane(PIC_VO(thePic),UGW_LLL_old,UGW_LUR_old,UGW_LLL(theUgW),UGW_LUR(theUgW))) return (PE_OTHER);
    }
    else
    {
      /* set pixel range of windows */
      qw = (COORD)(UGW_LUR(theUgW)[0]-UGW_LLL(theUgW)[0])/(COORD)(UGW_LUR_old[0]-UGW_LLL_old[0]);
      qh = (COORD)(UGW_LUR(theUgW)[1]-UGW_LLL(theUgW)[1])/(COORD)(UGW_LUR_old[1]-UGW_LLL_old[1]);
      if (qw==1.0 || qh==1.0) scaling = 1.0;
      else scaling = MIN(qw,qh);
      for (thePic=GetFirstPicture(theUgW); thePic!=NULL; thePic=GetNextPicture(thePic))
      {
        V2_SUBTRACT(PIC_GLL(thePic),UGW_LLL_old,Offset)
        V2_SCALE(scaling,Offset)
        V2_ADD(UGW_LLL(theUgW),Offset,PIC_GLL(thePic))
        V2_SUBTRACT(PIC_GUR(thePic),UGW_LLL_old,Offset)
        V2_SCALE(scaling,Offset)
        V2_ADD(UGW_LLL(theUgW),Offset,PIC_GUR(thePic))
      }
    }
    break;
  case DOC_CHANGETOOL :
    /* change tool */
    WinID = theEvent.DocChangeTool.win;
    theUgW = WinID2UgWindow(WinID);
    if (UGW_CURRTOOL(theUgW)==theEvent.DocChangeTool.Tool)
      switch (theEvent.DocChangeTool.Tool)
      {
      case arrowTool :
        switch (ArrowToolState)
        {
        case POINTER :
          SetToolName(arrowTool,"drag tool");
          ArrowToolState = DRAG;
          break;
        case DRAG :
          SetToolName(arrowTool,"zoom tool");
          ArrowToolState = ZOOM;
          break;
        default :
          SetToolName(arrowTool,arrowToolName);
          ArrowToolState = POINTER;
          break;
        }
        break;

      case gnoedelTool :
        switch (MarkToolState)
        {
        case NO_REFINEMENT :
          SetToolName(gnoedelTool,"mark red tool");
          MarkToolState = RED;
          break;
        case RED :
          SetToolName(gnoedelTool,"mark blue tool");
          MarkToolState = BLUE;
          break;
        case BLUE :
          SetToolName(gnoedelTool,"mark copy tool");
          MarkToolState = COPY;
          break;
        case COPY :
          SetToolName(gnoedelTool,"unmark tool");
          MarkToolState = NO_REFINEMENT;
          break;
        }
        break;

      case handTool :
        switch (HandToolState)
        {
        case SEL_NODE :
          SetToolName(handTool,"select vec tool");
          HandToolState = SEL_VECTOR;
          break;
        case SEL_VECTOR :
          SetToolName(handTool,"select nd tool");
          HandToolState = SEL_NODE;
          break;
        }
        break;
      }
    UGW_CURRTOOL(theUgW) = theEvent.DocChangeTool.Tool;
    InvalidateUgWindow(theUgW);
    if (PIC_UGW(currPicture)==theUgW)
      UpdateUgWindow(theUgW,currPicture);
    else
      UpdateUgWindow(theUgW,NULL);
    break;
  case DOC_CONTENTCLICK :
    WinID = theEvent.DocDrag.win;
    theUgW = WinID2UgWindow(WinID);
    V2_COPY(theEvent.DocContentClick.MousePosition,MousePosition)
    thePic = Mouse2Picture(theUgW,MousePosition);
    if (thePic == NULL) break;
    if (currPicture != thePic)
    {
      SetCurrentPicture(thePic);
      break;
    }

    switch (UGW_CURRTOOL(theUgW))
    {
    case arrowTool :
      switch (ArrowToolState)
      {
      case ZOOM :
        ZoomPicture(currPicture,MousePosition);
        break;
      case DRAG :
        DragPicture(currPicture,MousePosition);
        break;
      }
      break;
    case crossTool :
      W_ID(&theWork) = INSERTBNDNODE_WORK;
      W_INSERTBNDNODE_WORK(&theWork)->PixelX = MousePosition[0];
      W_INSERTBNDNODE_WORK(&theWork)->PixelY = MousePosition[1];
      if (WorkOnPicture(currPicture,&theWork)) return (PE_OTHER);
      break;
    case choiceTool :
      W_ID(&theWork) = MOVENODE_WORK;
      W_MOVENODE_WORK(&theWork)->PixelX = MousePosition[0];
      W_MOVENODE_WORK(&theWork)->PixelY = MousePosition[1];
      if (WorkOnPicture(currPicture,&theWork)) return (PE_OTHER);
      break;
    case circleTool :
      W_ID(&theWork) = INSERTNODE_WORK;
      W_INSERTNODE_WORK(&theWork)->PixelX = MousePosition[0];
      W_INSERTNODE_WORK(&theWork)->PixelY = MousePosition[1];
      if (WorkOnPicture(currPicture,&theWork)) return (PE_OTHER);
      break;
    case handTool :
      switch (HandToolState)
      {
      case SEL_NODE :
        W_ID(&theWork) = SELECTNODE_WORK;
        W_SELECTNODE_WORK(&theWork)->PixelX = MousePosition[0];
        W_SELECTNODE_WORK(&theWork)->PixelY = MousePosition[1];
        break;
      case SEL_VECTOR :
        W_ID(&theWork) = SELECTVECTOR_WORK;
        W_SELECTVECTOR_WORK(&theWork)->PixelX = MousePosition[0];
        W_SELECTVECTOR_WORK(&theWork)->PixelY = MousePosition[1];
        break;
      default :
        return (PE_OTHER);
      }
      if (WorkOnPicture(currPicture,&theWork)) return (PE_OTHER);
      break;
    case heartTool :
      W_ID(&theWork) = SELECTELEMENT_WORK;
      W_SELECTELEMENT_WORK(&theWork)->PixelX = MousePosition[0];
      W_SELECTELEMENT_WORK(&theWork)->PixelY = MousePosition[1];
      if (WorkOnPicture(currPicture,&theWork)) return (PE_OTHER);
      break;
    case gnoedelTool :
      W_ID(&theWork) = MARKELEMENT_WORK;
      W_MARKELEMENT_WORK(&theWork)->PixelX = MousePosition[0];
      W_MARKELEMENT_WORK(&theWork)->PixelY = MousePosition[1];
      W_MARKELEMENT_WORK(&theWork)->rule   = MarkToolState;
      if (WorkOnPicture(currPicture,&theWork)) return (PE_OTHER);
      break;
    default :
      break;
    }
    break;
  case DOC_UPDATE :
    WinID = theEvent.DocDrag.win;
    theUgW = WinID2UgWindow(WinID);
    if (InvalidatePicturesOfUgWindow(theUgW)) return(PE_OTHER);
    InvalidateUgWindow(theUgW);
    for (thePic=GetFirstPicture(theUgW); thePic!=NULL; thePic=GetNextPicture(thePic))
      if (thePic==currPicture) DrawPictureFrame(thePic,WOP_ACTIVE);
      else DrawPictureFrame(thePic,WOP_NOT_ACTIVE);

    break;
  }

  /* return */
  if (EVENT_TYPE(theEvent) == TERM_STRING)
    return (PE_STRING);

  return (PE_OTHER);
}

/****************************************************************************/
/*D
        UserInterrupt - check whether a user interrupt event was encounterd

        SYNOPSIS:
        INT UserInterrupt (const char *text)

        PARAMETERS:
   .   text - if an interrupt event was found and if 'text==NULL' 'YES' will be returned
                        otherwise a promt "### user-interrupt in <text>? will be prompted and 'YES'
                        will be returned only if a 'y' was entered into the shell

        DESCRIPTION:
        Check whether a user interrupt event was encounterd and return 'YES' or 'NO' correspondingly.
        If yes the mutelevel is set to 0 if it was < 0.

        RETURN VALUE:
        INT
   .n   YES: a user interrupt was encountered
   .n   NO:  no interrupt
   D*/
/****************************************************************************/

INT UserInterrupt (const char *text)
{
  INT Code,EventMask,mutelevel;
  char buffer[128];
        #ifdef ModelP
  int status;
  int fanout=1;
        #endif

        #ifndef STDIF
    #ifdef ModelP
  if (me == master)
  {
    #endif

  EventMask = TERM_CMDKEY;

  Code = ProcessEvent(buffer,EventMask);

        #ifdef ModelP
  /* if UserInterrupt called in InterpretString() then the interrupt  */
  /* is related to the master only !!									*/
  /* TODO: check this condition in newer releases						*/
  if (strcmp(text,"InterpretString")==0) fanout=0;

  if (fanout) Broadcast(&Code,sizeof(INT));
        #endif

  if (Code==PE_INTERRUPT)
  {
    if (text==NULL)
    {
                        #ifdef ModelP
      status = YES;
      if (fanout) Broadcast(&status,sizeof(int));
                        #endif
      return (YES);
    }
    else
    {
      mutelevel = GetMuteLevel();
      if (GetMuteLevel()<0)
        SetMuteLevel(0);
      sprintf(buffer,"### user-interrupt in '%s'?",text);
      UserWrite(buffer);
      UserRead(buffer);
      if (buffer[0]=='y')
      {
                #ifdef ModelP
        status = YES;
        if (fanout) Broadcast(&status,sizeof(int));
                                #endif
        return (YES);
      }
      else
      {
                #ifdef ModelP
        status = NO;
        if (fanout) Broadcast(&status,sizeof(int));
                                #endif
        SetMuteLevel(mutelevel);
        return (NO);
      }
    }
  }
    #ifdef ModelP
}
else
{
  Broadcast(&Code,sizeof(INT));
  if (Code==PE_INTERRUPT)
  {
    Broadcast(&status,sizeof(int));
    return (status);
  }
}
    #endif
        #endif /* STDIF */

  return (NO);
}

/****************************************************************************/
/*                                                                          */
/* Function:  ParExecCommand                                                */
/*                                                                          */
/* Purpose:   Broadcast a command line string to all processors, execute    */
/*            command on each processor and collect global status after     */
/*            termination.                                                  */
/*                                                                          */
/* Input:     pointer to the command line string                            */
/*                                                                          */
/* Output:    maximum of all return values on the different processors      */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP

int ParExecCommand (char *s)
{
  int error;
  int l,n;

  PRINTDEBUG(ui,4,("%d: ParExecCommand(%.30s)...\n",me,s))

  /* broadcast command line to all processors */
  PRINTDEBUG(ui,4,("%d:         Broadcast(%.30s)...\n",me,s))
  s[MAXCMDSIZE-1] = (char) 0;
  Broadcast(s,MAXCMDSIZE);

  /* execute command on each processor */
  PRINTDEBUG(ui,4,("%d:         ExecCommand(%.30s)...\n",me,s))
  error = ExecCommand(s);

  /* collect result code */
  PRINTDEBUG(ui,4,("%d:         (Get)Concentrate(%.30s)...\n",me,s))
  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(int));
    error = MAX(error,n);
  }
  Concentrate(&error,sizeof(int));

  /* fanout error code */
  PRINTDEBUG(ui,4,("%d:         Broadcast(%d)...\n",me,error))
  Broadcast(&error,sizeof(int));

  PRINTDEBUG(ui,4,("%d: ...end ParExecCommand(%.30s)...\n",me,s))

  /* return global status */
  return(error);
}

#endif

/****************************************************************************/
/*D
        UserIn - call process event until a string is enterd into the shell and return
                                the string (`all` events will be handled in the meantime)

        SYNOPSIS:
        INT UserIn (char *String)

        PARAMETERS:
   .   String - store shell string to be entered here

        DESCRIPTION:
        Call process event until a string is enterd into the shell and return
        the string (`all` events will be handled in the meantime, i.e it is possible
        for example to resize a graphics window etc.). It is called by the --> 'CommandLoop'.

        RETURN VALUE:
        INT
   .n   1: a process event error occured
   .n   0: ok
   D*/
/****************************************************************************/

INT UserIn (char *String)
{
  INT Code,EventMask;

  EventMask = EVERY_EVENT;

  /* loop till string is entered */
  while (TRUE)
  {
    Code = ProcessEvent(String,EventMask);
    if (Code == PE_ERROR) return (1);
    if (Code == PE_STRING)
    {
      WriteLogFile(String);
      return (0);
    }
  }
}

/****************************************************************************/
/*D
        UserRead - call process event until a string is enterd into the shell and return
                                the string (only 'TERM_STRING' events will be handled in the meantime)

        SYNOPSIS:
        INT UserRead (char *String)

        PARAMETERS:
   .   String - store shell string to be entered here

        DESCRIPTION:
        Call process event until a string is enterd into the shell and return
        the string (only 'TERM_STRING' events will be handled in the meantime, i.e it is possible
        for example to resize a graphics window etc.). It is called by --> 'InterpretString'
        and by 'UserInterrupt'.

        RETURN VALUE:
        INT
   .n   1: a process event error occured
   .n   0: ok
   D*/
/****************************************************************************/

INT UserRead (char *String)
{
  INT Code,EventMask;

  EventMask = TERM_STRING;

  /* loop till string is entered */
  while (TRUE)
  {
    Code = ProcessEvent(String,EventMask);
    assert (Code!=PE_ERROR);
    if (Code == PE_ERROR) return (1);
    if (Code == PE_STRING)
    {
      WriteLogFile(String);
      return (0);
    }
  }
}

/****************************************************************************/
/*D
        SetRefreshState -  determines wether pictures in 'UGWINDOW's will be updated

        SYNOPSIS:
        INT SetRefreshState (INT status)

        PARAMETERS:
   .   status - 'TRUE' or 'FALSE'

        DESCRIPTION:
        If autoRefresh is on then invalid pictures will be updated automatically
        during the command loop. In the other case nothing is done.

        We recommend to set refresh on only for displaying small grids or on very fast
        machines (especially be careful in 3D!) since every update event encountered
        for the window will force ug to replot all pictures contained.

        RETURN VALUE:
        INT
   .n   0: ok
   D*/
/****************************************************************************/

INT SetRefreshState (INT status)
{
  autoRefresh = status;
  return (0);
}

/****************************************************************************/
/*D
        InitUgInterface - initialize 'uginterface.c'

        SYNOPSIS:
        INT InitUgInterface ()

        PARAMETERS:
        --

        DESCRIPTION:
        The command key environment directory '/Cmd Keys' is created and the default
        output device is set locally.

        RETURN VALUE:
        INT
   .n   '__LINE__': failed to create the '/Cmd Keys' dir or the DefaultDevice is 'NULL'
   .n   0: ok

   D*/
/****************************************************************************/

INT InitUgInterface ()
{
  /* install the /Cmd Keys directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitUgInterface","could not changedir to root");
    return(__LINE__);
  }
  theCmdKeyDirID = GetNewEnvDirID();
  if (MakeEnvItem("Cmd Keys",theCmdKeyDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitUgInterface","could not install '/Cmd Keys' dir");
    return(__LINE__);
  }
  theCmdKeyVarID = GetNewEnvVarID();

  DefaultDevice = GetDefaultOutputDevice();

  return (0);
}
