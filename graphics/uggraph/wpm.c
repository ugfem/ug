// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  wpm.c                                                                                                                 */
/*																			*/
/* Purpose:   window picture manager										*/
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

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compiler.h"
#include "misc.h"
#include "evm.h"
#include "gm.h"
#include "num.h"
#include "wpm.h"
#include "general.h"

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
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static DOUBLE ex[3] = {1.0, 0.0, 0.0};
static DOUBLE ey[3] = {0.0, 1.0, 0.0};
static DOUBLE ez[3] = {0.0, 0.0, 1.0};

static INT theUgWindowsDirID;
static INT theUgWinDirID;

static INT thePicVarID;

static INT thePlotObjTypesDirID;
static INT thePlotObjTypesVarID;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*D
   CreatePicture - Allocate a new PICTURE

   SYNOPSIS:
   PICTURE *CreatePicture (const char *PictureName, UGWINDOW *theUgWindow,
   const INT *Global_LL, const INT *Global_UR);

   PARAMETERS:
   .  PictureName - name of the 'PICTURE' to be created
   .  theUgWindow - pointer to the 'UGWINDOW' in which the picture will be created
   .  Global_LL - LowerLeft corner of the 'PICTURE' in the pixelspace of the interior of
   the 'UGWINDOW'
   .  Global_UR - UpperRight corner of the 'PICTURE' in the pixelspace of the interior of
   the 'UGWINDOW'

   DESCRIPTION:
   This function allocates a new PICTURE in the specified 'UGWINDOW' with size and position
   specified by `Global_LL` and `Global_UR`

   RETURN VALUE:
   PICTURE *
   .n   pointer to 'PICTURE'
   .n   NULL if not created
   D*/
/****************************************************************************/

PICTURE *CreatePicture (const char *PictureName, UGWINDOW *theUgWindow, const INT *Global_LL, const INT *Global_UR)
{
  PICTURE *thePicture;
  INT sign;

  /* check if window exists */
  if (theUgWindow == NULL) return(NULL);

  /* allocate Image envItem */
  if (ChangeEnvDir("/UgWindows") == NULL) return (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theUgWindow)) == NULL) return (NULL);
  if (strlen(PictureName)>=NAMESIZE || strlen(PictureName)<1) return (NULL);
  if ((thePicture = (PICTURE *) MakeEnvItem(PictureName,thePicVarID,sizeof(PICTURE))) == NULL)
  {
    UserWrite("error: cannot create picture\n");
    return (NULL);
  }
  theUgWindow->NbPictures++;

  /* init picture */
  ENVITEM_LOCKED(thePicture)                      = NO;
  PIC_UGW(thePicture)                             = theUgWindow;
  PIC_VALID(thePicture)                           = NO;
  VO_STATUS(PIC_VO(thePicture))           = NOT_INIT;
  PO_POT(PIC_PO(thePicture))                      = NULL;

  sign = PIC_SIGN_X(thePicture) = SIGNUM(UGW_LUR(theUgWindow)[0] - UGW_LLL(theUgWindow)[0]);
  if (sign==0) return (NULL);
  PIC_GLL(thePicture)[0] = UGW_LLL(theUgWindow)[0] + sign*Global_LL[0];
  PIC_GUR(thePicture)[0] = UGW_LLL(theUgWindow)[0] + sign*Global_UR[0];
  sign = PIC_SIGN_Y(thePicture) = SIGNUM(UGW_LUR(theUgWindow)[1] - UGW_LLL(theUgWindow)[1]);
  if (sign==0) return (NULL);
  PIC_GLL(thePicture)[1] = UGW_LLL(theUgWindow)[1] + sign*Global_LL[1];
  PIC_GUR(thePicture)[1] = UGW_LLL(theUgWindow)[1] + sign*Global_UR[1];

  return (thePicture);
}

/****************************************************************************/
/*D
   DisposePicture - Dispose Picture

   SYNOPSIS:
   INT DisposePicture (PICTURE *thePicture);

   PARAMETERS:
   .  thePicture - the 'PICTURE' to dispose

   DESCRIPTION:
   This function disposes the 'PICTURE' and removes it from the 'UGWINDOW'-List.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/****************************************************************************/

INT DisposePicture (PICTURE *thePicture)
{
  UGWINDOW *theUgWindow;

  /* check if the image exists */
  if (thePicture == NULL) return (1);

  /* find ugwindow structure */
  if ((theUgWindow = PIC_UGW(thePicture)) == NULL) return (1);
  if (UGW_NPIC(theUgWindow) <= 0) return (1);

  /* dispose image */
  if (ChangeEnvDir("/UgWindows") == NULL) return (0);
  if (ChangeEnvDir(UGW_NAME(theUgWindow)) == NULL) return (0);
  if (RemoveEnvItem((ENVITEM*)thePicture)) return (1);
  UGW_NPIC(theUgWindow) -= 1;

  return (0);
}

void ResetToolBoxState (UGWINDOW *ugw)
{
  UGW_CURRTOOL(ugw)               = arrowTool;
  UGW_CURRFUNC(ugw)               = 0;
  UGW_INFOTEXT(ugw)[0]    = '\0';
  UGW_BOXSTATE(ugw)               = BOX_INVALID;
}

/****************************************************************************/
/*D
   CreateUgWindow - Allocate a new UGWINDOW

   SYNOPSIS:
   UGWINDOW *CreateUgWindow (OUTPUTDEVICE *theOutputDevice,
   const char *UgWindowName, INT x, INT y, INT width, INT height);

   PARAMETERS:
   .  theOutputDevice - 'OUTPUTDEVICE' on which the 'UGWINDOW' is opened (for example
   `screen` or `meta`)
   .  UgWindowName - the name of the 'UGWINDOW'
   .  x - pixel x-position of lowerleft-outter corner on the 'OUTPUTDEVICE'-screen
   .  y - pixel y-position of lowerleft-outter corner on the 'OUTPUTDEVICE'-screen
   .  width - pixel-width of the 'UGWINDOW' on the 'OUTPUTDEVICE'-screen
   .  height - pixel-height of the 'UGWINDOW' on the 'OUTPUTDEVICE'-screen

   DESCRIPTION:
   This function allocates a new UGWINDOW on the 'OUTPUTDEVICE'. The 'OUTPUTDEVICE'-
   screen is usually the monitore (so a quite finite pixelspace), for the meta-'OUTPUTDEVICE'
   it is the infinite pixel-space. The 'OUTPUTDEVICE'-screen of the monitore is defined to
   have its zero-point at its lower-left corner

   RETURN VALUE:
   UGWINDOW *
   .n   pointer to  'UGWINDOW'
   .n   NULL if cannot be created.
   D*/
/****************************************************************************/

UGWINDOW *CreateUgWindow (OUTPUTDEVICE *theOutputDevice, const char *UgWindowName, INT x, INT y, INT width, INT height)
{
  UGWINDOW *theWindow;
  WINDOWID winID;
  INT error;

  /* check outputdevice */
        #ifdef ModelP
  if (me == master)
        #endif
  if (theOutputDevice == NULL) return (NULL);

  /* allocate UgWindow envItem */
  if (ChangeEnvDir("/UgWindows") == NULL) return (NULL);
  if (strlen(UgWindowName)>=NAMESIZE || strlen(UgWindowName)<=1) return (NULL);
  if ((theWindow = (UGWINDOW *) MakeEnvItem(UgWindowName,theUgWinDirID,sizeof(UGWINDOW))) == NULL) return (NULL);

  /* open window on device and set sizes */
        #ifdef ModelP
  {
    INT size,l;
    INT *data;
    size = 8*sizeof(INT);
    data = (INT *)malloc(size);

    if (me == master)
    {
        #endif
  winID = (*theOutputDevice->OpenOutput)(UgWindowName, x, y, width, height, UGW_GLL(theWindow), UGW_GUR(theWindow), UGW_LLL(theWindow), UGW_LUR(theWindow), &error);
  if (error)
  {
    if (DisposeUgWindow(theWindow))
    {
      UserWrite("cannot open IFWindow: datastructure corrupted\n");
      return (NULL);
    }
    UserWrite("cannot open IFWindow\n");
    return (NULL);
  }
        #ifdef ModelP
  memcpy((void *)data,(void *)UGW_GLL(theWindow),size);

  /* spread coordinates to slaves */
  for(l=0; l<degree; l++)
    Spread(l,(void *)data,size);
}
else {
  GetSpread((void *)data,size);

  memcpy((void *)UGW_GLL(theWindow),(void *)data,size);

  /* get coordinates from master */
  for(l=0; l<degree; l++)
    Spread(l,(void *)data,size);
}
free(data);
}
        #endif

  /* set the other stuff */
  ENVITEM_LOCKED(theWindow)       = NO;
  UGW_NPIC(theWindow)             = 0;
  UGW_OUTPUTDEV(theWindow)        = theOutputDevice;
  UGW_VALID(theWindow)            = NO;
  UGW_IFWINDOW(theWindow)         = winID;

  return (theWindow);
}

/****************************************************************************/
/*D
    OpenPlacedPictures - Plot toolbox and if necessary infobox

   SYNOPSIS:
   INT UpdateUgWindow (UGWINDOW *theUgWindow, const PICTURE *EvalPicture);

   PARAMETERS:
   .  theUgWindow - update that 'UGWINDOW'
   .  EvalPicture - see DESCRIPTION

   DESCRIPTION:
   This function plots the toolbox of the 'UGWINDOW'. If the EvalPicture is a
   'PICTURE' of the theUgWindow and if it is initialized the info-box is
   plotted (the information comes from the corrisponding 'MULTIGRID').
   Has only effect on windows on the monitor ('screen').


   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

/* not used. see placer.c */
/*
   INT PlacePictures (PLACEMENT_TASK *task, PLACEMENT_REAL *real)
   {
    INT w,h,nx,ny,p,i,j,npic;

    w = task->winUR[0] - task->winLL[0] - 10;
    h = task->winUR[1] - task->winLL[1] - 10;;

    for (p=1000; p>90; p-=10)
    {
        nx = w/p;
        ny = h/p;
        if (nx*ny>=task->n) break;
        }
    if (p==90) return (1);

    npic=0;
    for (i=0; i<nx; i++)
       for (j=0; j<ny; j++)
       {
           real->picLL[npic][0] = 10+i*p;
           real->picLL[npic][1] = 10+j*p;
           real->picUR[npic][0] = (i+1)*p;
           real->picUR[npic][1] = (j+1)*p;
           npic++;
           if (npic>=task->n) return (0);
           }

        return (0);
   }
 */

INT OpenPlacedPictures (OUTPUTDEVICE *theOutputDevice, PLACEMENT_TASK *task)
{
  INT i,j;
  PLACEMENT_REAL real;
  UGWINDOW *theWin;
  PICTURE *thePic[WPM_PLM_PMAX];

  /* check */
  if (task->n<1) return (1);

  /* place pictures */
  if (PlacePictures(task,&real)) return (1);

  /* realize pictures */
  theWin = CreateUgWindow(theOutputDevice,task->win_name,real.winLL[0],real.winLL[1],real.winUR[0]-real.winLL[0],real.winUR[1]-real.winLL[1]);
  if (theWin==NULL) return (1);
  for (i=0; i<task->n; i++)
  {
    thePic[i] = CreatePicture (task->pic_name[i],theWin,real.picLL[i],real.picUR[i]);
    if (thePic[i]==NULL)
    {
      for (j=0; j<i; j++)
        DisposePicture(thePic[j]);
      return (1);
    }
  }

  return (0);
}

/****************************************************************************/
/*D
   UpdateUgWindow - Plot toolbox and if necessary infobox

   SYNOPSIS:
   INT UpdateUgWindow (UGWINDOW *theUgWindow, const PICTURE *EvalPicture);

   PARAMETERS:
   .  theUgWindow - update that 'UGWINDOW'
   .  EvalPicture - see DESCRIPTION

   DESCRIPTION:
   This function plots the toolbox of the 'UGWINDOW'. If the EvalPicture is a
   'PICTURE' of the theUgWindow and if it is initialized the info-box is
   plotted (the information comes from the corrisponding 'MULTIGRID').
   Has only effect on windows on the monitor ('screen').


   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT UpdateUgWindow (UGWINDOW *theUgWindow, const PICTURE *EvalPicture)
{
  if (theUgWindow==NULL) return (0);

  /* update ugwindow */

  if ((*UGW_OUTPUTDEV(theUgWindow)->UpdateOutput)(UGW_IFWINDOW(theUgWindow),UGW_CURRTOOL(theUgWindow)))
    return (1);

  /* window is valid */
  UGW_VALID(theUgWindow) = YES;

  return (0);
}

/****************************************************************************/
/*D
   DisposeUgWindow - Dispose a window

   SYNOPSIS:
   INT DisposeUgWindow (UGWINDOW *theUgWindow);

   PARAMETERS:
   .  theUgWindow - the 'UGWINDOW' to dispose

   DESCRIPTION:
   This function disposes a window and excludes it from window list of
   document, closes associated output. If there is still a 'PICTURE' open
   on the 'UGWINDOW' it returns an error-code.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT DisposeUgWindow (UGWINDOW *theUgWindow)
{
  OUTPUTDEVICE *OutputDevice;

  /* check if there are no pictures on the UgWindow */
  if (UGW_NPIC(theUgWindow) != 0) return (1);

  /* find output device */
  OutputDevice = UGW_OUTPUTDEV(theUgWindow);
  if (OutputDevice == NULL) return (1);

  /* close associated IFWindow */
  if ((*OutputDevice->CloseOutput)(theUgWindow->theIFWindow)) return (1);

  /* dispose window */
  if (ChangeEnvDir("/UgWindows") == NULL) return (1);
  if (RemoveEnvItem((ENVITEM*)theUgWindow)) return (1);

  return (0);
}

/****************************************************************************/
/*D
   GetUgPicture - Get picture of UgWindow by name

   SYNOPSIS:
   PICTURE *GetUgPicture (const UGWINDOW *theUgWindow, const char *name);

   PARAMETERS:
   .  theUgWindow - searches a 'PICTURE' on that 'UGWINDOW'
   .  name - searches 'PICTURE' with that name

   DESCRIPTION:
   This function gets picture of UgWindow by name.

   RETURN VALUE:
   PICTURE *
   .n     pointer to PICTURE
   .n     NULL if there is no with that name.
   D*/
/****************************************************************************/

PICTURE *GetUgPicture (const UGWINDOW *theUgWindow, const char *name)
{
  if (ChangeEnvDir("/UgWindows")==NULL) return(NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theUgWindow)) == NULL) return (NULL);
  return((PICTURE*) SearchEnv(name,".",thePicVarID,SEARCHALL));
}

/****************************************************************************/
/*D
   GetFirstPicture - Get the first picture of UgWindow

   SYNOPSIS:
   PICTURE *GetFirstPicture (const UGWINDOW *theUgWindow);

   PARAMETERS:
   .  theUgWindow - the 'UGWINDOW'

   DESCRIPTION:
   This function gets the first picture in the UgWindow-List.

   RETURN VALUE:
   PICTURE *
   .n      pointer to PICTURE
   .n      NULL if there is no
   D*/
/****************************************************************************/

PICTURE *GetFirstPicture (const UGWINDOW *theUgWindow)
{
  ENVITEM *thePicture;

  if (theUgWindow == NULL) return (NULL);
  for(thePicture=((ENVDIR*)theUgWindow)->down; thePicture!=NULL; thePicture = thePicture->v.next)
    if (thePicture->v.type == thePicVarID)
      return ((PICTURE*)thePicture);
  return (NULL);
}

/****************************************************************************/
/*D
   GetNextPicture - Get the next picture of UgWindow

   SYNOPSIS:
   PICTURE *GetNextPicture (const PICTURE *thePicture);

   PARAMETERS:
   .  thePicture -

   DESCRIPTION:
   This function gets the next 'PICTURE' in the picture-list of the UgWindow
   of 'thePicture'.

   RETURN VALUE:
   PICTURE *
   .n      pointer to PICTURE *
   .n      NULL if there is no
   D*/
/****************************************************************************/

PICTURE *GetNextPicture (const PICTURE *thePicture)
{
  ENVITEM *theNextPicture;

  if (thePicture == NULL) return (NULL);
  for(theNextPicture=NEXT_ENVITEM((ENVITEM*)thePicture); theNextPicture!=NULL; theNextPicture = NEXT_ENVITEM(theNextPicture))
    if (ENVITEM_TYPE(theNextPicture) == thePicVarID)
      return ((PICTURE*)theNextPicture);
  return (NULL);
}

/****************************************************************************/
/*D
   ListWindowPictureHeader - print information head

   SYNOPSIS:
   void ListWindowPictureHeader (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function prints (on ug shell) the information-head for UgWindow- and Picture-List.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

#define WPL_FORMAT                              "%-2.1s%-15.12s%-15.12s%-15.12s%-15.120s%-15.12s%-15.12s\n"

void ListWindowPictureHeader (void)
{
  UserWriteF(WPL_FORMAT,"","UgWindow","Picture","VO_Status","PlotObjType","PO_Status","Multigrid");
  UserWriteF(WPL_FORMAT,"","--------","-------","---------","-----------","---------","---------");
  return;
}
/****************************************************************************/
/*D
   ListUgWindow - List information about UgWindow

   SYNOPSIS:
   void ListUgWindow (const UGWINDOW *theUgWindow, INT current);

   PARAMETERS:
   .  theUgWindow - Information about that UgWindow
   .  current - 0 if the UgWindow is not current

   DESCRIPTION:
   This function lists the information about UgWindow and plots an hash-symbol if
   it is told to be current.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/
void ListUgWindow (const UGWINDOW *theUgWindow, INT current)
{
  if (current) UserWriteF(WPL_FORMAT,"#",ENVITEM_NAME(theUgWindow),"","","","","");
  else UserWriteF(WPL_FORMAT,"",ENVITEM_NAME(theUgWindow),"","","","","");
  return;
}
/****************************************************************************/
/*D
   ListPicture - List information and Picture

   SYNOPSIS:
   void ListPicture (const PICTURE *thePicture, INT current);

   PARAMETERS:
   .  thePicture - Information about that picture.
   .  current - 0 if the picture is not current

   DESCRIPTION:
   This function lists information the Picture: Its dimension, its 'STATUS'
   ('NOT_INIT', 'NOT_ACTIVE' or 'ACTIVE'), name of the 'FORMAT', name
   of the 'UGWINDOW' and its own name are plotted. If 'current' is 1 an
   asterisk is plotted.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/
void ListPicture (const PICTURE *thePicture, INT current)
{
  UGWINDOW *theUgW;
  char b1[2], b2[11], b3[30], b4[30], b5[30];
  INT VO_Status, PO_Status;
  int PO_Dim;

  theUgW = PIC_UGW(thePicture);
  VO_Status = VO_STATUS(PIC_VO(thePicture));
  PO_Status = PO_STATUS(PIC_PO(thePicture));
  if (current) sprintf(b1,"%s","*");
  else sprintf(b1,"%s","");
  switch (VO_Status)
  {
  case NOT_INIT :
    sprintf(b2,"%s","NOT_INIT");
    break;
  case NOT_ACTIVE :
    sprintf(b2,"%s","NOT_ACTIVE");
    break;
  case ACTIVE :
    sprintf(b2,"%s","ACTIVE");
    break;
  default :
    return;
  }
  switch (PO_DIM(PIC_PO(thePicture)))
  {
  case NOT_DEFINED :
    break;
  case TYPE_2D :
    PO_Dim = 2;
    break;
  case TYPE_3D :
    PO_Dim = 3;
    break;
  default :
    return;
  }
  switch (PO_Status)
  {
  case NOT_INIT :
    sprintf(b3,"---");
    sprintf(b4,"%s","NOT_INIT");
    sprintf(b5,"---");
    break;
  case NOT_ACTIVE :
    sprintf(b3,"%s",ENVITEM_NAME(PIC_POT(thePicture)));
    sprintf(b4,"%s:%dD","NOT_ACTIVE",PO_Dim);
    sprintf(b5,"%s",ENVITEM_NAME(PIC_MG(thePicture)));
    break;
  case ACTIVE :
    sprintf(b3,"%s",ENVITEM_NAME(PIC_POT(thePicture)));
    sprintf(b4,"%s:%dD","ACTIVE",PO_Dim);
    sprintf(b5,"%s",ENVITEM_NAME(PIC_MG(thePicture)));
    break;
  default :
    return;
  }
  UserWriteF(WPL_FORMAT,b1,ENVITEM_NAME(theUgW),ENVITEM_NAME(thePicture),b2,b3,b4,b5);
  return;
}

/****************************************************************************/
/*D
   WinID2UgWindow - Get UgWindow by its id

   SYNOPSIS:
   UGWINDOW *WinID2UgWindow (WINDOWID id);

   PARAMETERS:
   .  id - the ID, a reference the to 'DEVICES'-window-management

   DESCRIPTION:
   This function gets UgWindow by its ID.

   RETURN VALUE:
   UGWINDOW *
   .n      pointer to UGWINDOW
   .n      NULL if there is no with name
   D*/
/****************************************************************************/

UGWINDOW *WinID2UgWindow (WINDOWID id)
{
  ENVITEM *theUgWindow;
  ENVDIR *theUgWRoot;

  if ((theUgWRoot=ChangeEnvDir("/UgWindows"))==NULL) return(NULL);
  for (theUgWindow=ENVDIR_DOWN(theUgWRoot); theUgWindow!=NULL; theUgWindow=NEXT_ENVITEM(theUgWindow))
    if (ENVITEM_TYPE(theUgWindow) == theUgWinDirID)
      if (UGW_IFWINDOW((UGWINDOW*)theUgWindow) == id)
        return ((UGWINDOW*)theUgWindow);

  return (NULL);
}

/****************************************************************************/
/*D
   GetUgWindow - Get UgWindow by name

   SYNOPSIS:
   UGWINDOW *GetUgWindow (const char *name);

   PARAMETERS:
   .  name - find 'UGWINDOW' with that name

   DESCRIPTION:
   This function gets UgWindow by name.

   RETURN VALUE:
   UGWINDOW *
   .n     pointer to UGWINDOW *
   .n     NULL if there is no with name
   D*/
/****************************************************************************/

UGWINDOW *GetUgWindow (const char *name)
{
  if (ChangeEnvDir("/UgWindows")==NULL) return(NULL);
  return((UGWINDOW*) SearchEnv(name,".",theUgWinDirID,SEARCHALL));
}

/****************************************************************************/
/*D
   GetFirstUgWindow - Get first UgWindow

   SYNOPSIS:
   UGWINDOW *GetFirstUgWindow (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function gets the first UgWindow in the internal list of UgWindows, located
   in the 'ENVDIR' `UgWindows`

   RETURN VALUE:
   UGWINDOW *
   .n     pointer to UGWINDOW
   .n     NULL if there is no
   D*/
/****************************************************************************/

UGWINDOW *GetFirstUgWindow (void)
{
  ENVITEM *theUgWindow;

  if ((theUgWindow=(ENVITEM*)ChangeEnvDir("/UgWindows")) == NULL) return (NULL);

  for (theUgWindow=ENVITEM_DOWN(theUgWindow); theUgWindow!=NULL; theUgWindow=NEXT_ENVITEM(theUgWindow))
    if (ENVITEM_TYPE(theUgWindow) == theUgWinDirID)
      return ((UGWINDOW*)theUgWindow);
  return (NULL);
}

/****************************************************************************/
/*D
   GetNextUgWindow - Get next UgWindow

   SYNOPSIS:
   UGWINDOW *GetNextUgWindow (const UGWINDOW *theUgWindow);

   PARAMETERS:
   .  theUgWindow - get the next UgWindow

   DESCRIPTION:
   This function gets next UgWindow (of 'theUgWindow') in the internal list
   of UgWindows, located in the 'ENVDIR' `UgWindows`

   RETURN VALUE:
   UGWINDOW *
   .n     pointer to UGWINDOW
   .n     NULL if there is no
   D*/
/****************************************************************************/

UGWINDOW *GetNextUgWindow (const UGWINDOW *theUgWindow)
{
  ENVITEM *theNextUgWindow;

  for (theNextUgWindow=theUgWindow->d.next; theNextUgWindow!=NULL; theNextUgWindow=theNextUgWindow->d.next )
    if (theNextUgWindow->d.type == theUgWinDirID)
      return ((UGWINDOW*)theNextUgWindow);
  return (NULL);
}

/****************************************************************************/
/*D
   Mouse2Picture - Find picture in UgWindow

   SYNOPSIS:
   PICTURE *Mouse2Picture (const UGWINDOW *theUgWindow, INT *MousePosition);

   PARAMETERS:
   .  theUgWindow - find in that 'PICTURE'
   .  MousePosition - pixel coordinates of mouse

   DESCRIPTION:
   This function finds picture in UgWindow by mouse-position

   RETURN VALUE:
   PICTURE *
   .n     pointer to PICTURE
   .n     NULL if no picture hit by mouse
   D*/
/****************************************************************************/

PICTURE *Mouse2Picture (const UGWINDOW *theUgWindow, INT *MousePosition)
{
  PICTURE *thePicture;
  DOUBLE a;

  for (thePicture=GetFirstPicture(theUgWindow); thePicture!=NULL; thePicture=GetNextPicture(thePicture))
  {
    a = ((DOUBLE)(MousePosition[0]-PIC_GLL(thePicture)[0]))/((DOUBLE)(PIC_GUR(thePicture)[0]-PIC_GLL(thePicture)[0]));
    if (a>0.0 && a<1.0)
    {
      a = ((DOUBLE)(MousePosition[1]-PIC_GLL(thePicture)[1]))/((DOUBLE)(PIC_GUR(thePicture)[1]-PIC_GLL(thePicture)[1]));
      if (a>0.0 && a<1.0)
        return (thePicture);
    }
  }

  return (NULL);
}

/****************************************************************************/
/*D
   InvalidatePicture -  invalidate a picture

   SYNOPSIS:
   INT InvalidatePicture (PICTURE *thePicture);

   PARAMETERS:

   .  thePicture - the picture to invalidate

   DESCRIPTION:
   This function marks thePicture for update, i.e. tells that the content of the
   picture has changed.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT InvalidatePicture (PICTURE *thePicture)
{
  PIC_VALID(thePicture) = NO;

  return (0);
}

/****************************************************************************/
/*D
   InvalidatePicturesOfUgWindow	- invalidate pictures of UgWindow

   SYNOPSIS:
   INT InvalidatePicturesOfUgWindow (UGWINDOW *theUgW);

   PARAMETERS:
   .  theUgW - invalidate pictures of that 'UGWINDOW'

   DESCRIPTION:
   This function marks all pictures of ugwindow for update i.e. tells that the
   content of the pictures has changed.
   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT InvalidatePicturesOfUgWindow (UGWINDOW *theUgW)
{
  PICTURE *thePicture;

  for (thePicture=GetFirstPicture(theUgW); thePicture!=NULL; thePicture=GetNextPicture(thePicture))
    PIC_VALID(thePicture) = NO;

  return (0);
}

/****************************************************************************/
/*D
   InvalidatePicturesOfMG - Invalidate all pictures of a 'MULTIGRID'

   SYNOPSIS:
   INT InvalidatePicturesOfMG (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - pointer to multigrid

   DESCRIPTION:
   This function marks all pictures of that multigrid for update, i.e. tells
   that the content of the pictures has changed.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/************************************************************************/

INT InvalidatePicturesOfMG (MULTIGRID *theMG)
{
  UGWINDOW *theUgWindow;
  PICTURE *thePicture;

  for (theUgWindow=GetFirstUgWindow(); theUgWindow!=NULL; theUgWindow=GetNextUgWindow(theUgWindow))
    for (thePicture=GetFirstPicture(theUgWindow); thePicture!=NULL; thePicture=GetNextPicture(thePicture))
      if (PO_MG(PIC_PO(thePicture)) == theMG)
        PIC_VALID(thePicture) = NO;

  return (0);
}

/****************************************************************************/
/*D
   InvalidateUgWindow - Invalidate UgWindow

   SYNOPSIS:
   INT InvalidateUgWindow (UGWINDOW *theUgWindow);

   PARAMETERS:
   .  theUgWindow -

   DESCRIPTION:
   This function invalidates a UgWindow, i.e. tell that its 'TOOLBOX' or 'Infobox'
   has changed.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/******************************************************************************/

INT InvalidateUgWindow (UGWINDOW *theUgWindow)
{
  UGW_VALID(theUgWindow) = NO;
  return (0);
}

/****************************************************************************/
/*D
   InvalidateUgWindowsOfMG - invalidate UgWindow

   SYNOPSIS:
   INT InvalidateUgWindowsOfMG (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG -

   DESCRIPTION:
   This function invalidates all UgWindows of the multigrid, i.e. tells that the
   'TOOLBOX' or 'InfoBox'has changed.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT InvalidateUgWindowsOfMG (MULTIGRID *theMG)
{
  UGWINDOW *theUgW;
  PICTURE *thePic;
  INT found;

  if (theMG==NULL) return (0);

  for (theUgW=GetFirstUgWindow(); theUgW!=NULL; theUgW=GetNextUgWindow(theUgW))
  {
    found=0;
    for (thePic=GetFirstPicture(theUgW); thePic!=NULL; thePic=GetNextPicture(thePic))
      if (PIC_MG(thePic) == theMG)
      {
        found=1;
        break;
      }
    if (found)
      UGW_VALID(theUgW) = NO;
  }

  return (0);
}

/****************************************************************************/
/*D
   MovePictureToNewWindow - create a new window and move the picture to it

   SYNOPSIS:
   INT MovePictureToNewWindow (PICTURE *pic)

   PARAMETERS:
   .  pic - picture to be moved

   DESCRIPTION:
   This function creates a new window with appropriate size and moves the picture to it.

   RETURN VALUE:
   INT
   .n     0: ok
   .n     1: could not create new window
   D*/
/****************************************************************************/

INT MovePictureToNewWindow (PICTURE *pic)
{
  UGWINDOW *oldWin,*newWin;
  INT x,y,w,h;

  oldWin = PIC_UGW(pic);
  x = 10;
  y = 10;
  w = fabs(PIC_GUR(pic)[0] - PIC_GLL(pic)[0]);
  h = fabs(PIC_GUR(pic)[1] - PIC_GLL(pic)[1]);
  if ((newWin=CreateUgWindow(UGW_OUTPUTDEV(oldWin),PIC_NAME(pic),x,y,w,h))==NULL)
    return (1);

  /* move picture to new window */
  MoveEnvItem((ENVITEM*)pic,(ENVDIR*)oldWin,(ENVDIR*)newWin);
  PIC_UGW(pic) = newWin;
  UGW_NPIC(oldWin)--;
  UGW_NPIC(newWin)++;

  /* set new coordinates of the picture */
  V2_COPY(UGW_LLL(newWin),PIC_GLL(pic));
  V2_COPY(UGW_LUR(newWin),PIC_GUR(pic));

  /* remove old window if empty */
  if (UGW_NPIC(oldWin)==NULL)
    if (DisposeUgWindow(oldWin))
      return (2);

  return (0);
}

/****************************************************************************/
/*D
   GetPlotObjType - Get PLOTOBJTYPE from name

   SYNOPSIS:
   PLOTOBJTYPE *GetPlotObjType (const char *PlotObjTypeName);

   PARAMETERS:
   .  PlotObjTypeName -

   DESCRIPTION:
   This function gets the PLOTOBJTYPE from name. Used in the routines for
   initializing 'PLOTOBJ's

   RETURN VALUE:
   PLOTOBJTYPE *
   .n     pointer to PLOTOBJTYPE
   .n     NULL if not existing
   D*/
/****************************************************************************/

PLOTOBJTYPE *GetPlotObjType (const char *PlotObjTypeName)
{
  if (ChangeEnvDir("/PlotObjTypes")==NULL) return(NULL);
  return((PLOTOBJTYPE*) SearchEnv(PlotObjTypeName,".",thePlotObjTypesVarID,SEARCHALL));
}

/****************************************************************************/
/*D
   GetFirstPlotObjType - Get first plot object

   SYNOPSIS:
   PLOTOBJTYPE *GetFirstPlotObjType (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function gets first plot object in the internal list of the ug,
   located in the 'ENVDIR' `PlotObjTypes`

   RETURN VALUE:
   PLOTOBJTYPE *
   .n      pointer to PLOTOBJTYPE
   .n      NULL if there is no
   D*/
/****************************************************************************/

PLOTOBJTYPE *GetFirstPlotObjType (void)
{
  ENVITEM *thePlotObj;

  if ((thePlotObj=(ENVITEM*)ChangeEnvDir("/PlotObjTypes")) == NULL) return (NULL);

  for (thePlotObj=ENVITEM_DOWN(thePlotObj); thePlotObj!=NULL; thePlotObj=NEXT_ENVITEM(thePlotObj))
    if (ENVITEM_TYPE(thePlotObj) == thePlotObjTypesVarID)
      return ((PLOTOBJTYPE*)thePlotObj);
  return (NULL);
}

/****************************************************************************/
/*D
   GetNextPlotObjType - Get next plotobjecttype following thePlotObjType

   SYNOPSIS:
   PLOTOBJTYPE *GetNextPlotObjType (const PLOTOBJTYPE *thePlotObj);

   PARAMETERS:
   .  thePlotObjType -

   DESCRIPTION:
   This function gets next plot object following thePlotObjType in the internal
   list of 'PLOTOBJTYPES', located in the 'ENVDIR' `PlotObjTypes`

   RETURN VALUE:
   PLOTOBJTYPE
   .n       pointer to PLOTOBJTYPE
   .n       NULL if there is no
   D*/
/****************************************************************************/

PLOTOBJTYPE *GetNextPlotObjType (const PLOTOBJTYPE *thePlotObjType)
{
  ENVITEM *theNextPlotObj;

  for (theNextPlotObj=thePlotObjType->v.next; theNextPlotObj!=NULL; theNextPlotObj=theNextPlotObj->v.next )
    if (theNextPlotObj->v.type == thePlotObjTypesVarID)
      return ((PLOTOBJTYPE*)theNextPlotObj);
  return (NULL);
}

/****************************************************************************/
/*D
   CreatePlotObjType - Create PLOTOBJTYPE with name

   SYNOPSIS:
   PLOTOBJTYPE *CreatePlotObjType (const char *PlotObjTypeName, INT size);

   PARAMETERS:
   .  PlotObjTypeName - cretae PlotObjType with that name
   .  size - total size of the PlotObjType

   DESCRIPTION:
   This function creates a struct 'PLOTOBJTYPE' with specified name, its size is variable
   according to an extension to a 'PLOTOBJHANDLING'.

   RETURN VALUE:
   PLOTOBJTYPE *
   .n     pointer to PLOTOBJTYPE
   .n     NULL if not existing
   D*/
/****************************************************************************/

PLOTOBJTYPE *CreatePlotObjType (const char *PlotObjTypeName, INT size)
{
  PLOTOBJTYPE *pot;

  /* change to directory */
  if (ChangeEnvDir("/PlotObjTypes")==NULL)
    return(NULL);

  if (size<sizeof(PLOTOBJTYPE))
    return (NULL);

  /* allocate structure */
  pot = (PLOTOBJTYPE*) MakeEnvItem (PlotObjTypeName,thePlotObjTypesVarID,size);
  if (pot==NULL)
    return (NULL);

  return (pot);
}

/****************************************************************************/
/*D
   CheckViewPoint - check view point

   SYNOPSIS:
   static INT CheckViewPoint (VIEWEDOBJ *theViewedObj, INT adjust, INT *viewpointcorrect);

   PARAMETERS:
   .  theViewedObj - checks that theViewedObj
   .  adjust - YES if it should be adjusted in the case of `not correct`
   .  viewpointcorrect - YES or NO

   DESCRIPTION:
   Checks if the ViewedObj is adjusted correctly, i.e. the sphere around the 'PLOTOBJ'
   lies completely in the front-half-space of the observer. If 'adjust' equals YES, the function
   tries to move the observer backward in order to find a correct position

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

static INT      CheckViewPoint (VIEWEDOBJ *theViewedObj, INT adjust, INT *viewpointcorrect)
{
  PLOTOBJ *thePlotObj;
  DOUBLE ViewDirection[3], scalarPrd, help[3];

  if (theViewedObj == NULL) return (1);
  if (adjust!=YES && adjust!=NO) return (1);
  if (VO_DIM(theViewedObj)!=TYPE_3D) return (1);
  thePlotObj = VO_PO(theViewedObj);
  *viewpointcorrect = YES;

  /* check position of viewpoint */
  V3_SUBTRACT(VO_VP(theViewedObj),VO_VT(theViewedObj),ViewDirection)
  if (V3_Normalize(ViewDirection))
  {
    UserWrite("ViewPoint and ViewTarget are identical\n");
    *viewpointcorrect = NO;
    return (0);
  }
  V3_LINCOMB(1.0,PO_MIDPOINT(thePlotObj),PO_RADIUS(thePlotObj),ViewDirection,help)
  V3_SUBTRACT(VO_VP(theViewedObj),help,help)
  V3_SCALAR_PRODUCT(ViewDirection,help,scalarPrd)
  if (scalarPrd <= SMALL_C)
  {
    UserWrite("parts of the object lies behind the observer\n");
    VO_STATUS(theViewedObj) = NOT_ACTIVE;
    *viewpointcorrect = NO;
    if (adjust==YES)
    {
      UserWrite("viewpoint has been adjusted\n");
      V3_LINCOMB(1.0,VO_VP(theViewedObj),SMALL_C-scalarPrd,ViewDirection,VO_VP(theViewedObj))
      *viewpointcorrect = YES;
    }
  }

  return (0);
}

/****************************************************************************/
/*D
   SetCutPlane - Initialization/change of cut plane (3D-View)

   SYNOPSIS:
   static INT SetCutPlane (CUT *theCut, INT RemoveCut, const DOUBLE *cutPoint, const DOUBLE *cutNormal)

   PARAMETERS:
   .  theCut - the structure 'CUT' to be initialized
   .  RemoveCut - remove (previously defined) cut
   .  cutPoint - 3-vector for point in plane (NULL if not changed)
   .  cutNormal - 3-vector for normal to plane (NULL if not changed)

   DESCRIPTION:
   This function initializes or changes the definition of the cut plane (only 3D-View).

   If the cutplane is not initialized, both options have to be specified together,
   if it is initialized, specifiing one option is possible to change the cutplane.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

static INT SetCutPlane (CUT *theCut, INT RemoveCut, const DOUBLE *cutPoint, const DOUBLE *cutNormal)
{
  INT popt, nopt;

  if (!(RemoveCut || (cutPoint!=NULL) || (cutNormal!=NULL)))
    /* nothing to do */
    return (0);

  if (RemoveCut)
  {
    CUT_STATUS(theCut) = NOT_INIT;
    return (0);
  }

  /* check if initialized */
  popt = nopt = 0;
  if (CUT_STATUS(theCut) != NOT_INIT)
    popt = nopt = 1;

  /* cut plane point option */
  if (cutPoint!=NULL)
  {
    popt = 1;
    V3_COPY(cutPoint,CUT_PP(theCut));
  }

  /* cut normal direction option */
  if (cutNormal!=NULL)
  {
    nopt = 1;
    V3_COPY(cutNormal,CUT_PN(theCut));
  }

  if (CUT_STATUS(theCut)==NOT_INIT)
    if (!(popt && nopt))
    {
      CUT_STATUS(theCut) = NOT_INIT;
      PrintErrorMessage('W',"SetCutPlane","for initializing cut define plane point AND normal\n");
      return (0);
    }

  /* check how cut plane can now be (re)defined */
  CUT_STATUS(theCut) = NOT_INIT;
  if (popt && nopt)
  {
    if (V3_ISZERO(CUT_PN(theCut)))
    {
      PrintErrorMessage('W',"SetCutPlane","cutting normal is (nearly) zero\n");
      CUT_STATUS(theCut) = NOT_ACTIVE;
    }
    else
      CUT_STATUS(theCut) = ACTIVE;
  }

  return (0);
}

/****************************************************************************/
/*D
   DisplayCutPlane - Display specification of the cut plane (3D)

   SYNOPSIS:
   static INT DisplayCutPlane (CUT *theCut);

   PARAMETERS:
   .  theCut - the 'CUT'plane to be displayed

   DESCRIPTION:
   This function displays specification of cut plane (only 3D). It is used
   by the function by the display-functions for the special content of the
   'PLOTOBJ'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

static INT DisplayCutPlane (const CUT *theCut)
{
  UserWrite("\n");

  /* display content */
  switch (CUT_STATUS(theCut))
  {
  case NOT_INIT :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"CUT STATUS","NOT_INIT");
    return (0);
  case NOT_ACTIVE :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"CUT STATUS","NOT_ACTIVE");
    break;
  case ACTIVE :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"CUT STATUS","ACTIVE");
    break;
  }
  UserWriteF(DISPLAY_PO_FORMAT_SFFF,"PlanePoint",(float)CUT_PP(theCut)[0],(float)CUT_PP(theCut)[1],(float)CUT_PP(theCut)[2]);
  UserWriteF(DISPLAY_PO_FORMAT_SFFF,"PlaneNormal",(float)CUT_PN(theCut)[0],(float)CUT_PN(theCut)[1],(float)CUT_PN(theCut)[2]);

  return (0);
}

/****************************************************************************/
/*D
   SetView - Set the view

   SYNOPSIS:
   INT SetView (PICTURE *thePicture, const DOUBLE *viewPoint, const DOUBLE *targetPoint,
                                const DOUBLE *xAxis, const INT *perspective,
                                INT RemoveCut, const DOUBLE *cutPoint, const DOUBLE *cutNormal);

   PARAMETERS:
   .  thePicture - set view of that picture
   .  viewPoint - new view point
   .  targetPoint - new target point
   .  xAxis - new xAxis
   .  perspective - change of perspective
   .  RemoveCut - remove (previously defined) cut
   .  cutPoint - point in cutting plane
   .  cutNormal - normal direction to cutting plane

   DESCRIPTION:
   This function initializes or changes the view.

   .  2D-view - The function can provide a full initialization from default-values.
   'viewPoint' and 'perspective' have to be NULL, they have no sense in 2D. If the
   'VIEWEDOBJ' is not initialized, the 'xAxis' resp. 'perspective' can be specified,
   if it is not, the default-value `phys. x-axis`(with the size of the 'PROJECTIONPLANE'
   chosen big enough to display the hole 'DOMAIN') resp. 'perspective view' is chosen.
   If the 'VIEWEDOBJ' is initialized, the values that are specified are taken and the old
   are kept.
   .  3D-view - The function can provide a full initialization from default-values, apart of
   the 'viewPoint', which has to be specified to initialize the 'VIEWEDOBJ'. As before
   if for initialization, some parameters are not specified (i.e. ptr==NULL) default values
   are taken: for the 'targetPoint' the midpoint of the 'DOMAIN', for the 'xAxis' a linearcombination
   of the `physical x-axis` and the viewdirection (i.e. to line from target to viewPoint)
   is chosen which is perpendicular to the viewdirection. The y-axis is chosen to build up
   a right-handed system of 'xAxis', 'yAxis' and viewdirection. As a default, the 'perspective'
   is chosen to be 'PERSPECTIVE'. If the view is initialized, new values are taken and old
   values are overtaken if they are not specified.
   .  Remark - If the 'xAxis' is specified, the length of the vector defines the distance from
   the midpoint of the 'PROJECTIONPLANE' to the right side, the hight of the 'PROJECTIONPLANE'
   is allways chosen that the image of the 'PLOTOBJ' on the 'PICTURE' is not distorted, since
   the measures of the 'PICTURE' have to be specified before, this is allways possible.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT SetView (PICTURE *thePicture, const DOUBLE *viewPoint, const DOUBLE *targetPoint, const DOUBLE *xAxis, const INT *perspective,
             INT RemoveCut, const DOUBLE *cutPoint, const DOUBLE *cutNormal, DOUBLE *scale)
{
  VIEWEDOBJ *theViewedObj;
  PLOTOBJ *thePlotObj;
  DOUBLE DefaultVP[3], DefaultVT[3], DefaultVTOld[3], DefaultPMP[3],
         DefaultPXD[3], DefaultPYD[3], DefaultPJ, DefaultSXD[3], DefaultSYD[3], DefaultScale[3];
  DOUBLE ViewDirection[3], ViewDirectionOld[3], CanvasRatio, RotationAxis[3];
  DOUBLE angle, norm;
  INT ViewedObjNotInit, viewpointcorrect;


  /* basics */
  if (thePicture == NULL) return (1);

  /* some inits */
  theViewedObj = PIC_VO(thePicture);
  ViewedObjNotInit = (VO_STATUS(theViewedObj)==NOT_INIT);
  VO_STATUS(theViewedObj) = NOT_INIT;
  thePlotObj = VO_PO(theViewedObj);
  if (PO_STATUS(thePlotObj) == NOT_INIT)
  {
    UserWrite("specify object first\n");
    return (0);
  }
  CanvasRatio = ABS(((DOUBLE)(PIC_GLL(thePicture)[1]-PIC_GUR(thePicture)[1]))
                    / ((DOUBLE)(PIC_GLL(thePicture)[0]-PIC_GUR(thePicture)[0])));

  /* set values */
  switch (PO_DIM(thePlotObj))
  {
  case TYPE_2D :
    if (viewPoint != NULL || perspective != NULL) return (1);

    if (ViewedObjNotInit)
    {
      /* set default values */
      V2_COPY(PO_MIDPOINT(thePlotObj),DefaultVT)
      V2_COPY(PO_MIDPOINT(thePlotObj),DefaultPMP)
      V2_COPY(ex,DefaultPXD)
      V2_COPY(ey,DefaultPYD)
      if (CanvasRatio >= 1.0)
      {
        V2_SCALE(PO_RADIUS(thePlotObj),DefaultPXD)
        V2_SCALE(PO_RADIUS(thePlotObj)*CanvasRatio,DefaultPYD)
      }
      else
      {
        V2_SCALE(PO_RADIUS(thePlotObj)/CanvasRatio,DefaultPXD)
        V2_SCALE(PO_RADIUS(thePlotObj),DefaultPYD)
      }
      V2_COPY(ex,VO_SXD(theViewedObj))
      V2_COPY(ey,VO_SYD(theViewedObj))
      DefaultScale[0] = DefaultScale[1] = 1.0;
    }
    else
    {
      /* save initial values 2D */
      V2_COPY(VO_VT(theViewedObj),DefaultVT)
      V2_COPY(VO_PMP(theViewedObj),DefaultPMP)
      V2_COPY(VO_PXD(theViewedObj),DefaultPXD)
      V2_COPY(VO_PYD(theViewedObj),DefaultPYD)
      V2_COPY(VO_SXD(theViewedObj),DefaultSXD)
      V2_COPY(VO_SYD(theViewedObj),DefaultSYD)
      V2_COPY(VO_SCALE(theViewedObj),DefaultScale)
    }

    /* modify values */
    if (targetPoint != NULL) V2_COPY(targetPoint,DefaultVT)                                                                     /* modify values 2D                             */
      if (xAxis != NULL)
      {
        V2_COPY(xAxis,DefaultPXD)
        V2_COPY(DefaultPXD,DefaultPYD)
        V2_Rotate(DefaultPYD,PI/2.0);
        V2_SCALE(CanvasRatio,DefaultPYD)
      }
    if (scale!=NULL)
    {
      V2_COPY(scale,DefaultScale);
    }

    /* save values */
    V2_COPY(DefaultVT,VO_VT(theViewedObj))
    V2_COPY(DefaultVT,VO_PMP(theViewedObj))
    /*V2_COPY(DefaultPMP,VO_PMP(theViewedObj))*/
    V2_COPY(DefaultPXD,VO_PXD(theViewedObj))
    V2_COPY(DefaultPYD,VO_PYD(theViewedObj))
    V2_COPY(DefaultScale,VO_SCALE(theViewedObj))

    /* set status of ViewedObj */
    VO_STATUS(theViewedObj) = ACTIVE;
    if (VO_PXD(theViewedObj)[0]==0.0 && VO_PXD(theViewedObj)[1]==0.0)
      VO_STATUS(theViewedObj) = NOT_ACTIVE;
    break;

  case TYPE_3D :
    if (ViewedObjNotInit)
    {
      /* set default values 3D */
      if (viewPoint==NULL)
      {
        DefaultVP[_X_] = 1.;
        DefaultVP[_Y_] = 2.;
        DefaultVP[_Z_] = 4.;
        V3_SCALE(PO_RADIUS(thePlotObj),DefaultVP)
        V3_ADD(DefaultVP,PO_MIDPOINT(thePlotObj),DefaultVP)
      }
      else
      {
        V3_COPY(viewPoint,DefaultVP)
      }
      V3_COPY(PO_MIDPOINT(thePlotObj),DefaultVT)
      V3_COPY(PO_MIDPOINT(thePlotObj),DefaultPMP)
      V3_SUBTRACT(DefaultVP,DefaultVT,ViewDirection)
      V3_Orthogonalize(ex,ViewDirection,DefaultPXD);
      if (V3_Normalize(DefaultPXD))
      {
        V3_Orthogonalize(ey,ViewDirection,DefaultPXD);
        if (V3_Normalize(DefaultPXD)) return (1);
      }
      V3_VECTOR_PRODUCT(ViewDirection,DefaultPXD,DefaultPYD)
      if (V3_Normalize(DefaultPYD))
      {
        /* ViewDirection is zero */
        V3_COPY(ey,DefaultPYD);
      }
      if (CanvasRatio >= 1.0)
      {
        V3_SCALE(PO_RADIUS(thePlotObj),DefaultPXD)
        V3_SCALE(PO_RADIUS(thePlotObj)*CanvasRatio,DefaultPYD)
      }
      else
      {
        V3_SCALE(PO_RADIUS(thePlotObj)/CanvasRatio,DefaultPXD)
        V3_SCALE(PO_RADIUS(thePlotObj),DefaultPYD)
      }
      DefaultPJ = YES;

      if (PO_USESCUT(thePlotObj))
        CUT_STATUS(VO_CUT(theViewedObj))                        = NOT_INIT;
    }
    else
    {
      V3_COPY(VO_VP(theViewedObj),DefaultVP)                                                                                                            /* save initial values 3D			*/
      V3_COPY(VO_VT(theViewedObj),DefaultVT)
      V3_COPY(VO_PMP(theViewedObj),DefaultPMP)
      V3_COPY(VO_PXD(theViewedObj),DefaultPXD)
      V3_COPY(VO_PYD(theViewedObj),DefaultPYD)
      DefaultPJ = VO_PERSPECTIVE(theViewedObj);
    }

    /* modify values */
    V3_SUBTRACT(DefaultVP,DefaultVT,ViewDirectionOld)                                                                                           /* modify view-and target point         */
    V3_COPY(DefaultVT,DefaultVTOld)
    if (viewPoint != NULL) V3_COPY(viewPoint,DefaultVP)
      if (targetPoint != NULL) V3_COPY(targetPoint,DefaultVT)
        V3_SUBTRACT(DefaultVP,DefaultVT,ViewDirection)
        V3_VECTOR_PRODUCT(ViewDirectionOld,ViewDirection,RotationAxis)
        if (V3_Normalize(RotationAxis)) V3_COPY(ex,RotationAxis)
          if (V3_Angle(ViewDirectionOld,ViewDirection,&angle)) return (1);
    if (V3_Rotate(DefaultPXD,RotationAxis,angle)) return (1);
    if (V3_Rotate(DefaultPYD,RotationAxis,angle)) return (1);
    V3_SUBTRACT(DefaultPMP,DefaultVTOld,DefaultPMP)
    if (V3_Rotate(DefaultPMP,RotationAxis,angle)) return (1);
    V3_ADD(DefaultPMP,DefaultVT,DefaultPMP)
    if (xAxis != NULL)
    {
      V3_Orthogonalize(xAxis,ViewDirection,DefaultPXD);                                                                                         /* modify xAxis of viewplane		*/
      V3_EUKLIDNORM(DefaultPXD,norm)
      V3_VECTOR_PRODUCT(ViewDirection,DefaultPXD,DefaultPYD)
      if (V3_Normalize(DefaultPYD)) return (1);
      V3_SCALE(CanvasRatio*norm,DefaultPYD)
    }
    if (perspective!=NULL)
      DefaultPJ = *perspective;

    /* save values */
    V3_COPY(DefaultVP,VO_VP(theViewedObj))
    V3_COPY(DefaultVT,VO_VT(theViewedObj))
    V3_COPY(DefaultPMP,VO_PMP(theViewedObj))
    V3_COPY(DefaultPXD,VO_PXD(theViewedObj))
    V3_COPY(DefaultPYD,VO_PYD(theViewedObj))
    VO_PERSPECTIVE(theViewedObj) = DefaultPJ;

    /* set status of ViewedObj */
    VO_STATUS(theViewedObj) = ACTIVE;
    if (V3_ISZERO(VO_PXD(theViewedObj)))
      VO_STATUS(theViewedObj) = NOT_ACTIVE;
    if (CheckViewPoint(theViewedObj,NO,&viewpointcorrect)) return (1);
    if (viewpointcorrect==NO)
      VO_STATUS(theViewedObj) = NOT_ACTIVE;

    if (PO_USESCUT(thePlotObj))
    {
      /* (re)define cut plane */
      if (SetCutPlane(VO_CUT(theViewedObj),RemoveCut,cutPoint,cutNormal)) return (1);

      /*if (CUT_STATUS(VO_CUT(theViewedObj))==NOT_ACTIVE)
              CUT_STATUS(VO_CUT(theViewedObj)) = NOT_INIT;*/

    }
    break;

  case NOT_DEFINED :
    break;

  default :
    return (1);
  }

  VO_STATUS(theViewedObj) = MIN(PO_STATUS(thePlotObj),VO_STATUS(theViewedObj));
  switch (VO_STATUS(theViewedObj))
  {
  case ACTIVE :
    return (0);
  case NOT_ACTIVE :
    UserWrite("viewed object is NOT_ACTIVE\n");
    return (0);
  case NOT_INIT :
    UserWrite("viewed object is NOT_INIT\n");
    return (0);
  default :
    return (1);
  }
}

INT CopyView (const PICTURE *mypic, INT all, INT cut)
{
  UGWINDOW *myUGW,*theUGW;
  MULTIGRID *myMG;
  PICTURE *pic;
  const VIEWEDOBJ *myvo;
  VIEWEDOBJ *vo;
  const PLOTOBJ *thePlotObj;
  INT myViewDim;

  /* basics */
  if (mypic == NULL) return (1);

  /* some inits */
  myvo = PIC_VO(mypic);
  if (VO_STATUS(myvo)!=ACTIVE)
  {
    UserWrite("view is not active\n");
    return (0);
  }
  thePlotObj = VO_PO(myvo);
  if (PO_STATUS(thePlotObj) == NOT_INIT)
  {
    UserWrite("specify object first\n");
    return (0);
  }
  myUGW = PIC_UGW(mypic);
  myMG  = PO_MG(PIC_PO(mypic));
  myViewDim = PO_DIM(thePlotObj);

  if (cut && !PO_USESCUT(VO_PO(myvo)))
    cut = FALSE;

  for (theUGW=GetFirstUgWindow(); theUGW!=NULL; theUGW=GetNextUgWindow(theUGW))
  {
    if (!all)
      theUGW = myUGW;
    for (pic=GetFirstPicture(theUGW); pic!=NULL; pic=GetNextPicture(pic))
    {
      if (pic==mypic) continue;

      vo = PIC_VO(pic);
      if ((PO_DIM(VO_PO(vo))==myViewDim) && (PO_MG(PIC_PO(pic))==myMG))
      {
        /* same multigrid and same dimension */
        switch (myViewDim)
        {
        case TYPE_2D :
          V2_COPY(VO_VT(myvo),VO_VT(vo))
          V2_COPY(VO_PMP(myvo),VO_PMP(vo))
          V2_COPY(VO_PXD(myvo),VO_PXD(vo))
          V2_COPY(VO_PYD(myvo),VO_PYD(vo))
          V2_COPY(VO_SCALE(myvo),VO_SCALE(vo))
          break;

        case TYPE_3D :
          V3_COPY(VO_VP(myvo),VO_VP(vo))
          V3_COPY(VO_VT(myvo),VO_VT(vo))
          V3_COPY(VO_PMP(myvo),VO_PMP(vo))
          V3_COPY(VO_PXD(myvo),VO_PXD(vo))
          V3_COPY(VO_PYD(myvo),VO_PYD(vo))
          VO_PERSPECTIVE(vo) = VO_PERSPECTIVE(myvo);
          if (cut && PO_USESCUT(VO_PO(vo)))
            if (SetCutPlane(VO_CUT(vo),NO,CUT_PP(VO_CUT(myvo)),CUT_PN(VO_CUT(myvo)))) return (1);
          break;
        }
        VO_STATUS(vo) = ACTIVE;
        PIC_VALID(pic) = NO;
      }
    }
    if (!all)
      break;
  }

  return (0);
}

/****************************************************************************/
/*D
   PrintViewSettings - print current view settings of a picture

   SYNOPSIS:
   INT PrintViewSettings (const PICTURE *thePicture)

   PARAMETERS:
   .  thePicture - display the view of that picture

   DESCRIPTION:
   This function displays the view ('viewPoint', 'targetPoint', 'perspective', 'xAxis') stored
   in the 'VIEWEDOBJ' as setview command which would yield the same result.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if not active
   D*/
/****************************************************************************/

INT PrintViewSettings (const PICTURE *thePicture)
{
  const VIEWEDOBJ *theViewedObj;
  const DOUBLE *obs,*tgt,*pxd,*pp,*pn;

  theViewedObj = PIC_VO(thePicture);

  if (VO_STATUS(theViewedObj)!=ACTIVE)
  {
    UserWrite("plotobject not active\n");
    return (1);
  }

  tgt = VO_VT(theViewedObj);
  pxd = VO_PXD(theViewedObj);
  obs = VO_VP(theViewedObj);

  switch (PO_DIM(PIC_PO(thePicture)))
  {
  case TYPE_2D :
    UserWriteF("setview $i $t %g %g $x %g %g\n",
               tgt[_X_],tgt[_Y_],
               pxd[_X_],pxd[_Y_]);
    break;

  case TYPE_3D :
    UserWriteF("setview $i\n\t\t$o %g %g %g\n\t\t$t %g %g %g\n\t\t$x %g %g %g\n\t\t$p %c",
               obs[_X_],obs[_Y_],obs[_Z_],
               tgt[_X_],tgt[_Y_],tgt[_Z_],
               pxd[_X_],pxd[_Y_],pxd[_Z_],
               (VO_PERSPECTIVE(theViewedObj)) ? '<' : '=');
    if (PO_USESCUT(VO_PO(theViewedObj)) && (CUT_STATUS(VO_CUT(theViewedObj))==ACTIVE))
    {
      pp = CUT_PP(VO_CUT(theViewedObj));
      pn = CUT_PN(VO_CUT(theViewedObj));
      UserWriteF("\n\t\t$P %g %g %g\n\t\t$N %g %g %g",
                 pp[_X_],pp[_Y_],pp[_Z_],
                 pn[_X_],pn[_Y_],pn[_Z_]);
    }
    UserWrite(";\n");
    break;
  }

  return (0);
}

/****************************************************************************/
/*D
   DisplayViewOfViewedObject - Display the view

   SYNOPSIS:
   INT DisplayViewOfViewedObject (const PICTURE *thePicture);

   PARAMETERS:
   .  thePicture - display the view of that picture

   DESCRIPTION:
   This function displays the view ('viewPoint', 'targetPoint', 'perspective', 'xAxis') stored
   in the 'VIEWEDOBJ'.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured
   D*/
/****************************************************************************/

INT DisplayViewOfViewedObject (const PICTURE *thePicture)
{
  DOUBLE width;

  UserWrite("-----------------------\n");
  UserWrite(" Display of View of VO \n");
  UserWrite("-----------------------\n");

  switch (VO_STATUS(PIC_VO(thePicture)))
  {
  case NOT_INIT :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"VO_STATUS","NOT_INIT");
    return (0);
    break;
  case NOT_ACTIVE :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"VO_STATUS","NOT_ACTIVE");
    break;
  case ACTIVE :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"VO_STATUS","ACTIVE");
    break;
  default :
    return (1);
  }

  switch (VO_DIM(PIC_VO(thePicture)))
  {
  case NOT_DEFINED :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"Dim","NOT_DEFINED");
    break;
  case TYPE_2D :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"Dim","TYPE_2D");
    UserWriteF(DISPLAY_PO_FORMAT_SFF,"Target",(float)(VO_VT(PIC_VO(thePicture))[0]),(float)(VO_VT(PIC_VO(thePicture))[1]));
    V2_EUKLIDNORM(VO_PXD(PIC_VO(thePicture)),width)
    width *= 2.0;
    UserWriteF(DISPLAY_PO_FORMAT_SF,"WinWidth",(float)width);
    break;
  case TYPE_3D :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"Dim","TYPE_3D");
    UserWriteF(DISPLAY_PO_FORMAT_SFFF,"Observer",(float)(VO_VP(PIC_VO(thePicture))[0]),(float)(VO_VP(PIC_VO(thePicture))[1]),(float)(VO_VP(PIC_VO(thePicture))[2]));
    UserWriteF(DISPLAY_PO_FORMAT_SFFF,"Target",(float)(VO_VT(PIC_VO(thePicture))[0]),(float)(VO_VT(PIC_VO(thePicture))[1]),(float)(VO_VT(PIC_VO(thePicture))[2]));
    V3_EUKLIDNORM(VO_PXD(PIC_VO(thePicture)),width)
    width *= 2.0;
    UserWriteF(DISPLAY_PO_FORMAT_SF,"WinWidth",(float)width);

    /* print content of cut plane */
    if (PO_USESCUT(PIC_PO(thePicture)))
      if (DisplayCutPlane(VO_CUT(PIC_VO(thePicture)))) return (1);
    break;
  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*D
   ResizeViewPlane - Resize physical view plane

   SYNOPSIS:
   INT ResizeViewPlane (VIEWEDOBJ *theVO, const INT *Pix_LL_old,
   const INT *Pix_UR_old, const INT *Pix_LL_new, const INT *Pix_UR_new);

   PARAMETERS:
   .  theVO - the ViewedObject
   .  'Pix_LL_old',~'Pix_UR_old' - size of the 'PICTURE' before resizing (pixel coordinates
   on the 'UGWINDOW')
   .  'Pix_LL_new',~'Pix_UR_new' - size of the 'PICTURE' after resizing (pixel coordinates
   on the 'UGWINDOW')

   DESCRIPTION:
   This function resizes physical view plane (i.e. 'xAxis', 'yAxis') to the change of the
   pix-coordinate in a way that the scale of the picture remains the same.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT ResizeViewPlane (VIEWEDOBJ *theVO, const INT *Pix_LL_old, const INT *Pix_UR_old, const INT *Pix_LL_new, const INT *Pix_UR_new)
{
  DOUBLE delta[3],q[2];

  if (VO_STATUS(theVO) == NOT_INIT) return (0);

  q[0] = 1.0/(DOUBLE)(Pix_UR_old[0]-Pix_LL_old[0]);
  q[1] = 1.0/(DOUBLE)(Pix_UR_old[1]-Pix_LL_old[1]);

  switch (VO_DIM(theVO))
  {
  case TYPE_2D :
    V2_LINCOMB(q[0]*(Pix_UR_new[0]-Pix_UR_old[0]+Pix_LL_new[0]-Pix_LL_old[0]),VO_PXD(theVO),
               q[1]*(Pix_UR_new[1]-Pix_UR_old[1]+Pix_LL_new[1]-Pix_LL_old[1]),VO_PYD(theVO),
               delta)
    V2_ADD(VO_PMP(theVO),delta,VO_PMP(theVO))
    V2_SCALE(q[0]*(Pix_UR_new[0]-Pix_LL_new[0]),VO_PXD(theVO))
    V2_SCALE(q[1]*(Pix_UR_new[1]-Pix_LL_new[1]),VO_PYD(theVO))
    break;
  case TYPE_3D :
    V3_LINCOMB(q[0]*(Pix_UR_new[0]-Pix_UR_old[0]+Pix_LL_new[0]-Pix_LL_old[0]),VO_PXD(theVO),
               q[1]*(Pix_UR_new[1]-Pix_UR_old[1]+Pix_LL_new[1]-Pix_LL_old[1]),VO_PYD(theVO),
               delta)
    V3_ADD(VO_PMP(theVO),delta,VO_PMP(theVO))
    V3_SCALE(q[0]*(Pix_UR_new[0]-Pix_LL_new[0]),VO_PXD(theVO))
    V3_SCALE(q[1]*(Pix_UR_new[1]-Pix_LL_new[1]),VO_PYD(theVO))
    break;
  default :
    return (1);
  }

  return (0);
}


/****************************************************************************/
/*D
   Walk	- Modify the view by walking

   SYNOPSIS:
   INT Walk (PICTURE *thePicture, const DOUBLE *vrsDelta);

   PARAMETERS:
   .  thePicture - the view of this 'PICTURE' will be changed
   .  vrsDelta - desired change in the viewpoint

   DESCRIPTION:
   This function modifies the view by walking, it adds to the viewPoint the 'vrsDelta'
   and changes according to this movement the 'xAxis' and 'yAxis'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT Walk (PICTURE *thePicture, const DOUBLE *vrsDelta)
{
  VIEWEDOBJ *theViewedObj;
  DOUBLE VP[3], XD[3], YD[3], ZD[3];

  /* basics */
  if (thePicture == NULL || vrsDelta==NULL) return (1);
  theViewedObj = PIC_VO(thePicture);
  if (VO_STATUS(theViewedObj)==NOT_INIT)
  {
    UserWrite("status of view: NOT_INIT\n");
    return (0);
  }

  /* walk */
  switch (VO_DIM(theViewedObj))
  {
  case TYPE_2D :
    V2_COPY(VO_PXD(theViewedObj),XD)
    if (V2_Normalize(XD)) return (1);
    V2_COPY(VO_PYD(theViewedObj),XD)
    if (V2_Normalize(YD)) return (1);
    V2_LINCOMB(vrsDelta[0],XD,vrsDelta[1],YD,VP)
    V2_ADD(VO_VP(theViewedObj),vrsDelta,VP)
    break;
  case TYPE_3D :
    V3_COPY(VO_PXD(theViewedObj),XD)
    if (V3_Normalize(XD)) return (1);
    V3_COPY(VO_PYD(theViewedObj),YD)
    if (V3_Normalize(YD)) return (1);
    V3_VECTOR_PRODUCT(YD,XD,ZD)
    V3_LINCOMB(vrsDelta[0],XD,vrsDelta[1],YD,VP)
    V3_LINCOMB(1.0,VP,vrsDelta[2],ZD,VP)
    V3_ADD(VO_VP(theViewedObj),VP,VP)
    break;
  default :
    return (1);
  }
  if (SetView(thePicture,VP,NULL,NULL,NULL,NO,NULL,NULL,NULL)) return (1);

  return (0);
}

/****************************************************************************/
/*D
   RunAroundTargetPoint	- Modify the view by running around the midpoint

   SYNOPSIS:
   INT RunAroundTargetPoint (PICTURE *thePicture, DOUBLE vrsDirectionAngle,
   DOUBLE vrsAngle);

   PARAMETERS:
   .  thePicture - the view of this 'PICTURE' will be changed
   .  vrsDirectionAngle - determines direction in which to run
   .  vrsAngle - run this angle

   DESCRIPTION:
   This function modifies the view by running around the midpoint. the 'vrsDirectionAngle'
   determines the direction in which to run around the targetPoint, it is taken w.r.t.
   the xAxis of the viewplane (0 degrees) and mathmatically positive, i.e. going in direction
   of yAxis means 'vrsDirectionAngle'=90 degrees, in the minus xAxis means 180  degrees etc.
   In this direction you run around the 'targetPoint' by an angle of 'vrsAngle'. The
   'xAxis' and 'yAxis' then are adjusted.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT RunAroundTargetPoint (PICTURE *thePicture, DOUBLE vrsDirectionAngle, DOUBLE vrsAngle)
{
  VIEWEDOBJ *theViewedObj;
  DOUBLE VP[3], RotationAxis[3], ViewDirection[3], TurnDirection[3];

  /* basics */
  if (thePicture == NULL) return (1);
  theViewedObj = PIC_VO(thePicture);
  if (VO_DIM(theViewedObj) != TYPE_3D)
  {
    UserWrite("dimension of view is not 3D\n");
    return (0);
  }

  /* run around */
  V3_SUBTRACT(VO_VP(theViewedObj),VO_VT(theViewedObj),ViewDirection)
  V3_COPY(VO_PXD(theViewedObj),TurnDirection)
  if (V3_Rotate(TurnDirection,ViewDirection,vrsDirectionAngle))
  {
    UserWrite("cannot run around target\n");
    return (0);
  }
  V3_VECTOR_PRODUCT(ViewDirection,TurnDirection,RotationAxis)
  if (V3_Rotate(ViewDirection,RotationAxis,vrsAngle))
  {
    UserWrite("cannot run around target\n");
    return (0);
  }
  V3_ADD(VO_VT(theViewedObj),ViewDirection,VP)
  if (SetView(thePicture,VP,NULL,NULL,NULL,NO,NULL,NULL,NULL)) return (1);

  return (0);
}

/****************************************************************************/
/*D
   Zoom	- zoom the view, i.e.: change size of projection plane

   SYNOPSIS:
   INT Zoom (PICTURE *thePicture, DOUBLE factor);

   PARAMETERS:
   .  thePicture - change size of projection plane of that 'PICTURE'
   .  factor - zoom factor

   DESCRIPTION:
   This function zooms the view, i.e.: change size of projection plane. A factor
   0.5 shrinks the 'PROJECTOINPLANE' by an factor of two, the image will be more
   detailled.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT Zoom (PICTURE *thePicture, DOUBLE factor)
{
  VIEWEDOBJ *theViewedObj;

  /* basics */
  if (thePicture == NULL) return (1);
  theViewedObj = PIC_VO(thePicture);
  if (VO_STATUS(theViewedObj)==NOT_INIT)
  {
    UserWrite("status of view: NOT_INIT\n");
    return (0);
  }
  if (factor <= 0.0)
  {
    UserWrite("zoom factor has to be positve\n");
    return (0);
  }

  /* zoom */
  switch (VO_DIM(theViewedObj))
  {
  case TYPE_2D :
    V2_SCALE(factor,VO_PXD(theViewedObj))
    V2_SCALE(factor,VO_PYD(theViewedObj))
    break;
  case TYPE_3D :
    V3_SCALE(factor,VO_PXD(theViewedObj))
    V3_SCALE(factor,VO_PYD(theViewedObj))
    break;
  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*D
   DragProjectionPlane - moves the 'PROJECTIONPLANE'

   SYNOPSIS:
   INT DragProjectionPlane (PICTURE *thePicture, DOUBLE vrsDeltaX,
   DOUBLE vrsDeltaY);

   PARAMETERS:
   .  thePicture - of this 'PICTURE'
   .  vrsDeltaX,~vrsDeltaY - displacement of the 'PROJECTIONPLANE'

   DESCRIPTION:
   This function moves the 'PROJECTIONPLANE' in its 'xAxis' and 'yAxis' direction,
   by adding 'vrsDeltaX' resp. 'vrsDeltaY'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT DragProjectionPlane (PICTURE *thePicture, DOUBLE vrsDeltaX, DOUBLE vrsDeltaY)
{
  VIEWEDOBJ *theViewedObj;
  DOUBLE DragVector[3], help[3];

  /* basics */
  if (thePicture == NULL) return (1);
  theViewedObj = PIC_VO(thePicture);
  if (VO_STATUS(theViewedObj)==NOT_INIT)
  {
    UserWrite("status of view: NOT_INIT\n");
    return (0);
  }

  /* drag */
  switch (VO_DIM(theViewedObj))
  {
  case TYPE_2D :
    V2_COPY(VO_PXD(theViewedObj),DragVector)
    V2_Normalize(DragVector);
    V2_COPY(VO_PYD(theViewedObj),help)
    V2_Normalize(help);
    V2_LINCOMB(vrsDeltaX,DragVector,vrsDeltaY,help,DragVector)
    V2_ADD(VO_PMP(theViewedObj),DragVector,VO_PMP(theViewedObj))
    break;
  case TYPE_3D :
    V3_COPY(VO_PXD(theViewedObj),DragVector)
    V3_Normalize(DragVector);
    V3_COPY(VO_PYD(theViewedObj),help)
    V3_Normalize(help);
    V3_LINCOMB(vrsDeltaX,DragVector,vrsDeltaY,help,DragVector)
    V3_ADD(VO_PMP(theViewedObj),DragVector,VO_PMP(theViewedObj))
    break;
  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*D
   RotateProjectionPlane -

   SYNOPSIS:
   INT RotateProjectionPlane (PICTURE *thePicture, DOUBLE vrsAngle);

   PARAMETERS:
   .  thePicture - change view of that 'PICTURE'
   .  vrsAngle - rotate that angle

   DESCRIPTION:
   This function rotates the 'PROJECTIONPLANE' by an agle 'vrsAngle', around
   the viewdirection (vector from viewPoint to viewTarget) in mathmatical
   positive sense (i.e. rotation by 90 degree moves the new 'xAxis' to be parallel
   to the old 'yAxis', and the new 'yAxis' happens to be parallel to the old
   minus-'xAxis'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT RotateProjectionPlane (PICTURE *thePicture, DOUBLE vrsAngle)
{
  VIEWEDOBJ *theViewedObj;
  DOUBLE ViewDirection[3];

  /* basics */
  if (thePicture == NULL) return (1);
  theViewedObj = PIC_VO(thePicture);
  if (VO_STATUS(theViewedObj)==NOT_INIT)
  {
    UserWrite("status of view: NOT_INIT\n");
    return (0);
  }

  /* drag */
  switch (VO_DIM(theViewedObj))
  {
  case TYPE_2D :
    V2_Rotate(VO_PXD(theViewedObj),vrsAngle);
    V2_Rotate(VO_PYD(theViewedObj),vrsAngle);
    break;
  case TYPE_3D :
    V3_SUBTRACT(VO_VP(theViewedObj),VO_VT(theViewedObj),ViewDirection)
    if (V3_Normalize(ViewDirection))
    {
      UserWrite("cannot rotate Projection plane\n");
      return (0);
    }
    V3_Rotate(VO_PXD(theViewedObj),ViewDirection,vrsAngle);
    V3_Rotate(VO_PYD(theViewedObj),ViewDirection,vrsAngle);
    break;
  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*D
   SpecifyPlotObject -

   SYNOPSIS:
   static INT SpecifyPlotObject (PLOTOBJ *thePlotObj, MULTIGRID *theMG,
   const char *thePlotObjTypeName, INT argc, char **argv);

   PARAMETERS:
   .  thePlotObj - specify this 'PLOTOBJ'
   .  theMG - the pointer to multigrid
   .  thePlotObjTypeName - name of the 'PLOTOBJTYPE'
   .  argc,~argv - parameters coming from 'COMMANDS', tunneled to init routines
                for 'PLOTOBJ' (stored in 'PLOTOBJTYPE')

   DESCRIPTION:
   This function initializes or changes the specification of a 'PLOTOBJ'. This
   procedure has two parts.
   .  General~part - If the 'PLOTOBJ' is not initialized, the 'thePlotObjTypeName'
                  has to be specified to chose the type of plotobject and the
                                  'theMG' has to be specified to specifiy the database.
                                  If the 'PLOTOBJ' is initialized, you can reinitialize the
                                  plotobject by specifiing the 'thePlotObjTypeName', only in this
                                  case the 'theMG' is used.
   .  Specific~part - According to the chosen 'PLOTOBJTYPE' the specific part is
                   initialized or changed (see 'setplotobject')

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

static INT SpecifyPlotObject (PLOTOBJ *thePlotObj, MULTIGRID *theMG, const char *thePlotObjTypeName, INT argc, char **argv)
{
  INT i, ret, def;

  if (thePlotObj==NULL || (theMG==NULL && thePlotObjTypeName!=NULL)) return (1);
  if (PO_STATUS(thePlotObj)==NOT_INIT && thePlotObjTypeName==NULL)
  {
    UserWrite("cannot initialize PlotObject\n");
    return (0);
  }

  /* find objecttype */
  if (thePlotObjTypeName!=NULL)
  {
    PO_POT(thePlotObj) = GetPlotObjType(thePlotObjTypeName);
    if (PO_POT(thePlotObj)==NULL)
    {
      UserWrite("cannot find specified PlotObjectType\n");
      return (0);
    }
    PO_STATUS(thePlotObj)=NOT_INIT;
    PO_MG(thePlotObj) = theMG;
  }

  /* set general information */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    def = YES;
  }
  else
  {
    def = PO_CBD(thePlotObj);
  }
  for (i=1; i<argc; i++)
  {
    if (strcmp(argv[i],"clearOn")==0)
      def = YES;
    if (strcmp(argv[i],"clearOff")==0)
      def = NO;
  }
  PO_CBD(thePlotObj) = def;

  /* init object */
  PO_USESCUT(thePlotObj) = NO;
  ret = (*PO_POT(thePlotObj)->SetPlotObjProc)(thePlotObj,argc,argv);
  switch(ret)
  {
  case ACTIVE :
    PO_STATUS(thePlotObj) = ACTIVE;
    break;
  case NOT_ACTIVE :
    PO_STATUS(thePlotObj) = NOT_ACTIVE;
    UserWrite("plot object is NOT_ACTIVE\n");
    break;
  case NOT_INIT :
    PO_STATUS(thePlotObj) = NOT_INIT;
    PO_POT(thePlotObj) = NULL;
    UserWrite("plot object is NOT_INIT\n");
    break;
  default :
    return (1);
  }

  return (0);
}

/****************************************************************************/
/*D
   SpecifyPlotObjOfViewedObject - specify the 'PLOTOBJ' of 'VIEWEDOBJ'

   SYNOPSIS:
   INT SpecifyPlotObjOfViewedObject (PICTURE *thePicture, MULTIGRID *theMG,
   const char *PlotObjTypeName, INT argc, char **argv);

   PARAMETERS:
   .  thePicture - specify the 'PLOTOBJ' of this 'PICTURE'
   .  theMG - the pointer to multigrid
   .  thePlotObjTypeName - name of the 'PLOTOBJTYPE'
   .  argc,~argv - parameters coming from 'COMMANDS', tunneled to init routines
                for 'PLOTOBJ' (stored in 'PLOTOBJTYPE')

   DESCRIPTION:
   This function initializes the 'PLOTOBJ' of the picture by calling 'SpecifyPlotObject'
   and sets the view. If the 'a'-option from 'COMMANDS' is pecified a 3D view will be adjusted,
   i.e. the observer is moved away from the 'PLOTOBJ' if he/she is in the bounding sphere
   of the 'PLOTOBJ' in order to find a point outside.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT SpecifyPlotObjOfViewedObject (PICTURE *thePicture, MULTIGRID *theMG, const char *PlotObjTypeName, INT argc, char **argv)
{
  VIEWEDOBJ *theViewedObj;
  PLOTOBJ *thePlotObj;
  PLOTOBJTYPE *thePlotObjType;
  INT aopt, i, viewpointcorrect;

  /* basics */
  if (thePicture == NULL) return (1);
  theViewedObj = PIC_VO(thePicture);
  thePlotObj = VO_PO(theViewedObj);
  PO_PIC(thePlotObj) = thePicture;

  /* store */
  thePlotObjType = PO_POT(thePlotObj);

  /* init PlotObject */
  if (SpecifyPlotObject(thePlotObj,theMG,PlotObjTypeName,argc,argv)) return (1);
  VO_STATUS(theViewedObj) = MIN(VO_STATUS(theViewedObj),PO_STATUS(thePlotObj));

  /* check view */
  if (PO_POT(thePlotObj)!=thePlotObjType)
  {
    if (VO_STATUS(theViewedObj) != NOT_INIT)
      UserWrite("PlotObjectType has changed: view is reset now\n");
    VO_STATUS(theViewedObj) = NOT_INIT;
    return (0);
  }
  if (PO_DIM(thePlotObj)==TYPE_3D)
  {
    /* adjust view option */
    aopt = NO;
    for (i=1; i<argc; i++)
      if (argv[i][0]=='a')
      {
        aopt = YES;
        break;
      }

    if (CheckViewPoint(theViewedObj,aopt,&viewpointcorrect)) return (1);
    if (viewpointcorrect && VO_STATUS(theViewedObj)==ACTIVE)
      VO_STATUS(theViewedObj) = ACTIVE;
  }
  if (SetView(thePicture,NULL,NULL,NULL,NULL,NO,NULL,NULL,NULL)) return (1);

  return (0);
}

/****************************************************************************/
/*D
   DisplayObject - Display the description of the object

   SYNOPSIS:
   INT DisplayObject (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj - the 'PLOTOBJ' to display

   DESCRIPTION:
   This function displays the specification of the plotobject.

   .  General~part - Independent of the 'PLOTOBJTYPE'. Shows the name of the
                  'PLOTOBJ', the ame of the 'MULTIGRID', the 'status' of the
                                  'PLOTOBJ' ('INIT', 'NOT_ACTIVE' or 'ACTIVE') and its
                                  bounding sphere.
   .  Specific~part - Depends on the 'PLOTOBJTYPE'. The display routine for the
                   'PLOTOBJ' is used.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT DisplayObject (PLOTOBJ *thePlotObj)
{
  PLOTOBJTYPE *thePOT;

  if (thePlotObj==NULL) return (1);
  thePOT = PO_POT(thePlotObj);
  UserWrite("-----------------------\n");
  UserWrite(" Display of PlotObject \n");
  UserWrite("-----------------------\n");
  switch (PO_STATUS(thePlotObj))
  {
  case NOT_INIT :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"PO-NAME","---");
    UserWriteF(DISPLAY_PO_FORMAT_SS,"MG-NAME","---");
    UserWriteF(DISPLAY_PO_FORMAT_SS,"STATUS","NOT_INIT");
    return (0);
  case NOT_ACTIVE :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"PO-NAME",thePOT->v.name);
    UserWriteF(DISPLAY_PO_FORMAT_SS,"MG-NAME",MGNAME(PO_MG(thePlotObj)));
    if (PO_DIM(thePlotObj) == TYPE_2D)
      UserWriteF(DISPLAY_PO_FORMAT_SS,"STATUS","NOT_ACTIVE:2D");
    else
      UserWriteF(DISPLAY_PO_FORMAT_SS,"STATUS","NOT_ACTIVE:3D");
    break;
  case ACTIVE :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"PO-NAME",thePOT->v.name);
    UserWriteF(DISPLAY_PO_FORMAT_SS,"MG-NAME",MGNAME(PO_MG(thePlotObj)));
    if (PO_DIM(thePlotObj) == TYPE_2D)
      UserWriteF(DISPLAY_PO_FORMAT_SS,"STATUS","ACTIVE:2D");
    else
      UserWriteF(DISPLAY_PO_FORMAT_SS,"STATUS","ACTIVE:3D");
    break;
  }
  UserWriteF(DISPLAY_PO_FORMAT_SS,"CLEAR FIRST",PO_CBD(thePlotObj) ? "YES" : "NO");

  if (thePOT == NULL) return (0);

  switch (PO_DIM(thePlotObj))
  {
  case TYPE_2D :
    UserWriteF(DISPLAY_PO_FORMAT_SFF,"MIDPOINT",(float)PO_MIDPOINT(thePlotObj)[0],(float)PO_MIDPOINT(thePlotObj)[1]);
    UserWriteF(DISPLAY_PO_FORMAT_SF,"RADIUS",(float)PO_RADIUS(thePlotObj));
    break;
  case TYPE_3D :
    UserWriteF(DISPLAY_PO_FORMAT_SFFF,"MIDPOINT",(float)PO_MIDPOINT(thePlotObj)[0],(float)PO_MIDPOINT(thePlotObj)[1],(float)PO_MIDPOINT(thePlotObj)[2]);
    UserWriteF(DISPLAY_PO_FORMAT_SF,"RADIUS",(float)PO_RADIUS(thePlotObj));
    break;
  }
  UserWrite("\n");

  /* display special object description */
  if (PO_POT(thePlotObj)->DispPlotObjProc == NULL) return (1);
  if ((*PO_POT(thePlotObj)->DispPlotObjProc)(thePlotObj)) return (1);

  /* end */
  UserWrite("-----------------------\n");

  return (0);
}

/****************************************************************************/
/*D
   Matrix - plot object to show matrix entries in row/column style

   DESCRIPTION:
   The Matrix plot object shows the stiffness matrix (any mat data desc)
   in row/column style. The infobox shows the (row,column) index of the
   entry the mouse is pointing to. Colors indicate the size of the entries
   ('findrange'!). Zero entries are suppressed. If it is magnified enough the
   value of the entries is printed.

   Mandatory options are:~
   .    $e~<mat~eval~proc>		- either specify <mat eval proc>
   .    $M~<matdata~desc>		- or			 <matdata desc>

   Possible options:~
   .    $f~<from>				- from value
   .    $t~<to>				- to value
   .    $T~<thresh>			- supress entries smaller than <thresh>
   .    $l~0|1					- logarithmic scale off/on
   .    $B~0|1~<dash>~<space>	- block vectors off/on
   .    $r~0|1					- size relative to diagonal entry off/on
   .    $C~0|1					- regular entries off/on
   .    $E~0|1					- extra entries off/on

   KEYWORDS:
   graphics, plot, window, picture, plotobject, matdesc
   D*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  InitMatrixPlotObject											*/
/*																			*/
/* Purpose:   initialization of matrix object								*/
/*																			*/
/* Input:	  PLOTOBJ *thePlotObj, INT argc, char **argv					*/
/*																			*/
/* Output:	  INT NOT_INIT, NOT_ACTIVE or ACTIVE							*/
/*																			*/
/****************************************************************************/

static INT InitMatrixPlotObject (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  struct MatrixPlotObj *theMpo;
  GRID *theGrid;
  INT i;
  float fValue, f2Value;
  int iValue;
  char buffer[64];

  theMpo = &(thePlotObj->theMpo);
  theGrid = GRID_ON_LEVEL(PO_MG(thePlotObj),PO_MG(thePlotObj)->currentLevel);
  if (theGrid == NULL) return (NOT_INIT);
  PO_MIDPOINT(thePlotObj)[0] = PO_MIDPOINT(thePlotObj)[1] = NVEC(theGrid)/2.0;
  PO_RADIUS(thePlotObj) = NVEC(theGrid)/2.0;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theMpo->min                     =-4.0;
    theMpo->max                     = 4.0;
    theMpo->log                     = FALSE;
    theMpo->conn            = TRUE;
    theMpo->extra           = FALSE;
    theMpo->rel                     = FALSE;
    theMpo->EvalFct         = NULL;
    theMpo->Matrix          = NULL;
    theMpo->dash            = 0.0;
    theMpo->space           = 0.0;
  }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'f' :
      if (sscanf(argv[i],"f %f",&fValue)==1)
        theMpo->min = fValue;
      break;

    case 't' :
      if (sscanf(argv[i],"t %f",&fValue)==1)
        theMpo->max = fValue;
      break;

    case 'T' :
      if (sscanf(argv[i],"T %f",&fValue)!=1)
        PrintErrorMessage('E',"Matrix","specify value with T option");
      else
        theMpo->thresh = fValue;
      break;

    case 'l' :
      if (sscanf(argv[i],"l %d",&iValue)!=1)
        break;
      if (iValue==1)
        theMpo->log = YES;
      else if (iValue==0)
        theMpo->log = NO;
      break;

    case 'B' :
      iValue = 0;
      fValue = 0.0;
      f2Value = 0.0;
      if (sscanf(argv[i],"BV %d %f %f",&iValue,&fValue,&f2Value)==0)
        break;
      if (iValue==1)
        theMpo->BV = YES;
      else if (iValue==0)
        theMpo->BV = NO;
      theMpo->dash = fValue;
      theMpo->space = f2Value;
      break;

    case 'r' :
      if (sscanf(argv[i],"r %d",&iValue)!=1)
        break;
      if (iValue==1)
        theMpo->rel = YES;
      else if (iValue==0)
        theMpo->rel = NO;
      break;

    case 'C' :
      if (sscanf(argv[i],"C %d",&iValue)!=1)
        break;
      if (iValue==1)
        theMpo->conn = YES;
      else if (iValue==0)
        theMpo->conn = NO;
      break;

    case 'E' :
      if (sscanf(argv[i],"E %d",&iValue)!=1)
        break;
      if (iValue==1)
        theMpo->extra = YES;
      else if (iValue==0)
        theMpo->extra = NO;
      break;

    case 'e' :
      if (sscanf(argv[i],"e %s",buffer)!=1)
        break;
      theMpo->EvalFct = GetMatrixValueEvalProc(buffer);
      if (theMpo->EvalFct == NULL)
      {
        UserWrite("cannot find plot procedure\n");
        return (NOT_ACTIVE);
      }
      break;

    case 'M' :
      if (sscanf(argv[i],"M %s",buffer)!=1)
        break;
      theMpo->Matrix=GetMatDataDescByName(PO_MG(thePlotObj),buffer);
      if (theMpo->Matrix == NULL)
      {
        UserWrite("cannot find matrix symbol\n");
        return (NOT_ACTIVE);
      }
      break;
    }
  if ((theMpo->EvalFct == NULL) && (theMpo->Matrix == NULL))
  {
    UserWrite("specify a scalar matrix symbol or a matrix plot procedure\n");
    return (NOT_ACTIVE);
  }

  return (ACTIVE);
}

/****************************************************************************/
/*																			*/
/* Function:  DisplayMatrixPlotObject										*/
/*																			*/
/* Purpose:   display content of matrix object								*/
/*																			*/
/* Input:	  PLOTOBJ *thePlotObj											*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT DisplayMatrixPlotObject (PLOTOBJ *thePlotObj)
{
  struct MatrixPlotObj *theMpo;

  theMpo = &(thePlotObj->theMpo);

  /* print range */
  UserWriteF(DISPLAY_PO_FORMAT_SFF,"range",(float)theMpo->min,(float)theMpo->max);
  UserWriteF(DISPLAY_PO_FORMAT_SS,"regular conn.",(theMpo->conn) ? "YES" : "NO");
  UserWriteF(DISPLAY_PO_FORMAT_SS,"extra   conn.",(theMpo->extra) ? "YES" : "NO");
  UserWriteF(DISPLAY_PO_FORMAT_SS,"use log",(theMpo->log) ? "YES" : "NO");
  UserWriteF(DISPLAY_PO_FORMAT_SS,"rel values",(theMpo->rel) ? "YES" : "NO");
  UserWriteF(DISPLAY_PO_FORMAT_SF,"Thresh",(float)theMpo->thresh);
  UserWriteF(DISPLAY_PO_FORMAT_SS,"BV blocks",(theMpo->BV) ? "YES" : "NO");

  /* print procedure name */
  if (theMpo->EvalFct!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc",((ENVVAR*)theMpo->EvalFct)->name);

#ifdef __NP__
  if (theMpo->Matrix!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"Matrix",ENVITEM_NAME(theMpo->Matrix));
#endif

  return (0);
}

#ifdef __TWODIM__

/****************************************************************************/
/*
   InitGridPlotObject_2D - Initialization of 2D grid object

   SYNOPSIS:
   static INT InitGridPlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv);

   PARAMETERS:
   .  thePlotObj - the 'PLOTOBJ' to be initialized.
   .  argc, argv -

   DESCRIPTION:
   This function does initialization of 2D grid object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   Grid - plot object

   please lookup 'Grid2D' or 'Grid3D' depending on the space dimension you are
   working with (the options are not quite the same)
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   Grid2D - plot object for multigrid (omit 2D when using it)

   DESCRIPTION:
   The Grid plot object shows the multigrid with its elements.

   Possible options:~
   .    $c~0|1					- use color off/on
   .    $w~c|i|r|a				- show copy/irregular/regular/all elements
   .    $b~0|1					- show boundary off/on
   .    $r~0|1					- show refinement marks off/on
   .    $i~0|1					- show error indicator marks off/on
                                                                (only with r and c option off)
   .    $e~0|1					- show element IDs off/on
   .    $S~0|1					- show element subdomain IDs off/on
   .    $n~0|1					- show node IDs off/on
   .    $m~0|1					- plot node markers off/on
   .    $s~<shrink>			- factor to shrink elements
   .    $p~<shrink>			- parallel only: factor to shrink processor partition
   .    $free~<vd>				- vec data desc describing new global coordinates of free boundary

   KEYWORDS:
   graphics, plot, window, picture, plotobject, multigrid, elements
   D*/
/****************************************************************************/

static INT InitGridPlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  VECDATA_DESC *vd;
  struct GridPlotObj2D *theGpo;
  char buffer[VALUELEN];

  theGpo = &(thePlotObj->theGpo);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theGpo->ShrinkFactor            = 1.0;
                #ifdef ModelP
    theGpo->PartShrinkFactor        = 1.0;
                #endif
    theGpo->ElemColored             = 1;
    theGpo->WhichElem                       = PO_ALL;
    theGpo->PlotBoundary            = YES;
    theGpo->PlotElemID                      = NO;
    theGpo->PlotNodeID                      = NO;
    theGpo->PlotNodes                       = NO;
    theGpo->PlotRefMarks            = NO;
    theGpo->PlotIndMarks            = NO;
    theGpo->FreeBnd                 = NULL;
  }

  /* scan options */
  if (ReadArgvChar("w",buffer,argc,argv)==0)
    switch (buffer[0])
    {
    case 'c' : theGpo->WhichElem = PO_COPY; break;
    case 'i' : theGpo->WhichElem = PO_IRR; break;
    case 'r' : theGpo->WhichElem = PO_REG; break;
    case 'a' : theGpo->WhichElem = PO_ALL; break;
    default :  return (NOT_ACTIVE);
    }
  ReadArgvDOUBLE("s",&theGpo->ShrinkFactor,       argc,argv);
  ReadArgvINT   ("c",&theGpo->ElemColored,        argc,argv);
  ReadArgvINT   ("b",&theGpo->PlotBoundary,       argc,argv);
  ReadArgvINT   ("r",&theGpo->PlotRefMarks,       argc,argv);
  ReadArgvINT   ("i",&theGpo->PlotIndMarks,       argc,argv);
  ReadArgvINT   ("e",&theGpo->PlotElemID,         argc,argv);
  ReadArgvINT   ("S",&theGpo->PlotSubdomain,      argc,argv);
  ReadArgvINT   ("n",&theGpo->PlotNodeID,         argc,argv);
  ReadArgvINT   ("m",&theGpo->PlotNodes,          argc,argv);

  vd = ReadArgvVecDesc(PO_MG(thePlotObj),"free",argc,argv);
  if (vd!=NULL) theGpo->FreeBnd = vd;
        #ifdef ModelP
  ReadArgvDOUBLE("p",&theGpo->PartShrinkFactor,argc,argv);
  if (theGpo->PartShrinkFactor<=0.0 || theGpo->PartShrinkFactor>1.0)
    return (NOT_ACTIVE);
        #endif

  /* check validity */
  if (theGpo->ShrinkFactor<=0.0 || theGpo->ShrinkFactor>1.0)
    return (NOT_ACTIVE);
  if (theGpo->ElemColored<0 || theGpo->ElemColored>2) return (NOT_ACTIVE);
  if (theGpo->PlotIndMarks == YES)
  {
    if ((theGpo->ElemColored == YES) ||
        (theGpo->PlotRefMarks == YES))
    {
      UserWrite("use i option only without c and r option\n");
      return (NOT_ACTIVE);
    }
  }

  if (theGpo->FreeBnd!=NULL)
  {
    if (VD_ncmps_in_otype_mod(theGpo->FreeBnd,NODEVEC,NON_STRICT)!=DIM)
      return (NOT_ACTIVE);

    if (!VD_SUCC_COMP(theGpo->FreeBnd))
      return (NOT_ACTIVE);
  }

  return (ACTIVE);
}

/****************************************************************************/
/*
   DisplayGridPlotObject_2D - Display content of 2D grid object

   SYNOPSIS:
   static INT DisplayGridPlotObject_2D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 2D grid object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayGridPlotObject_2D (PLOTOBJ *thePlotObj)
{
  struct GridPlotObj2D *theGpo;

  theGpo = &(thePlotObj->theGpo);

  /* print content */
  UserWriteF(DISPLAY_PO_FORMAT_SF,"ShrinkFactor",(float)theGpo->ShrinkFactor);
        #ifdef ModelP
  UserWriteF(DISPLAY_PO_FORMAT_SF,"PartShrinkFactor",(float)theGpo->PartShrinkFactor);
        #endif
  if (theGpo->PlotBoundary == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"BND","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"BND","NO");

  UserWriteF(DISPLAY_PO_FORMAT_SS,"Node markers",         BOOL_2_YN(theGpo->PlotNodes));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"ref marks",            BOOL_2_YN(theGpo->PlotRefMarks));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"indicator marks",      BOOL_2_YN(theGpo->PlotIndMarks));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"ElemID",                       BOOL_2_YN(theGpo->PlotElemID));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"subdomID",                     BOOL_2_YN(theGpo->PlotSubdomain));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"NodeID",                       BOOL_2_YN(theGpo->PlotNodeID));

  switch (theGpo->WhichElem)
  {
  case PO_COPY :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"WHICH_Elem","COPY");
    break;
  case PO_IRR :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"WHICH_Elem","IRREGULAR");
    break;
  case PO_REG :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"WHICH_Elem","REGULAR");
    break;
  case PO_ALL :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"WHICH_Elem","ALL");
    break;
  }

  UserWriteF(DISPLAY_PO_FORMAT_SI,"COLORED",(int)theGpo->ElemColored);

  if (theGpo->FreeBnd!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"free bnd",ENVITEM_NAME(theGpo->FreeBnd));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"free bnd","NO");

  return (0);
}

/****************************************************************************/
/*
   InitVecMat_2D - Initialization of 2D vector-matrix graph object

   SYNOPSIS:
   static INT InitVecMat_2D (PLOTOBJ *thePlotObj, INT argc, char **argv);

   PARAMETERS:
   .  thePlotObj - the 'PLOTOBJ' to be initialized.
   .  argc, argv -

   DESCRIPTION:
   This function does initialization of 2D vector-matrix graph object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   VecMat - plot object

   please lookup 'VecMat2D' or 'VecMat3D' depending on the space dimension you are
   working with (the options are not quite the same)
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   VecMat2D - plot object for vector and matrix data (stencils) (omit 2D when using it)

   DESCRIPTION:
   The VecMat plot object shows vector and matrix data using stencils.

   Possible options:~
   .   $t~[nd][ed][el]		- display only specified vector types and their connections
   .   $o~0|1|2|3               - show VECTOR/BLOCKVECTOR order (see also below)
   .   $m~0|1                   - plot markers for vectors, meaning depending on order mode, see below
   .   $i~0|1                   - plot indices for vectors, meaning depending on order mode, see below
   .n							if order = 0 (default)
   .n							  marker:
   .n							    circle:    node-vector
   .n							    rhombus:   edge-vector
   .n							    square:    elem-vector
   .n							  color:
   .n							    magenta:   VCLASS 0
   .n							    black:     VCLASS 1
   .n							    yellow:    VCLASS 2
   .n							    red:       VCLASS 3
   .n							  index: the position of the vector in the vector list

   .n							if order = 1
   .n							  marker:
   .n							    filled circles
   .n							  color: indicates the BLOCKVECTOR a VECTOR belongs to
   .n							  index: indicates BLOCKVECTOR number

   .n							if order = 2 (fitting well for orderv or lineorderv)
   .n							  marker:
   .n							    circle:    FIRST-set
   .n							    rhombus:   LAST-set
   .n							    square:    CUT-set
   .n							  color: indicates the cycle in which the vector was found
   .n							  index: indicates gen (FCL) and cycle number (e.g. F1)

   .n							if order = 3 (fitting well for lineorderv)
   .n							  marker:
   .n							    circle:    FIRST-set
   .n							    rhombus:   LAST-set
   .n							    square:    CUT-set
   .n							  color: indicates the line in which the vector was found
   .n							  index: indicates gen (FCL), cycle number, line number and
   .n							         number in line (e.g. F_1,12^3)

   .   $d~0|1                   - show dependencies indicated by arrows
   .   $c~0|1                   - also show connections (black)
   .   $e~0|1                   - also show extra connections (cyan)
   .   $C~0|1                   - visualize doubly linked vector list ($e and $c are reset to 0)
                                                        vectors are connected by lines in the order they appear in the list
   .   $b~0|1                   - also plot boundary
   .   $M~[<matdata~desc>] - plot user data of <matdata desc> for the selected vectors
   .n							and their neighbours. If <matdata desc> is omitted no matrix
   .n							data will be plotted
   .   $V~[<vecdata~desc>] - plot user data of <vecdata desc> for the selected vectors.
   .n							If <vecdata desc> is omitted no vector data will be plotted

   KEYWORDS:
   graphics, plot, window, picture, plotobject, vecdesc, matdesc, stencils, order
   D*/
/****************************************************************************/

static INT InitVecMat_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  FORMAT *theFormat;
  struct VecMatPlotObj2D *theVmo;
  char name[NAMELEN];
  INT i,j,rt,ct;
  int iValue;

  theVmo = &(thePlotObj->theVmo);
  theFormat = MGFORMAT(PO_MG(thePlotObj));
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theVmo->Marker                          = NO;
    for (i=0; i<MAXVECTORS; i++)
      theVmo->Type[i]                 = (FMT_S_VEC_TP(theFormat,i)>0);
    theVmo->Connections                     = YES;
    theVmo->Extra                           = NO;
    theVmo->Idx                                     = NO;
    theVmo->Part                            = NO;
    theVmo->Order                           = 0;
    theVmo->Dependency                      = NO;
    theVmo->ConnectVectors          = NO;
    theVmo->Boundary                        = YES;
    theVmo->vd                                      = NULL;
    theVmo->md                                      = NULL;
  }

  /* color mode */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'm' :
      if (sscanf(argv[i],"m %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Marker = YES;
      else if (iValue==0) theVmo->Marker = NO;
      break;

    case 't' :
      for (j=0; j<NVECTYPES; j++)
        if (FMT_S_VEC_TP(theFormat,j)>0)
          if (strchr(argv[i]+1,FMT_VTYPE_NAME(theFormat,j))!=NULL)
            theVmo->Type[j] = YES;
          else
            theVmo->Type[j] = NO;
      break;

    case 'c' :
      if (sscanf(argv[i],"c %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Connections = YES;
      else if (iValue==0) theVmo->Connections = NO;
      break;

    case 'C' :
      if (sscanf(argv[i],"C %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->ConnectVectors = YES;
      else if (iValue==0) theVmo->ConnectVectors = NO;
      break;

    case 'e' :
      if (sscanf(argv[i],"e %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Extra = YES;
      else if (iValue==0) theVmo->Extra = NO;
      break;

    case 'i' :
      if (sscanf(argv[i],"i %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Idx = YES;
      else if (iValue==0) theVmo->Idx = NO;
      break;

    case 'p' :
      if (sscanf(argv[i],"p %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Part = YES;
      else if (iValue==0) theVmo->Part = NO;
      break;

    case 'd' :
      if (sscanf(argv[i],"d %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Dependency = YES;
      else if (iValue==0) theVmo->Dependency = NO;
      break;

    case 'o' :
      if (sscanf(argv[i],"o %d",&iValue)!=1)
        break;
      theVmo->Order = iValue;
      theVmo->Order = MAX(theVmo->Order,0);
      theVmo->Order = MIN(theVmo->Order,3);
      break;

    case 'b' :
      if (sscanf(argv[i],"b %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Boundary = YES;
      else if (iValue==0) theVmo->Boundary = NO;
      break;

    case 'V' :
      if ((sscanf(argv[i],"V %s",name)!=1) ||
          ((theVmo->vd =
              GetVecDataDescByName(PO_MG(thePlotObj),name))==NULL))
      {
        UserWrite("no vector specified, vec data switched off\n");
        theVmo->vd = NULL;
      }
      break;

    case 'M' :
      if ((sscanf(argv[i],"M %s",name)!=1) ||
          ((theVmo->md =
              GetMatDataDescByName(PO_MG(thePlotObj),name))==NULL))
      {
        UserWrite("no matrix specified, mat data switched off\n");
        theVmo->md = NULL;
      }
      break;
    }

  if (theVmo->ConnectVectors)
  {
    theVmo->Connections = NO;
    theVmo->Extra           = NO;
  }

  if (theVmo->vd && theVmo->md)
    /* check compatibility of vec and mat desc */
    for (rt=0; rt<NVECTYPES; rt++)
      if (theVmo->Type[rt])
        for (ct=0; ct<NVECTYPES; ct++)
          if (theVmo->Type[ct])
            if (VD_NCMPS_IN_TYPE(theVmo->vd,ct)!=MD_ROWS_IN_RT_CT(theVmo->md,rt,ct))
            {
              UserWrite("vec desc and mat desc incompatible\n");
              return (NOT_ACTIVE);
            }
  return (ACTIVE);
}

/****************************************************************************/
/*
   DisplayVecMat_2D - Display content of 2D vector-matrix graph object

   SYNOPSIS:
   static INT DisplayVecMat_2D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 2D vector-matrix graph object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayVecMat_2D (PLOTOBJ *thePlotObj)
{
  FORMAT *fmt;
  struct VecMatPlotObj2D *theVmo;
  char buffer[16];
  INT i;

  theVmo = &(thePlotObj->theVmo);

  /* print content */
  if (theVmo->Marker == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"marker","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"marker","NO");

  fmt = MGFORMAT(PO_MG(thePlotObj));
  for (i=0; i<MAXVECTORS; i++)
    if (FMT_S_VEC_TP(fmt,i)>0)
    {
      sprintf(buffer,"type %c",FMT_VTYPE_NAME(fmt,i));
      UserWriteF(DISPLAY_PO_FORMAT_SS,buffer,BOOL_2_YN(theVmo->Type[i]));
    }

  UserWriteF(DISPLAY_PO_FORMAT_SS,"index",                BOOL_2_YN(theVmo->Idx));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"connections",  BOOL_2_YN(theVmo->Connections));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"extra",                BOOL_2_YN(theVmo->Extra));

  switch (theVmo->Order)
  {
  case 1 : UserWriteF(DISPLAY_PO_FORMAT_SS,"order","1: blockvector order"); break;
  case 2 : UserWriteF(DISPLAY_PO_FORMAT_SS,"order","2: vector order, cuts black"); break;
  case 3 : UserWriteF(DISPLAY_PO_FORMAT_SS,"order","3: line order"); break;
  }

  UserWriteF(DISPLAY_PO_FORMAT_SS,"dependency",   BOOL_2_YN(theVmo->Dependency));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"connect",              BOOL_2_YN(theVmo->ConnectVectors));
  UserWriteF(DISPLAY_PO_FORMAT_SS,"boundary",             BOOL_2_YN(theVmo->Boundary));

  if (theVmo->vd!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"vec data",ENVITEM_NAME(theVmo->vd));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"vec data","NO");
  if (theVmo->md!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"mat data",ENVITEM_NAME(theVmo->md));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"mat data","NO");

  return (0);
}

/****************************************************************************/
/*
   InitScalarFieldPlotObject_2D	- Initialization of 2D scalar field object

   SYNOPSIS:
   static INT InitScalarFieldPlotObject_2D (PLOTOBJ *thePlotObj, INT argc,
   char **argv);

   PARAMETERS:
   .  thePlotObj -
   .  argc -
   .  argv -

   DESCRIPTION:
   This function makes an initialization of 2D scalar field object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   EScalar - plot object

   please lookup 'EScalar2D' or 'EScalar3D' depending on the space dimension you are
   working with (the options are not quite the same)
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   EScalar2D - plot object for scalar grid functions (omit 2D when using it)

   DESCRIPTION:
   The EScalar plot object shows a contour or color plot of a scalar grid function.

   Mandatory options are:~
   .    $e~<vec~eval~proc>		- either specify <vec eval proc>
   .    $s~<vecdata~desc>		- or			 <vecdata desc> (standard eval proc, nodal values only)

   Possible options:~
   .    $f~<from>				- from value
   .    $t~<to>				- to value
   .    $d~<depth>				- integer for recursive depth (default 0,
                                                          CAUTION: may slow down dramatically!)
   .    $m COLOR|CONTOURS_EQ	- mode: color or equidistant contour lines
   .    $n~<levels>			- number of levels for contour lines
   .    $g~0|1					- show grid off/on

   KEYWORDS:
   graphics, plot, window, picture, plotobject, vecdesc, function
   D*/
/****************************************************************************/

static INT InitScalarFieldPlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  struct ElemScalarPlotObj2D *theEspo;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEspo = &(thePlotObj->theEspo);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  ret = ACTIVE;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theEspo->min                    = 0.0;
    theEspo->max                    = 1.0;
    theEspo->mode                   = PO_COLOR;
    theEspo->PlotGrid               = NO;
    theEspo->depth                  = 0;
    theEspo->numOfContours  = 10;
  }

  /* set plot grid option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='g')
    {
      if (sscanf(argv[i],"g %d",&iValue)!=1)
        break;
      if (iValue==1)
        theEspo->PlotGrid = YES;
      else if (iValue==0)
        theEspo->PlotGrid = NO;
      break;
    }

  /* set mode option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='m')
    {
      if (sscanf(argv[i],"m %s",buffer)!=1)
        break;
      if (strcmp(buffer,"COLOR") == 0)
        theEspo->mode = PO_COLOR;
      else if (strcmp(buffer,"CONTOURS_EQ") == 0)
        theEspo->mode = PO_CONTOURS_EQ;
      break;
    }

  /* set from option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='f')
    {
      if (sscanf(argv[i],"f %g",&fValue)!=1)
        break;
      theEspo->min = fValue;
      break;
    }

  /* set to option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='t')
    {
      if (sscanf(argv[i],"t %g",&fValue)!=1)
        break;
      theEspo->max = fValue;
      break;
    }
  if (theEspo->min >= theEspo->max )
  {
    UserWrite("minValue is bigger than maxValue\n");
    ret = NOT_ACTIVE;
  }

  /* set depth option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='d')
    {
      if (sscanf(argv[i],"d %d",&iValue)!=1)
        break;
      theEspo->depth = iValue;
      break;
    }
  if (theEspo->depth<0 || theEspo->depth>4)
  {
    UserWrite("depth is not valid\n");
    ret = NOT_ACTIVE;
  }

  /* set n option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='n')
    {
      if (sscanf(argv[i],"n %d",&iValue)!=1)
        break;
      if (iValue>1)
        theEspo->numOfContours = iValue;
      break;
    }
  if (theEspo->numOfContours <= 1 )
  {
    UserWrite("number of contours is smaller than 1\n");
    ret = NOT_ACTIVE;
  }

  if (theEspo->numOfContours >= PO_MAXCONTOURS)
  {
    PrintErrorMessageF('E',"InitScalarFieldPlotObject_2D","number of contours is greater than the limit (%d)",PO_MAXCONTOURS);
    ret = NOT_ACTIVE;
  }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEspo),buffer);
      theEspo->EvalFct = GetElementValueEvalProc(buffer);
      break;
    }

  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEspo),buffer);
      if (theEspo->EvalFct == NULL)
        theEspo->EvalFct = GetElementValueEvalProc("nvalue");
      break;
    }

  if (theEspo->EvalFct == NULL)
  {
    UserWrite("cannot find plot procedure\n");
    ret = NOT_ACTIVE;
  }

  /* do what is to do */
  if (theEspo->mode == PO_CONTOURS_EQ && ret == ACTIVE)
    for (i=0; i<theEspo->numOfContours; i++)
      theEspo->contValues[i] = theEspo->min + (DOUBLE)i * (theEspo->max - theEspo->min) / (DOUBLE)(theEspo->numOfContours-1);

  /* return */
  return (ret);
}

/****************************************************************************/
/*
   DisplayScalarFieldPlotObject_2D - Display content of 2D scalar field object

   SYNOPSIS:
   static INT DisplayScalarFieldPlotObject_2D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 2D scalar field object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayScalarFieldPlotObject_2D (PLOTOBJ *thePlotObj)
{
  struct ElemScalarPlotObj2D *theEspo;

  theEspo = &(thePlotObj->theEspo);

  /* print content */
  if (theEspo->EvalFct != NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theEspo->EvalFct));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc","---");

  UserWriteF(DISPLAY_PO_FORMAT_SS,"name",PO_NAME(theEspo));

  if (theEspo->PlotGrid == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"Grid","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"Grid","NO");

  UserWriteF(DISPLAY_PO_FORMAT_SFF,"Range",(float)theEspo->min,(float)theEspo->max);
  UserWriteF(DISPLAY_PO_FORMAT_SI,"Depth",(int)theEspo->depth);
  switch (theEspo->mode)
  {
  case PO_COLOR :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"PlotMode","COLOR");
    break;
  case PO_CONTOURS_EQ :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"PlotMode","CONTOURS_EQ");
    UserWriteF(DISPLAY_PO_FORMAT_SI,"NbOfCont",(int)theEspo->numOfContours);
  }

  return (0);
}

/****************************************************************************/
/*
   InitVectorFieldPlotObject_2D - Initialization of 2D vector field object

   SYNOPSIS:
   static INT InitVectorFieldPlotObject_2D (PLOTOBJ *thePlotObj, INT argc,
   char **argv);

   PARAMETERS:
   .  thePlotObj -
   .  argc -
   .  argv -

   DESCRIPTION:
   This function initializes 2D vector field object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   EVector - plot object

   please lookup 'EVector2D' or 'EVector3D' depending on the space dimension you are
   working with (the options are not quite the same)
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   EVector2D - plot object for (two-)vector grid functions (omit 2D when using it)

   DESCRIPTION:
   The EVector plot object shows a vector plot of a grid function. The vectors
   are plotted in a regular raster.

   Mandatory options are:~
   .    $e~<vec~eval~proc>		- either specify <vec eval proc>
   .    $s~<vecdata~desc>		- or			 <vecdata desc> (standard eval proc, nodal values only)

   Possible options:~
   .    $t~<to>				- to value
   .    $r~<raster>			- raster size
   .    $l~<cut~len>			- cut off len relative to raster size (default 1 to avoid overlap)
   .    $c~0|1					- cut off vectors off/on
   .    $g~0|1					- show grid off/on

   KEYWORDS:
   graphics, plot, window, picture, plotobject, vecdesc, function
   D*/
/****************************************************************************/

static INT InitVectorFieldPlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  PICTURE *pic;
  struct ElemVectorPlotObj2D *theEvpo;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEvpo = &(thePlotObj->theEvpo);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  pic = PO_PIC(thePlotObj);
  ret = ACTIVE;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theEvpo->PlotGrid               = NO;
    theEvpo->max                    = 1.0;
    theEvpo->CutVectors     = YES;
    theEvpo->RasterSize     = PO_RADIUS(thePlotObj)/10.0;
    theEvpo->CutLenFactor   = 1.0;
  }

  /* set plot grid option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='g')
    {
      if (sscanf(argv[i],"g %d",&iValue)!=1)
        break;
      if (iValue==1)
        theEvpo->PlotGrid = YES;
      else if (iValue==0)
        theEvpo->PlotGrid = NO;
      break;
    }

  /* set to option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='t')
    {
      if (sscanf(argv[i],"t %g",&fValue)!=1)
        break;
      theEvpo->max = fValue;
      break;
    }
  if (theEvpo->max <= 0.0)
  {
    UserWrite("maxValue is smaller than zero\n");
    ret =NOT_ACTIVE;
  }

  /* set rastersize option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='r')
    {
      if (sscanf(argv[i],"r %g",&fValue)!=1)
        break;
      if (fValue<3.)
      {
        /* TODO (HRR 971010): replace this by {PrintErrorMessage(); ret =NOT_ACTIVE;} once all scripts are changed */
        printf("ERROR: die Rasterweite von EVector in 2D muss in --> PIXELN <-- angegeben werden\n");
        ASSERT(FALSE);
      }
      if (fValue>MIN(fabs(PIC_GLL(pic)[_X_]-PIC_GUR(pic)[_X_]),fabs(PIC_GLL(pic)[_Y_]-PIC_GUR(pic)[_Y_]))/2)
      {
        PrintErrorMessage('E',"InitVectorFieldPlotObject_2D","rastersize > half picture size");
        ret =NOT_ACTIVE;
      }
      theEvpo->RasterSize = fValue;
      break;
    }
  if (theEvpo->RasterSize <= 0.0)
  {
    UserWrite("RasterSize is smaller than zero\n");
    ret =NOT_ACTIVE;
  }

  /* set CutLenFactor option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='l')
    {
      if (sscanf(argv[i],"l %g",&fValue)!=1)
        break;
      theEvpo->CutLenFactor = fValue;
      break;
    }
  if (theEvpo->CutLenFactor < 0.1 || theEvpo->CutLenFactor > 10)
  {
    UserWrite("CutLenFactor is not in [0.1,10]\n");
    ret =NOT_ACTIVE;
  }

  /* set cut vector option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='c')
    {
      if (sscanf(argv[i],"c %d",&iValue)!=1)
        break;
      if (iValue==1)
        theEvpo->CutVectors = YES;
      else if (iValue==0)
        theEvpo->CutVectors = NO;
      break;
    }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEvpo),buffer);
      theEvpo->EvalFct = GetElementVectorEvalProc(buffer);
      break;
    }

  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEvpo),buffer);
      if (theEvpo->EvalFct == NULL)
        theEvpo->EvalFct = GetElementVectorEvalProc("nvector");
      break;
    }

  if (theEvpo->EvalFct == NULL)
  {
    UserWrite("cannot find plot procedure\n");
    ret = NOT_ACTIVE;
  }

  /* return */
  return (ret);
}

/****************************************************************************/
/*
   DisplayVectorFieldPlotObject_2D - Display content of 2D vector field object

   SYNOPSIS:
   static INT DisplayVectorFieldPlotObject_2D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 2D vector field object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayVectorFieldPlotObject_2D (PLOTOBJ *thePlotObj)
{
  struct ElemVectorPlotObj2D *theEvpo;

  theEvpo = &(thePlotObj->theEvpo);

  UserWriteF(DISPLAY_PO_FORMAT_SS,"name",PO_NAME(theEvpo));

  /* print content */
  if (theEvpo->EvalFct != NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theEvpo->EvalFct));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc","---");
  UserWriteF(DISPLAY_PO_FORMAT_SS,"Grid",(theEvpo->PlotGrid) ? "YES" : "NO");
  UserWriteF(DISPLAY_PO_FORMAT_SF,"maxValue",(float)theEvpo->max);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"RasterSize",(float)theEvpo->RasterSize);
  if (theEvpo->CutVectors == YES)
  {
    UserWriteF(DISPLAY_PO_FORMAT_SS,"CutVectors","YES");
    UserWriteF(DISPLAY_PO_FORMAT_SF,"CutLenFactor",(float)theEvpo->CutLenFactor);
  }
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"CutVectors","NO");

  return (0);
}

/****************************************************************************/
/*
   InitLinePlotObject_2D	- Initialization of 2D line object

   SYNOPSIS:
   static INT InitLinePlotObject_2D (PLOTOBJ *thePlotObj, INT argc,
   char **argv);

   PARAMETERS:
   .  thePlotObj -
   .  argc -
   .  argv -

   DESCRIPTION:
   This function makes an initialization of 2D line object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   Line - plot object to show a section through a 2D scalar grid function

   DESCRIPTION:
   The Line plot object shows a section through a 2D scalar grid function.

   Mandatory options are:~
   .    $e~<vec~eval~proc>		- either specify <vec eval proc>
   .    $s~<vecdata~desc>		- or			 <vecdata desc> (standard eval proc, nodal values only)

   Possible options:~
   .    $f~<from>				- from value
   .    $t~<to>				- to value
   .    $d~<depth>				- integer for recursive depth (default 0,
                                                          CAUTION: may slow down dramatically!)
   .    $a~<ratio>				- use anisotropic scaling of x- and y-axis
   .    $c~<value>				- use this color
   .    $l~<x>~<y>				- starting point
   .    $r~<x>~<y>				- end point
   .    $Ly~0|1				- use logarithmic scaale for values off/on

   KEYWORDS:
   graphics, plot, window, picture, plotobject, vecdesc, section
   D*/
/****************************************************************************/

static INT InitLinePlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  struct LinePlotObj2D *theLpo;
  INT i, ret;
  int iValue;
  float fValue[2];
  DOUBLE dist;
  char buffer[128];

  theLpo = &(thePlotObj->theLpo);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  PO_MIDPOINT(thePlotObj)[0] = PO_MIDPOINT(thePlotObj)[1] = 0.5;
  PO_RADIUS(thePlotObj) = 0.70711;
  ret = ACTIVE;

  theLpo->nHit = 0;
  theLpo->xmin = 1.0;
  theLpo->xmax = 0.0;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theLpo->min                                                             = 0.0;
    theLpo->max                                                             = 1.0;
    theLpo->yLog                                                    = 0;
    theLpo->left[0] = theLpo->left[1]               = 0.0;
    theLpo->right[0] = theLpo->right[1]     = 1.0;
    theLpo->color                                                   = 0.0;
    theLpo->aspectratio                                             = 1.0;
    theLpo->EvalFct                                                 = NULL;
  }

  /* set from option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='f')
    {
      if (sscanf(argv[i],"f %g",fValue)!=1)
        break;
      theLpo->min = fValue[0];
      break;
    }

  /* set to option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='t')
    {
      if (sscanf(argv[i],"t %g",fValue)!=1)
        break;
      theLpo->max = fValue[0];
      break;
    }
  if (theLpo->min >= theLpo->max )
  {
    UserWrite("minValue is bigger than maxValue\n");
    ret = NOT_ACTIVE;
  }

  /* set left option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='l')
    {
      if (sscanf(argv[i],"l %g %g",fValue,fValue+1)!=2)
        break;
      V2_COPY(fValue,theLpo->left)
      break;
    }

  /* set right option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='r')
    {
      if (sscanf(argv[i],"r %g %g",fValue,fValue+1)!=2)
        break;
      V2_COPY(fValue,theLpo->right)
      break;
    }
  V2_EUKLIDNORM_OF_DIFF(theLpo->left,theLpo->right,dist)
  if (dist==0.0)
  {
    UserWrite("left and right have to be different\n");
    ret = NOT_ACTIVE;
  }

  /* set color option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='c')
    {
      if (sscanf(argv[i],"c %g",fValue)!=1)
        break;
      theLpo->color = fValue[0];
      break;
    }
  if (theLpo->color<0.0 || theLpo->color>1.0)
  {
    UserWrite("color is not valid\n");
    ret = NOT_ACTIVE;
  }

  /* set aspectratio option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='a')
    {
      if (sscanf(argv[i],"a %g",fValue)!=1)
        break;
      theLpo->aspectratio = fValue[0];
      break;
    }
  if (theLpo->aspectratio<=0.0)
  {
    UserWrite("aspect ratio is not valid\n");
    ret = NOT_ACTIVE;
  }

  /* log-y option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='L')
      if (sscanf(argv[i],"Ly %d",&iValue)==1)
      {
        theLpo->yLog = iValue;
        break;
      }

  /* set depth option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='d')
    {
      if (sscanf(argv[i],"d %d",&iValue)!=1)
        break;
      theLpo->depth = iValue;
      break;
    }
  if (theLpo->depth<0 || theLpo->depth>4)
  {
    UserWrite("depth is not valid\n");
    ret = NOT_ACTIVE;
  }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theLpo),buffer);
      theLpo->EvalFct = GetElementValueEvalProc(buffer);
      break;
    }
  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theLpo),buffer);
      if (theLpo->EvalFct == NULL)
        theLpo->EvalFct = GetElementValueEvalProc("nvalue");
      break;
    }
  if (theLpo->EvalFct == NULL)
  {
    UserWrite("cannot find plot procedure\n");
    ret = NOT_ACTIVE;
  }

  /* midpoint and radius */
  PO_MIDPOINT(thePlotObj)[0] = 0.5; PO_MIDPOINT(thePlotObj)[1] = 0.5*theLpo->aspectratio;
  PO_RADIUS(thePlotObj) = 0.5 * SQRT(1.0 + theLpo->aspectratio*theLpo->aspectratio);

  return (ret);
}

/****************************************************************************/
/*
   DisplayLinePlotObject_2D - Display content of 2D line object

   SYNOPSIS:
   static INT DisplayLinePlotObject_2D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 2D line object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayLinePlotObject_2D (PLOTOBJ *thePlotObj)
{
  struct LinePlotObj2D *theLpo;

  theLpo = &(thePlotObj->theLpo);

  /* print content */
  if (theLpo->EvalFct != NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theLpo->EvalFct));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc","---");

  UserWriteF(DISPLAY_PO_FORMAT_SS,"name",PO_NAME(theLpo));

  UserWriteF(DISPLAY_PO_FORMAT_SFF,"Range",(float)theLpo->min,(float)theLpo->max);
  UserWriteF(DISPLAY_PO_FORMAT_SFF,"left",(float)theLpo->left[0],(float)theLpo->left[1]);
  UserWriteF(DISPLAY_PO_FORMAT_SFF,"right",(float)theLpo->right[0],(float)theLpo->right[1]);
  UserWriteF(DISPLAY_PO_FORMAT_SI,"y-log",(int)theLpo->yLog);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"color",(float)theLpo->color);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"asp.ratio",(float)theLpo->aspectratio);
  UserWriteF(DISPLAY_PO_FORMAT_SI,"Depth",(int)theLpo->depth);
  UserWrite("\ncomputed values:\n");
  UserWriteF(DISPLAY_PO_FORMAT_SI,"nHit",(int)theLpo->nHit);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"x-min",(float)theLpo->xmin);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"x-max",(float)theLpo->xmax);

  return (0);
}

#endif

#ifdef __THREEDIM__

/****************************************************************************/
/*																			*/
/* Function:  InitDomainPlotObject_3D										*/
/*																			*/
/* Purpose:   initialization of 3D domain object							*/
/*																			*/
/* Input:	  PLOTOBJ *thePlotObj, INT argc, char **argv					*/
/*																			*/
/* Output:	  INT NOT_INIT, NOT_ACTIVE or ACTIVE							*/
/*																			*/
/****************************************************************************/

static INT InitDomainPlotObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  struct DomainPlotObj3D *theDpo;
  INT i;
  int iValue;

  theDpo = &(thePlotObj->theDpo);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V3_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theDpo->NbOfSteps               = 1;
  }

  /* set n option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='n')
    {
      if (sscanf(argv[i],"n %d",&iValue)!=1)
        break;
      if (iValue>1)
        theDpo->NbOfSteps = iValue;
      break;
    }
  if (theDpo->NbOfSteps <= 0 || theDpo->NbOfSteps>1000)
  {
    UserWrite("number of steps must ly in (1,1000)\n");
    return (NOT_ACTIVE);
  }

  return (ACTIVE);
}

/****************************************************************************/
/*																			*/
/* Function:  DisplayDomainPlotObject_3D									*/
/*																			*/
/* Purpose:   display content of 3D domain field object                                         */
/*																			*/
/* Input:	  PLOTOBJ *thePlotObj											*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT DisplayDomainPlotObject_3D (PLOTOBJ *thePlotObj)
{
  struct DomainPlotObj3D *theDpo;

  theDpo = &(thePlotObj->theDpo);

  /* print content */
  UserWriteF("NbOfSteps  =%d\n",(int)theDpo->NbOfSteps);

  return (0);
}

/****************************************************************************/
/*
   InitVecMat_3D - Initialization of 2D vector-matrix graph object

   SYNOPSIS:
   static INT InitVecMat_3D (PLOTOBJ *thePlotObj, INT argc, char **argv);

   PARAMETERS:
   .  thePlotObj - the 'PLOTOBJ' to be initialized.
   .  argc, argv -

   DESCRIPTION:
   This function does initialization of 2D vector-matrix graph object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   VecMat3D - plot object for vector and matrix data (stencils) (omit 3D when using it)

   DESCRIPTION:
   The VecMat plot object shows vector and matrix data using stencils.

   Possible options:~
   .   $t~[nd][ed][el][sd]	- display only specified vector types and their connections
   .   $i~0|1                   - plot indices for vectors
   .   $M~[<matdata~desc>] - plot user data of <matdata desc> for the selected vectors
   .n							and their neighbours. If <matdata desc> is omitted no matrix
   .n							data will be plotted
   .   $V~[<vecdata~desc>] - plot user data of <vecdata desc> for the selected vectors.
   .n							If <vecdata desc> is omitted no vector data will be plotted

   KEYWORDS:
   graphics, plot, window, picture, plotobject, vecdesc, matdesc, stencils, order
   D*/
/****************************************************************************/

static INT InitVecMat_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  FORMAT *theFormat;
  struct VecMatPlotObj3D *theVmo;
  char name[NAMELEN];
  INT i,j,rt,ct;
  int iValue;

  theVmo = &(thePlotObj->theVmo);
  theFormat = MGFORMAT(PO_MG(thePlotObj));
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    for (i=0; i<MAXVECTORS; i++)
      theVmo->Type[i]                 = (FMT_S_VEC_TP(theFormat,i)>0);
    theVmo->vd                                      = NULL;
    theVmo->md                                      = NULL;
  }

  /* color mode */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 't' :
      for (j=0; j<NVECTYPES; j++)
        if (FMT_S_VEC_TP(theFormat,j)>0)
          if (strchr(argv[i]+1,FMT_VTYPE_NAME(theFormat,j))!=NULL)
            theVmo->Type[j] = YES;
          else
            theVmo->Type[j] = NO;
      break;

    case 'i' :
      if (sscanf(argv[i],"i %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Idx = YES;
      else if (iValue==0) theVmo->Idx = NO;
      break;

    case 'V' :
      if ((sscanf(argv[i],"V %s",name)!=1) ||
          ((theVmo->vd =
              GetVecDataDescByName(PO_MG(thePlotObj),name))==NULL))
      {
        UserWrite("no vector specified, vec data switched off\n");
        theVmo->vd = NULL;
      }
      break;

    case 'M' :
      if ((sscanf(argv[i],"M %s",name)!=1) ||
          ((theVmo->md =
              GetMatDataDescByName(PO_MG(thePlotObj),name))==NULL))
      {
        UserWrite("no matrix specified, mat data switched off\n");
        theVmo->md = NULL;
      }
      break;
    }

  /* check compatibility of vec and mat desc */
  if (theVmo->vd || theVmo->md)
    for (rt=0; rt<NVECTYPES; rt++)
      if (theVmo->Type[rt])
      {
        if (theVmo->vd)
          if (!VD_ISDEF_IN_TYPE(theVmo->vd,rt))
          {
            UserWrite("vec desc does not include types of specified types\n");
            return (NOT_ACTIVE);
          }
        if (theVmo->md)
          for (ct=0; ct<NVECTYPES; ct++)
            if (theVmo->Type[ct])
            {
              if (!MD_ISDEF_IN_RT_CT(theVmo->md,rt,ct))
              {
                UserWrite("mat desc does not include column types of specified types\n");
                return (NOT_ACTIVE);
              }
              if (theVmo->vd)
                if (VD_NCMPS_IN_TYPE(theVmo->vd,ct)!=MD_ROWS_IN_RT_CT(theVmo->md,rt,ct))
                {
                  UserWrite("vec desc and mat desc incompatible\n");
                  return (NOT_ACTIVE);
                }
            }
      }

  return (ACTIVE);
}

/****************************************************************************/
/*
   DisplayVecMat_3D - Display content of 2D vector-matrix graph object

   SYNOPSIS:
   static INT DisplayVecMat_3D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 2D vector-matrix graph object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayVecMat_3D (PLOTOBJ *thePlotObj)
{
  FORMAT *fmt;
  struct VecMatPlotObj3D *theVmo;
  INT i;
  char buffer[128];

  theVmo = &(thePlotObj->theVmo);

  /* print content */
  fmt = MGFORMAT(PO_MG(thePlotObj));
  for (i=0; i<MAXVECTORS; i++)
    if (FMT_S_VEC_TP(fmt,i)>0)
    {
      sprintf(buffer,"type %c",FMT_VTYPE_NAME(fmt,i));
      UserWriteF(DISPLAY_PO_FORMAT_SS,buffer,BOOL_2_YN(theVmo->Type[i]));
    }

  UserWriteF(DISPLAY_PO_FORMAT_SS,"index",(theVmo->Idx) ? "YES" : "NO");

#ifdef __NP__
  if (theVmo->vd!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"vec data",ENVITEM_NAME(theVmo->vd));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"vec data","NO");
  if (theVmo->md!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"mat data",ENVITEM_NAME(theVmo->md));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"mat data","NO");
#endif

  return (0);
}

/****************************************************************************/
/*
   InitGridObject_3D - Initialization of 3D grid object

   SYNOPSIS:
   static INT InitGridObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv);

   PARAMETERS:
   .  thePlotObj -
   .  argc -
   .  argv -

   DESCRIPTION:
   This function makes an initialization of 3D grid object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   Grid3D - plot object for multigrid (omit 3D when using it)

   DESCRIPTION:
   The Grid plot object shows the multigrid with its elements.

   Possible options:~
   .    $c~0|1					- use color off/on
   .    $w~c|i|r|a				- show copy/irregular/regular/all elements
   .    $n[i]~0|1				- plot node markers (and IDs)
   .    $v[i]~0|1				- plot vector markers (and indices)
   .    $t~<type~list>			- only vectors of specified types
   .    <type~list>			- a list composed by any of nd, ed, el, si, seperated by blanks
   .    $s~<shrink>			- factor to shrink elements
   .    $a~0..1                - contribution of ambient light to face intensity
   .    $p~<shrink>			- parallel only: factor to shrink processor partition

   KEYWORDS:
   graphics, plot, window, picture, plotobject, multigrid, elements
   D*/
/****************************************************************************/

static INT InitGridObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  FORMAT *theFormat;
  MULTIGRID *theMG;
  struct GridPlotObj3D *theGpo;
  INT i,j;
  int iValue;
  float fValue;
  char buffer[64],c;

  theGpo = &(thePlotObj->theGpo);
  theMG = PO_MG(thePlotObj);
  theFormat = MGFORMAT(theMG);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V3_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  PO_USESCUT(thePlotObj) = YES;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theGpo->ShrinkFactor            = 1.0;
                #ifdef ModelP
    theGpo->PartShrinkFactor        = 1.0;
                #endif
    theGpo->NodeMarkers                     = NO;
    theGpo->NodeIndex                       = NO;
    theGpo->Vectors                         = NO;
    theGpo->VecIndex                        = NO;
    for (i=0; i<MAXVOBJECTS; i++)
      theGpo->OType[i]                = VEC_DEF_IN_OBJ_OF_MG(theMG,i);
    theGpo->ElemColored             = NO;
    theGpo->WhichElem                       = PO_ALL;
    theGpo->PlotSelection           = 0;
    theGpo->AmbientLight        = 1.0;
  }

  /* set shrink option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %f",&fValue)!=1)
        break;
      theGpo->ShrinkFactor = fValue;
      break;
    }
  if (theGpo->ShrinkFactor<=0.0 || theGpo->ShrinkFactor>1.0)
    return (NOT_ACTIVE);

        #ifdef ModelP
  /* set shrink option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='p')
    {
      if (sscanf(argv[i],"p %f",&fValue)!=1)
        break;
      theGpo->PartShrinkFactor = fValue;
      break;
    }
  if (theGpo->PartShrinkFactor<=0.0 || theGpo->PartShrinkFactor>1.0)
    return (NOT_ACTIVE);
        #endif

  /* color mode */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='c')
    {
      if (sscanf(argv[i],"c %d",&iValue)!=1) break;
      theGpo->ElemColored = iValue;
      break;
    }
  if (theGpo->ElemColored<0 || theGpo->ElemColored>2) return (NOT_ACTIVE);

  /* selection option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='S')
    {
      theGpo->PlotSelection = 1;
      break;
    }

  /* node markers and indices */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='n')
    {
      if (sscanf(argv[i],"n%c %d",&c,&iValue)!=2) break;
      theGpo->NodeMarkers = iValue;
      if (argv[i][1]=='i')
        theGpo->NodeIndex = iValue;
      break;
    }

  /* vector options */
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'v' :
      if (sscanf(argv[i],"v%c %d",&c,&iValue)!=2)
        break;
      if              (iValue==1) {theGpo->Vectors = YES; theGpo->NodeMarkers = NO;}
      else if (iValue==0) theGpo->Vectors = NO;
      if (argv[i][1]=='i')
        theGpo->VecIndex = iValue;
      break;

    case 't' :
      for (j=0; j<MAXVOBJECTS; j++)
        if (strstr(argv[i]+1,ObjTypeName[j])!=NULL)
        {
          if (!VEC_DEF_IN_OBJ_OF_MG(theMG,i))
            PrintErrorMessageF('W',"InitGridObject_3D",
                               "no degrees of freedom in %s-vectors",ObjTypeName[j]);
          else
            theGpo->OType[j] = YES;
        }
        else
          theGpo->OType[j] = NO;
      break;
    }

  for (i=1; i<argc; i++)
  {
    if (argv[i][0]=='w')
    {
      sscanf(argv[i],"w %s",buffer);
      if (buffer[0] == 'c')
        theGpo->WhichElem = PO_COPY;
      else if (buffer[0] == 'i')
        theGpo->WhichElem = PO_IRR;
      else if (buffer[0] == 'r')
        theGpo->WhichElem = PO_REG;
      else if (buffer[0] == 'a')
        theGpo->WhichElem = PO_ALL;
      break;
    }
  }

  for (i=1; i<argc; i++)
    if (argv[i][0]=='a')
    {
      if (sscanf(argv[i], "a %f", &fValue) != 1) break;
      theGpo->AmbientLight = fValue;
      break;
    }
  if (theGpo->AmbientLight < 0.0 || theGpo->AmbientLight > 1.0)
    theGpo->AmbientLight = 1.0;

  return (ACTIVE);
}

/****************************************************************************/
/*
   DisplayGridPlotObject_3D - Display content of 3D grid field object

   SYNOPSIS:
   static INT DisplayGridPlotObject_3D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 3D grid field object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayGridPlotObject_3D (PLOTOBJ *thePlotObj)
{
  struct GridPlotObj3D *theGpo;
  INT i;
  char buffer[128];

  theGpo = &(thePlotObj->theGpo);
  PO_USESCUT(thePlotObj) = YES;

  /* print content */
  UserWriteF(DISPLAY_PO_FORMAT_SF,"ShrinkFactor",(float)theGpo->ShrinkFactor);
        #ifdef ModelP
  UserWriteF(DISPLAY_PO_FORMAT_SF,"PartShrinkFactor",(float)theGpo->PartShrinkFactor);
        #endif
  UserWriteF(DISPLAY_PO_FORMAT_SI,"colered elems",(int)theGpo->ElemColored);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"AmbientLight", (float)theGpo->AmbientLight);

  switch (theGpo->WhichElem)
  {
  case PO_COPY :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"WHICH_Elem","COPY");
    break;
  case PO_IRR :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"WHICH_Elem","IRREGULAR");
    break;
  case PO_REG :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"WHICH_Elem","REGULAR");
    break;
  case PO_ALL :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"WHICH_Elem","ALL");
    break;
  }
  UserWriteF(DISPLAY_PO_FORMAT_SI,"node markers",(int)theGpo->NodeMarkers);
  UserWriteF(DISPLAY_PO_FORMAT_SI,"node indices",(int)theGpo->NodeIndex);
  UserWriteF(DISPLAY_PO_FORMAT_SI,"vector markers",(int)theGpo->Vectors);
  UserWriteF(DISPLAY_PO_FORMAT_SI,"vector indices",(int)theGpo->VecIndex);
  for (i=0; i<MAXVOBJECTS; i++)
  {
    sprintf(buffer,"vobject %s",ObjTypeName[i]);
    UserWriteF(DISPLAY_PO_FORMAT_SS,buffer,BOOL_2_YN(theGpo->OType[i]));
  }
  UserWrite("\n");

  return (0);
}

/****************************************************************************/
/*
   InitScalarFieldPlotObject_3D - Initialization of 3D scalar field object

   SYNOPSIS:
   static INT InitScalarFieldPlotObject_3D (PLOTOBJ *thePlotObj, INT argc,
   char **argv);

   PARAMETERS:
   .  thePlotObj -
   .  argc -
   .  argv -

   DESCRIPTION:
   This function makes an initialization of 3D scalar field object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

/****************************************************************************/
/*D
   EScalar3D - plot object for scalar grid functions (omit 3D when using it)

   DESCRIPTION:
   The EScalar plot object shows a contour or color plot of a scalar grid function.

   Mandatory options are:~
   .    $e~<vec~eval~proc>		- either specify <vec eval proc>
   .    $s~<vecdata~desc>		- or			 <vecdata desc> (standard eval proc, nodal values only)

   Possible options:~
   .    $f~<from>				- from value
   .    $t~<to>				- to value
   .    $d~<depth>				- integer for recursive depth (default 0,
                                                          CAUTION: may slow down dramatically!)
   .    $m COLOR|CONTOURS_EQ	- mode: color or equidistant contour lines
   .    $n~<levels>			- number of levels for contour lines
   .    $a~0..1                - contribution of ambient light to face intensity (back grid)



   KEYWORDS:
   graphics, plot, window, picture, plotobject, vecdesc, function
   D*/
/****************************************************************************/

static INT InitScalarFieldPlotObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  struct ElemScalarPlotObj3D *theEspo;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEspo = &(thePlotObj->theEspo);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V3_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  PO_USESCUT(thePlotObj) = YES;
  ret = ACTIVE;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theEspo->min                    = 0.0;
    theEspo->max                    = 1.0;
    theEspo->mode                   = PO_COLOR;
    theEspo->numOfContours  = 10;
    theEspo->AmbientLight   = 1.0;
  }

  /* set mode option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='m')
    {
      if (sscanf(argv[i],"m %s",buffer)!=1)
        break;
      if (strcmp(buffer,"COLOR") == 0)
        theEspo->mode = PO_COLOR;
      else if (strcmp(buffer,"CONTOURS_EQ") == 0)
        theEspo->mode = PO_CONTOURS_EQ;
      break;
    }

  /* set depth option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='d')
    {
      if (sscanf(argv[i],"d %d",&iValue)!=1)
        break;
      theEspo->depth = iValue;
      break;
    }
  if (theEspo->depth<0 || theEspo->depth>4)
  {
    UserWrite("depth is not valid\n");
    ret = NOT_ACTIVE;
  }

  /* set from option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='f')
    {
      if (sscanf(argv[i],"f %g",&fValue)!=1)
        break;
      theEspo->min = fValue;
      break;
    }

  /* set to option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='t')
    {
      if (sscanf(argv[i],"t %g",&fValue)!=1)
        break;
      theEspo->max = fValue;
      break;
    }
  if (theEspo->min >= theEspo->max )
  {
    UserWrite("minValue is bigger than maxValue\n");
    ret = NOT_ACTIVE;
  }

  /* set n option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='n')
    {
      if (sscanf(argv[i],"n %d",&iValue)!=1)
        break;
      if (iValue>1)
        theEspo->numOfContours = iValue;
      break;
    }
  if (theEspo->numOfContours <= 1 )
  {
    UserWrite("number of contours is smaller than 1\n");
    ret = NOT_ACTIVE;
  }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEspo),buffer);
      theEspo->EvalFct = GetElementValueEvalProc(buffer);
      break;
    }

  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEspo),buffer);
      if (theEspo->EvalFct == NULL)
        theEspo->EvalFct = GetElementValueEvalProc("nvalue");
      break;
    }

  for (i=1; i<argc; i++)
    if (argv[i][0]=='a')
    {
      if (sscanf(argv[i], "a %f", &fValue) != 1) break;
      theEspo->AmbientLight = fValue;
      break;
    }
  if (theEspo->AmbientLight < 0.0 || theEspo->AmbientLight > 1.0)
    theEspo->AmbientLight = 1.0;

  if (theEspo->EvalFct == NULL)
  {
    UserWrite("cannot find plot procedure\n");
    ret = NOT_ACTIVE;
  }

  /* do what is to do */
  if (theEspo->mode == PO_CONTOURS_EQ && ret == ACTIVE)
    for (i=0; i<theEspo->numOfContours; i++)
      theEspo->contValues[i] = theEspo->min + (DOUBLE)i * (theEspo->max - theEspo->min) / (DOUBLE)(theEspo->numOfContours-1);

  /* return */
  return (ret);
}

/****************************************************************************/
/*
   DisplayScalarFieldPlotObject_3D - Display content of 3D scalar field object

   SYNOPSIS:
   static INT DisplayScalarFieldPlotObject_3D (PLOTOBJ *thePlotObj);

   PARAMETERS:
   .  thePlotObj -

   DESCRIPTION:
   This function displays content of 3D scalar field object.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

static INT DisplayScalarFieldPlotObject_3D (PLOTOBJ *thePlotObj)
{
  struct ElemScalarPlotObj3D *theEspo;

  theEspo = &(thePlotObj->theEspo);

  /* print content */
  if (theEspo->EvalFct != NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theEspo->EvalFct));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc","---");

  UserWriteF(DISPLAY_PO_FORMAT_SS,"name",PO_NAME(theEspo));

  UserWriteF(DISPLAY_PO_FORMAT_SFF,"Range",(float)theEspo->min,(float)theEspo->max);
  UserWriteF(DISPLAY_PO_FORMAT_SI,"Depth",(int)theEspo->depth);
  switch (theEspo->mode)
  {
  case PO_COLOR :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"PlotMode","COLOR");
    break;
  case PO_CONTOURS_EQ :
    UserWriteF(DISPLAY_PO_FORMAT_SS,"PlotMode","CONTOURS_EQ");
    UserWriteF(DISPLAY_PO_FORMAT_SI,"NbOfCont",(int)theEspo->numOfContours);
  }
  UserWrite("\n");

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  InitVectorFieldPlotObject_3D									*/
/*																			*/
/* Purpose:   initialization of 3D vector field object						*/
/*																			*/
/* Input:	  PLOTOBJ *thePlotObj, INT argc, char **argv					*/
/*																			*/
/* Output:	  INT NOT_INIT, NOT_ACTIVE or ACTIVE							*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   EVector3D - plot object for (three-)vector grid functions (omit 3D when using it)

   DESCRIPTION:
   The EVector plot object shows a vector plot of a grid function. The vectors
   are plotted in a regular raster.

   Mandatory options are:~
   .    $e~<vec~eval~proc>		- either specify <vec eval proc>
   .    $s~<vecdata~desc>		- or			 <vecdata desc> (standard eval proc, nodal values only)

   Possible options:~
   .    $t~<to>				- to value
   .    $r~<raster>			- raster size
   .    $l~<cut~len>			- cut off len relative to raster size (default 1 to avoid overlap)
   .    $c~0|1					- cut off vectors off/on
   .    $p~0|1					- project vector off/on (?)
   .    $a~0..1                - contribution of ambient light to face intensity (back grid)

   KEYWORDS:
   graphics, plot, window, picture, plotobject, vecdesc, function
   D*/
/****************************************************************************/

static INT InitVectorFieldPlotObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC *theBVPDesc;
  struct ElemVectorPlotObj3D *theEvpo;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEvpo = &(thePlotObj->theEvpo);
  theBVPDesc = MG_BVPD(PO_MG(thePlotObj));
  V3_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  PO_USESCUT(thePlotObj) = YES;
  ret = ACTIVE;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theEvpo->max                    = 1.0;
    theEvpo->CutVector              = YES;
    theEvpo->ProjectVector  = YES;
    theEvpo->RasterSize     = PO_RADIUS(thePlotObj)/10.0;
    theEvpo->EvalFct                = NULL;
    theEvpo->CutLenFactor   = 0.9;
    theEvpo->AmbientLight   = 1.0;
  }

  /* set to option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='t')
    {
      if (sscanf(argv[i],"t %g",&fValue)!=1)
        break;
      theEvpo->max = fValue;
      break;
    }
  if (theEvpo->max <= 0.0)
  {
    UserWrite("maxValue is smaller than zero\n");
    ret =NOT_ACTIVE;
  }

  /* set CutLenFactor option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='l')
    {
      if (sscanf(argv[i],"l %g",&fValue)!=1)
        break;
      theEvpo->CutLenFactor = fValue;
      break;
    }
  if (theEvpo->CutLenFactor < 0.1 || theEvpo->CutLenFactor > 10)
  {
    UserWrite("CutLenFactor is not in [0.1,10]\n");
    ret =NOT_ACTIVE;
  }

  /* set rastersize option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='r')
    {
      if (sscanf(argv[i],"r %g",&fValue)!=1)
        break;
      theEvpo->RasterSize = fValue;
      break;
    }
  if (theEvpo->RasterSize <= 0.0)
  {
    UserWrite("RasterSize is smaller than zero\n");
    ret =NOT_ACTIVE;
  }

  /* set cut vector option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='c')
    {
      if (sscanf(argv[i],"c %d",&iValue)!=1)
        break;
      if (iValue==1)
        theEvpo->CutVector = YES;
      else if (iValue==0)
        theEvpo->CutVector = NO;
      break;
    }

  /* set project vector option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='p')
    {
      if (sscanf(argv[i],"p %d",&iValue)!=1)
        break;
      if (iValue==1)
        theEvpo->ProjectVector = YES;
      else if (iValue==0)
        theEvpo->ProjectVector = NO;
      break;
    }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEvpo),buffer);
      theEvpo->EvalFct = GetElementVectorEvalProc(buffer);
      break;
    }

  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      strcpy(PO_NAME(theEvpo),buffer);
      if (theEvpo->EvalFct == NULL)
        theEvpo->EvalFct = GetElementVectorEvalProc("nvector");
      break;
    }

  for (i=1; i<argc; i++)
    if (argv[i][0]=='a')
    {
      if (sscanf(argv[i], "a %f", &fValue) != 1) break;
      theEvpo->AmbientLight = fValue;
      break;
    }
  if (theEvpo->AmbientLight < 0.0 || theEvpo->AmbientLight > 1.0)
    theEvpo->AmbientLight = 1.0;

  if (theEvpo->EvalFct == NULL)
  {
    UserWrite("cannot find plot procedure\n");
    ret = NOT_ACTIVE;
  }

  /* return */
  return (ret);
}

/****************************************************************************/
/*																			*/
/* Function:  DisplayVectorFieldPlotObject_3D								*/
/*																			*/
/* Purpose:   display content of 3D vector field object                                         */
/*																			*/
/* Input:	  PLOTOBJ *thePlotObj											*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

static INT DisplayVectorFieldPlotObject_3D (PLOTOBJ *thePlotObj)
{
  struct ElemVectorPlotObj3D *theEvpo;

  theEvpo = &(thePlotObj->theEvpo);

  /* print content */
  if (theEvpo->EvalFct != NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theEvpo->EvalFct));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc","---");

  UserWriteF(DISPLAY_PO_FORMAT_SS,"name",PO_NAME(theEvpo));

  UserWriteF(DISPLAY_PO_FORMAT_SFF,"Range",0.0,(float)theEvpo->max);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"RasterSize",(float)theEvpo->RasterSize);
  UserWriteF(DISPLAY_PO_FORMAT_SF,"CutLenFactor",(float)theEvpo->CutLenFactor);

  if (theEvpo->CutVector == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"CutVector","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"CutVector","NO");
  if (theEvpo->ProjectVector == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"ProjectVector","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"ProjectVector","NO");
  UserWrite("\n");

  return (0);
}

#endif

/****************************************************************************/
/*D
   InitPlotObjTypes - Initialization of all available 'PLOTOBJTYPE's

   SYNOPSIS:
   INT InitPlotObjTypes (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   The function enrols all 'PLOTOBJTYPE's.

   .  2D-'PLOTOBJ' - 'EScalar', 'EVector' and 'Grid'. Options for the initialization,
                  see 'setplotobject'.
   .  3D-'PLOTOBJ' - 'EScalar', 'EVector' and 'Grid'. Options for the initialization,
                  see 'setplotobject'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT InitPlotObjTypes (void)
{
  PLOTOBJTYPE *thePOT;


  if ((thePOT=GetPlotObjType("Matrix"))  == NULL) return (1);
  thePOT->Dimension                               = TYPE_2D;
  thePOT->SetPlotObjProc                  = InitMatrixPlotObject;
  thePOT->DispPlotObjProc                 = DisplayMatrixPlotObject;

  /* set data and procedures of PLOTOBJTYPE */
        #ifdef __TWODIM__
  if ((thePOT=GetPlotObjType("EScalar")) == NULL) return (1);
  thePOT->Dimension                               = TYPE_2D;
  thePOT->SetPlotObjProc                  = InitScalarFieldPlotObject_2D;
  thePOT->DispPlotObjProc                 = DisplayScalarFieldPlotObject_2D;

  if ((thePOT=GetPlotObjType("EVector")) == NULL) return (1);
  thePOT->Dimension                               = TYPE_2D;
  thePOT->SetPlotObjProc                  = InitVectorFieldPlotObject_2D;
  thePOT->DispPlotObjProc                 = DisplayVectorFieldPlotObject_2D;

  if ((thePOT=GetPlotObjType("Grid"))    == NULL) return (1);
  thePOT->Dimension                               = TYPE_2D;
  thePOT->SetPlotObjProc                  = InitGridPlotObject_2D;
  thePOT->DispPlotObjProc                 = DisplayGridPlotObject_2D;

  if ((thePOT=GetPlotObjType("VecMat"))  == NULL) return (1);
  thePOT->Dimension                               = TYPE_2D;
  thePOT->SetPlotObjProc                  = InitVecMat_2D;
  thePOT->DispPlotObjProc                 = DisplayVecMat_2D;

  if ((thePOT=GetPlotObjType("Line")) == NULL) return (1);
  thePOT->Dimension                               = TYPE_2D;
  thePOT->SetPlotObjProc                  = InitLinePlotObject_2D;
  thePOT->DispPlotObjProc                 = DisplayLinePlotObject_2D;
        #endif

        #ifdef __THREEDIM__
  if ((thePOT=GetPlotObjType("EScalar")) == NULL) return (1);
  thePOT->Dimension                               = TYPE_3D;
  thePOT->SetPlotObjProc                  = InitScalarFieldPlotObject_3D;
  thePOT->DispPlotObjProc                 = DisplayScalarFieldPlotObject_3D;

  if ((thePOT=GetPlotObjType("EVector")) == NULL) return (1);
  thePOT->Dimension                               = TYPE_3D;
  thePOT->SetPlotObjProc                  = InitVectorFieldPlotObject_3D;
  thePOT->DispPlotObjProc                 = DisplayVectorFieldPlotObject_3D;

  if ((thePOT=GetPlotObjType("VecMat"))  == NULL) return (1);
  thePOT->Dimension                               = TYPE_3D;
  thePOT->SetPlotObjProc                  = InitVecMat_3D;
  thePOT->DispPlotObjProc                 = DisplayVecMat_3D;

  if ((thePOT=GetPlotObjType("Grid"))    == NULL) return (1);
  thePOT->Dimension                               = TYPE_3D;
  thePOT->SetPlotObjProc                  = InitGridObject_3D;
  thePOT->DispPlotObjProc                 = DisplayGridPlotObject_3D;

  /* maybe later: KJ
     if ((thePOT=GetPlotObjType("Domain"))  == NULL) return (1);
     thePOT->Dimension				= TYPE_3D;
     thePOT->SetPlotObjProc			= InitDomainPlotObject_3D;
     thePOT->DispPlotObjProc            = DisplayDomainPlotObject_3D;
   */
#endif

  return (0);
}




/****************************************************************************/
/*D
   InitWPM - Initialization

   SYNOPSIS:
   INT InitWPM (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function creates the 'ENVDIR'ectories `PlotObjTypes` and `UgWindows`.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT InitWPM (void)
{
  /* install the /PlotObjTypes directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitWPM","could not changedir to root");
    return(__LINE__);
  }
  thePlotObjTypesDirID = GetNewEnvDirID();
  if (MakeEnvItem("PlotObjTypes",thePlotObjTypesDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitWPM","could not install '/PlotObjTypes' dir");
    return(__LINE__);
  }
  thePlotObjTypesVarID = GetNewEnvVarID();

  /* install the /UgWindows directory */
  if (ChangeEnvDir("/")==NULL)
  {
    PrintErrorMessage('F',"InitWPM","could not changedir to root");
    return(__LINE__);
  }
  theUgWindowsDirID = GetNewEnvDirID();
  if (MakeEnvItem("UgWindows",theUgWindowsDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitWPM","could not install '/UgWindows' dir");
    return(__LINE__);
  }
  theUgWinDirID = GetNewEnvDirID();
  thePicVarID   = GetNewEnvVarID();

  return (0);
}
