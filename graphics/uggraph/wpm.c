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

static COORD ex[3] = {1.0, 0.0, 0.0};
static COORD ey[3] = {0.0, 1.0, 0.0};
static COORD ez[3] = {0.0, 0.0, 1.0};

static INT theUgWindowsDirID;
static INT theUgWinDirID;

static INT thePicVarID;

static INT thePlotObjTypesDirID;
static INT thePlotObjTypesVarID;

/* data for CVS */
static char rcsid[] = "$Header$";

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
  if (theOutputDevice == NULL) return (NULL);

  /* allocate UgWindow envItem */
  if (ChangeEnvDir("/UgWindows") == NULL) return (NULL);
  if (strlen(UgWindowName)>=NAMESIZE || strlen(UgWindowName)<=1) return (NULL);
  if ((theWindow = (UGWINDOW *) MakeEnvItem(UgWindowName,theUgWinDirID,sizeof(UGWINDOW))) == NULL) return (NULL);

  /* open window on device and set sizes */
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

  /* set the other stuff */
  ENVITEM_LOCKED(theWindow)       = NO;
  UGW_NPIC(theWindow)             = 0;
  UGW_OUTPUTDEV(theWindow)        = theOutputDevice;
  UGW_CURRTOOL(theWindow)         = arrowTool;
  UGW_VALID(theWindow)            = NO;
  UGW_IFWINDOW(theWindow)         = winID;

  return (theWindow);
}

/****************************************************************************/
/*D
   UpdateUgWindow - Plot toolbox and if neccessary infobox

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
  MULTIGRID *theMG;
  INT size;
  char s[64];

  if (theUgWindow==NULL) return (0);

  /* update ugwindow */
  s[0] = '\0';
  if (EvalPicture==NULL)
  {
    strcpy(s,"---");
  }
  else
  if (PIC_UGW(EvalPicture)==theUgWindow)
  {
    if (VO_STATUS(PIC_VO(EvalPicture))==NOT_INIT)
      strcpy(s,"---");
    else
    {
      theMG = PIC_MG(EvalPicture);
      size = (HeapSize(MGHEAP(theMG))-HeapUsed(MGHEAP(theMG)))/1024;
      sprintf(s,"%d, %d k",(int)CURRENTLEVEL(theMG),(int)size);
    }
  }

  if ((*UGW_OUTPUTDEV(theUgWindow)->UpdateOutput)(UGW_IFWINDOW(theUgWindow),s,UGW_CURRTOOL(theUgWindow)))
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
  char buffer[128];

  sprintf(buffer,WPL_FORMAT,"","UgWindow","Picture","VO_Status","PlotObjType","PO_Status","Multigrid");
  UserWrite(buffer);
  sprintf(buffer,WPL_FORMAT,"","--------","-------","---------","-----------","---------","---------");
  UserWrite(buffer);
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
  char buffer[128];

  if (current) sprintf(buffer,WPL_FORMAT,"#",ENVITEM_NAME(theUgWindow),"","","","","");
  else sprintf(buffer,WPL_FORMAT,"",ENVITEM_NAME(theUgWindow),"","","","","");
  UserWrite(buffer);
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
  char b1[2], b2[11], b3[30], b4[30], b5[30], buffer[128];
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
  sprintf(buffer,WPL_FORMAT,b1,ENVITEM_NAME(theUgW),ENVITEM_NAME(thePicture),b2,b3,b4,b5);
  UserWrite(buffer);
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
  COORD a;

  for (thePicture=GetFirstPicture(theUgWindow); thePicture!=NULL; thePicture=GetNextPicture(thePicture))
  {
    a = ((COORD)(MousePosition[0]-PIC_GLL(thePicture)[0]))/((COORD)(PIC_GUR(thePicture)[0]-PIC_GLL(thePicture)[0]));
    if (a>0.0 && a<1.0)
    {
      a = ((COORD)(MousePosition[1]-PIC_GLL(thePicture)[1]))/((COORD)(PIC_GUR(thePicture)[1]-PIC_GLL(thePicture)[1]));
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
  /* change to directory */
  if (ChangeEnvDir("/PlotObjTypes")==NULL)
    return(NULL);

  if (size<sizeof(PLOTOBJTYPE))
    return (NULL);

  /* allocate structure */
  return ((PLOTOBJTYPE*) MakeEnvItem (PlotObjTypeName,thePlotObjTypesVarID,size));
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
  COORD ViewDirection[3], scalarPrd, help[3];

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
   SetView - Set the view

   SYNOPSIS:
   INT SetView (PICTURE *thePicture, const COORD *viewPoint,
   const COORD *targetPoint, const COORD *xAxis, const INT *perspective);

   PARAMETERS:
   .  thePicture - set view of that picture
   .  viewPoint - new view point
   .  targetPoint - new target point
   .  xAxis - new xAxis
   .  perspective - change of perspective

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

INT SetView (PICTURE *thePicture, const COORD *viewPoint, const COORD *targetPoint, const COORD *xAxis, const INT *perspective)
{
  VIEWEDOBJ *theViewedObj;
  PLOTOBJ *thePlotObj;
  COORD DefaultVP[3], DefaultVT[3], DefaultVTOld[3], DefaultPMP[3], DefaultPXD[3], DefaultPYD[3], DefaultPJ;
  COORD ViewDirection[3], ViewDirectionOld[3], CanvasRatio, RotationAxis[3];
  COORD angle, norm;
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
  CanvasRatio = ABS(((COORD)(PIC_GLL(thePicture)[1]-PIC_GUR(thePicture)[1]))/((COORD)(PIC_GLL(thePicture)[0]-PIC_GUR(thePicture)[0])));

  /* set values */
  switch (PO_DIM(thePlotObj))
  {
  case TYPE_2D :
    if (viewPoint != NULL || perspective != NULL) return (1);

    /* set values */
    if (ViewedObjNotInit)
    {
      V2_COPY(PO_MIDPOINT(thePlotObj),DefaultVT)                                                                                                /* set default values 2D	*/
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
    }
    else
    {
      V2_COPY(VO_VT(theViewedObj),DefaultVT)                                                                                                    /* save initial values 2D	*/
      V2_COPY(VO_PMP(theViewedObj),DefaultPMP)
      V2_COPY(VO_PXD(theViewedObj),DefaultPXD)
      V2_COPY(VO_PYD(theViewedObj),DefaultPYD)
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

    /* save values */
    V2_COPY(DefaultVT,VO_VT(theViewedObj))                                                                                              /* save values 2D					*/
    V2_COPY(DefaultVT,VO_PMP(theViewedObj))
    /*V2_COPY(DefaultPMP,VO_PMP(theViewedObj))*/
    V2_COPY(DefaultPXD,VO_PXD(theViewedObj))
    V2_COPY(DefaultPYD,VO_PYD(theViewedObj))

    /* set status of ViewedObj */
    VO_STATUS(theViewedObj) = ACTIVE;
    if (VO_PXD(theViewedObj)[0]==0.0 && VO_PXD(theViewedObj)[1]==0.0)
      VO_STATUS(theViewedObj) = NOT_ACTIVE;
    break;

  case TYPE_3D :
    if (viewPoint == NULL && ViewedObjNotInit)
    {
      UserWrite("initializing 3D-View, ViewPoint has to be specified\n");
      break;
    }

    /* set values */
    if (ViewedObjNotInit)
    {
      V3_COPY(viewPoint,DefaultVP)                                                                                                                              /* set default values 3D			*/
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
    V3_COPY(DefaultVP,VO_VP(theViewedObj))                                                                                                                      /* save values 3D					*/
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

/****************************************************************************/
/*D
   DisplayViewOfViewedObject - Display the view

   SYNOPSIS:
   INT DisplayViewOfViewedObject (PICTURE *thePicture);

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

INT DisplayViewOfViewedObject (PICTURE *thePicture)
{
  COORD width;
  char buffer[128];

  UserWrite("-----------------------\n");
  UserWrite(" Display of View of VO \n");
  UserWrite("-----------------------\n");

  switch (VO_STATUS(PIC_VO(thePicture)))
  {
  case NOT_INIT :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"VO_STATUS","NOT_INIT");
    UserWrite(buffer);
    return (0);
    break;
  case NOT_ACTIVE :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"VO_STATUS","NOT_ACTIVE");
    UserWrite(buffer);
    break;
  case ACTIVE :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"VO_STATUS","ACTIVE");
    UserWrite(buffer);
    break;
  default :
    return (1);
  }

  switch (VO_DIM(PIC_VO(thePicture)))
  {
  case NOT_DEFINED :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"Dim","NOT_DEFINED");
    UserWrite(buffer);
    break;
  case TYPE_2D :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"Dim","TYPE_2D");
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"Target",(float)(VO_VT(PIC_VO(thePicture))[0]),(float)(VO_VT(PIC_VO(thePicture))[1]));
    UserWrite(buffer);
    V2_EUKLIDNORM(VO_PXD(PIC_VO(thePicture)),width)
    width *= 2.0;
    sprintf(buffer,DISPLAY_PO_FORMAT_SF,"WinWidth",(float)width);
    UserWrite(buffer);
    break;
  case TYPE_3D :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"Dim","TYPE_3D");
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SFFF,"Observer",(float)(VO_VP(PIC_VO(thePicture))[0]),(float)(VO_VP(PIC_VO(thePicture))[1]),(float)(VO_VP(PIC_VO(thePicture))[2]));
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SFFF,"Target",(float)(VO_VT(PIC_VO(thePicture))[0]),(float)(VO_VT(PIC_VO(thePicture))[1]),(float)(VO_VT(PIC_VO(thePicture))[2]));
    UserWrite(buffer);
    V3_EUKLIDNORM(VO_PXD(PIC_VO(thePicture)),width)
    width *= 2.0;
    sprintf(buffer,DISPLAY_PO_FORMAT_SF,"WinWidth",(float)width);
    UserWrite(buffer);
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
  COORD delta[3],q[2];

  if (VO_STATUS(theVO) == NOT_INIT) return (0);

  q[0] = 1.0/(COORD)(Pix_UR_old[0]-Pix_LL_old[0]);
  q[1] = 1.0/(COORD)(Pix_UR_old[1]-Pix_LL_old[1]);

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
   INT Walk (PICTURE *thePicture, const COORD *vrsDelta);

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

INT Walk (PICTURE *thePicture, const COORD *vrsDelta)
{
  VIEWEDOBJ *theViewedObj;
  COORD VP[3], XD[3], YD[3], ZD[3];

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
  if (SetView(thePicture,VP,NULL,NULL,NULL)) return (1);

  return (0);
}

/****************************************************************************/
/*D
   RunAroundTargetPoint	- Modify the view by running around the midpoint

   SYNOPSIS:
   INT RunAroundTargetPoint (PICTURE *thePicture, COORD vrsDirectionAngle,
   COORD vrsAngle);

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

INT RunAroundTargetPoint (PICTURE *thePicture, COORD vrsDirectionAngle, COORD vrsAngle)
{
  VIEWEDOBJ *theViewedObj;
  COORD VP[3], RotationAxis[3], ViewDirection[3], TurnDirection[3];

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
  if (SetView(thePicture,VP,NULL,NULL,NULL)) return (1);

  return (0);
}

/****************************************************************************/
/*D
   Zoom	- zoom the view, i.e.: change size of projection plane

   SYNOPSIS:
   INT Zoom (PICTURE *thePicture, COORD factor);

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

INT Zoom (PICTURE *thePicture, COORD factor)
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
   INT DragProjectionPlane (PICTURE *thePicture, COORD vrsDeltaX,
   COORD vrsDeltaY);

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

INT DragProjectionPlane (PICTURE *thePicture, COORD vrsDeltaX, COORD vrsDeltaY)
{
  VIEWEDOBJ *theViewedObj;
  COORD DragVector[3], help[3];

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
   INT RotateProjectionPlane (PICTURE *thePicture, COORD vrsAngle);

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

INT RotateProjectionPlane (PICTURE *thePicture, COORD vrsAngle)
{
  VIEWEDOBJ *theViewedObj;
  COORD ViewDirection[3];

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
  if (SetView(thePicture,NULL,NULL,NULL,NULL)) return (1);

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
  char buffer[128];
  PLOTOBJTYPE *thePOT;

  if (thePlotObj==NULL) return (1);
  thePOT = PO_POT(thePlotObj);
  UserWrite("-----------------------\n");
  UserWrite(" Display of PlotObject \n");
  UserWrite("-----------------------\n");
  switch (PO_STATUS(thePlotObj))
  {
  case NOT_INIT :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"PO-NAME","---");
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"MG-NAME","---");
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"STATUS","NOT_INIT");
    UserWrite(buffer);
    return (0);
  case NOT_ACTIVE :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"PO-NAME",thePOT->v.name);
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"MG-NAME",MGNAME(PO_MG(thePlotObj)));
    UserWrite(buffer);
    if (PO_DIM(thePlotObj) == TYPE_2D)
      sprintf(buffer,DISPLAY_PO_FORMAT_SS,"STATUS","NOT_ACTIVE:2D");
    else
      sprintf(buffer,DISPLAY_PO_FORMAT_SS,"STATUS","NOT_ACTIVE:3D");
    UserWrite(buffer);
    break;
  case ACTIVE :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"PO-NAME",thePOT->v.name);
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"MG-NAME",MGNAME(PO_MG(thePlotObj)));
    UserWrite(buffer);
    if (PO_DIM(thePlotObj) == TYPE_2D)
      sprintf(buffer,DISPLAY_PO_FORMAT_SS,"STATUS","ACTIVE:2D");
    else
      sprintf(buffer,DISPLAY_PO_FORMAT_SS,"STATUS","ACTIVE:3D");
    UserWrite(buffer);
    break;
  }
  if (thePOT == NULL) return (0);

  switch (PO_DIM(thePlotObj))
  {
  case TYPE_2D :
    sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"MIDPOINT",(float)PO_MIDPOINT(thePlotObj)[0],(float)PO_MIDPOINT(thePlotObj)[1]);
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SF,"RADIUS",(float)PO_RADIUS(thePlotObj));
    UserWrite(buffer);
    break;
  case TYPE_3D :
    sprintf(buffer,DISPLAY_PO_FORMAT_SFFF,"MIDPOINT",(float)PO_MIDPOINT(thePlotObj)[0],(float)PO_MIDPOINT(thePlotObj)[1],(float)PO_MIDPOINT(thePlotObj)[2]);
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SF,"RADIUS",(float)PO_RADIUS(thePlotObj));
    UserWrite(buffer);
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
  float fValue;
  int iValue;
  char buffer[64];

  theMpo = &(thePlotObj->theMpo);
  theGrid = PO_MG(thePlotObj)->grids[PO_MG(thePlotObj)->currentLevel];
  if (theGrid == NULL) return (NOT_INIT);
  PO_MIDPOINT(thePlotObj)[0] = PO_MIDPOINT(thePlotObj)[1] = theGrid->nVector/2;
  PO_RADIUS(thePlotObj) = theGrid->nVector/2;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theMpo->min                     =-4.0;
    theMpo->max                     = 4.0;
    theMpo->log                     = FALSE;
    theMpo->conn            = TRUE;
    theMpo->extra           = FALSE;
    theMpo->EvalFct         = NULL;
    theMpo->Matrix          = NULL;
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
      theMpo->Matrix = GetMatSymbol(ENVITEM_NAME(MGFORMAT(PO_MG(thePlotObj))),buffer);
      if (theMpo->Matrix!=NULL)
        if (!MD_IS_SCALAR(SYM_MAT_DESC(theMpo->Matrix)))
          theMpo->Matrix = NULL;
      if (theMpo->Matrix == NULL)
      {
        UserWrite("cannot find scalar matrix symbol\n");
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
  UserWriteF(DISPLAY_PO_FORMAT_SF,"Thresh",(float)theMpo->thresh);

  /* print procedure name */
  if (theMpo->EvalFct!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"EvalProc",((ENVVAR*)theMpo->EvalFct)->name);

  if (theMpo->Matrix!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"Matrix",ENVITEM_NAME(theMpo->Matrix));

  return (0);
}

/****************************************************************************/
/*D
   SetCutPlane - Initialization/change of cut plane (3D-View)

   SYNOPSIS:
   static INT SetCutPlane (CUT *theCut, INT argc, char **argv);

   PARAMETERS:
   .  theCut - the structure 'CUT' to be initialized
   .  argc,~argv - arguments coming from the 'COMMANDS'

   DESCRIPTION:
   This function initializes or changes of cut plane (only 3D-View).You can specify
   the following 'COMMAND'-options

   .  P - the planepoint, three double-values are needed to specify a point on the
       cutplane
   .  N - the normal vector, three double-values are needed to specify the normale
       vector

   If the cutplane is not initialized, both options have to be specified together,
   if it is initialized, specifiing one option is possible to change the cutplane.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

static INT SetCutPlane (CUT *theCut, INT argc, char **argv)
{
  INT i, res, cppopt, cnpopt;
  float help[3];

  /* check if initialized */
  cppopt = cnpopt = 0;
  if (CUT_STATUS(theCut) != NOT_INIT)
    cppopt = cnpopt = 1;

  /* cut plane point option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='P')
    {
      res = sscanf(argv[i],"P %g %g %g",help, help+1, help+2);
      if (res!=3)
      {
        UserWrite("specify three values for cut plane point");
        return(1);
      }
      cppopt = 1;
      V3_COPY(help,theCut->PlanePoint);
      break;
    }

  /* cut normal point option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='N')
    {
      res = sscanf(argv[i],"N %g %g %g",help, help+1, help+2);
      if (res!=3)
      {
        UserWrite("specify three values for cut normal point");
        return(1);
      }
      cnpopt = 1;
      V3_COPY(help,theCut->PlaneNormal);
      break;
    }

  /* check how cut plane can now be (re)defined */
  CUT_STATUS(theCut) = NOT_INIT;
  if (cppopt && cnpopt)
  {
    if (V3_ISZERO(theCut->PlaneNormal))
    {
      UserWrite("cutting normal is (nearly) zero\n");
      CUT_STATUS(theCut) = NOT_ACTIVE;
    }
    else
      CUT_STATUS(theCut) = ACTIVE;
  }

  /* reset cut plane */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='R')
    {
      CUT_STATUS(theCut) = NOT_INIT;
      break;
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

static INT DisplayCutPlane (CUT *theCut)
{
  char buffer[128];

  /* display content */
  switch (CUT_STATUS(theCut))
  {
  case NOT_INIT :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"CUT STATUS","NOT_INIT");
    UserWrite(buffer);
    return (0);
  case NOT_ACTIVE :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"CUT STATUS","NOT_ACTIVE");
    UserWrite(buffer);
    break;
  case ACTIVE :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"CUT STATUS","ACTIVE");
    UserWrite(buffer);
    break;
  }
  sprintf(buffer,DISPLAY_PO_FORMAT_SFFF,"PlanePoint",(float)theCut->PlanePoint[0],(float)theCut->PlanePoint[1],(float)theCut->PlanePoint[2]);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SFFF,"PlaneNormal",(float)theCut->PlaneNormal[0],(float)theCut->PlaneNormal[1],(float)theCut->PlaneNormal[2]);
  UserWrite(buffer);

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

static INT InitGridPlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  struct GridPlotObj2D *theGpo;
  char buffer[64];
  INT i;
  int iValue;

  theGpo = &(thePlotObj->theGpo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theGpo->ElemColored     = NO;
    theGpo->WhichElem               = PO_ALL;
    theGpo->PlotBoundary    = YES;
    theGpo->PlotSegmentIDs  = NO;
    theGpo->PlotElemID              = NO;
    theGpo->PlotNodeID              = NO;
    theGpo->PlotNodes               = NO;
  }

  /* color mode */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='c')
    {
      if (sscanf(argv[i],"c %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->ElemColored = YES;
      else if (iValue==0)
        theGpo->ElemColored = NO;
      break;
    }
  for (i=1; i<argc; i++)
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


  /* set boundary option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='b')
    {
      if (sscanf(argv[i],"b %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->PlotBoundary = YES;
      else if (iValue==0)
        theGpo->PlotBoundary = NO;
      break;
    }

  /* set segID option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->PlotSegmentIDs = YES;
      else if (iValue==0)
        theGpo->PlotSegmentIDs = NO;
      break;
    }

  /* set refinement mark option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='r')
    {
      if (sscanf(argv[i],"r %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->PlotRefMarks = YES;
      else if (iValue==0)
        theGpo->PlotRefMarks = NO;
      break;
    }

  /* set elem id option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->PlotElemID = YES;
      else if (iValue==0)
        theGpo->PlotElemID = NO;
      break;
    }

  /* set node id option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='n')
    {
      if (sscanf(argv[i],"n %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->PlotNodeID = YES;
      else if (iValue==0)
        theGpo->PlotNodeID = NO;
      break;
    }

  /* plot node marker option */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='m')
    {
      if (sscanf(argv[i],"m %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->PlotNodes = YES;
      else if (iValue==0)
        theGpo->PlotNodes = NO;
      break;
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
  char buffer[128];

  theGpo = &(thePlotObj->theGpo);

  /* print content */
  if (theGpo->PlotBoundary == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"BND","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"BND","NO");
  UserWrite(buffer);

  if (theGpo->PlotSegmentIDs == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"SegID","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"SegID","NO");
  UserWrite(buffer);

  if (theGpo->PlotNodes == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"Node markers","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"Node markers","NO");
  UserWrite(buffer);

  if (theGpo->PlotRefMarks == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"ref marks","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"ref marks","NO");
  UserWrite(buffer);

  if (theGpo->PlotElemID == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"ElemID","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"ElemID","NO");
  UserWrite(buffer);

  if (theGpo->PlotNodeID == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"NodeID","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"NodeID","NO");
  UserWrite(buffer);

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
  UserWrite(buffer);

  if (theGpo->ElemColored == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"COLORED","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"COLORED","NO");
  UserWrite(buffer);

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

static INT InitVecMat_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  FORMAT *theFormat;
  struct VecMatPlotObj2D *theVmo;
  char name[NAMELEN];
  INT i,rt,ct;
  int iValue;

  theVmo = &(thePlotObj->theVmo);
  theFormat = MGFORMAT(PO_MG(thePlotObj));
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theVmo->Marker                          = NO;
    theVmo->Type[NODEVECTOR]        = YES;
    theVmo->Type[EDGEVECTOR]        = YES;
    theVmo->Type[ELEMVECTOR]        = YES;
    theVmo->Connections                     = YES;
    theVmo->Extra                           = NO;
    theVmo->Idx                                     = NO;
    theVmo->Order                           = 0;
    theVmo->Dependency                      = NO;
    theVmo->ConnectVectors          = NO;
    theVmo->Boundary                        = YES;
    theVmo->vs                                      = NULL;
    theVmo->ms                                      = NULL;
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
      if (strstr(argv[i],"nd")!=NULL) theVmo->Type[NODEVECTOR] = YES;
      else theVmo->Type[NODEVECTOR] = NO;
      if (strstr(argv[i],"ed")!=NULL) theVmo->Type[EDGEVECTOR] = YES;
      else theVmo->Type[EDGEVECTOR] = NO;
      if (strstr(argv[i],"el")!=NULL) theVmo->Type[ELEMVECTOR] = YES;
      else theVmo->Type[ELEMVECTOR] = NO;
      break;

    case 'c' :
      if (strcmp(argv[i],"connect")==0)
      {
        theVmo->ConnectVectors = YES;
        break;
      }
      if (sscanf(argv[i],"c %d",&iValue)!=1)
        break;
      if              (iValue==1) theVmo->Connections = YES;
      else if (iValue==0) theVmo->Connections = NO;
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
      if ((sscanf(argv[i],"V %s",name)!=1) || ((theVmo->vs=GetVecSymbol(ENVITEM_NAME(theFormat),name))==NULL))
      {
        UserWrite("no vector symbol specified, vec data switched off\n");
        theVmo->vs = NULL;
      }
      break;

    case 'M' :
      if ((sscanf(argv[i],"M %s",name)!=1) || ((theVmo->ms=GetMatSymbol(ENVITEM_NAME(theFormat),name))==NULL))
      {
        UserWrite("no matrix symbol specified, mat data switched off\n");
        theVmo->ms = NULL;
      }
      break;
    }

  if (theVmo->ConnectVectors)
  {
    theVmo->Connections = NO;
    theVmo->Extra           = NO;
  }

  /* check compatibility of vec and mat desc */
  if (theVmo->vs || theVmo->ms)
    for (rt=0; rt<NVECTYPES; rt++)
      if (theVmo->Type[rt])
      {
        if (theVmo->vs)
          if (!VD_ISDEF_IN_TYPE(SYM_VEC_DESC(theVmo->vs),rt))
          {
            UserWrite("vec desc does not include types of specified types\n");
            return (NOT_ACTIVE);
          }
        if (theVmo->ms)
          for (ct=0; ct<NVECTYPES; ct++)
            if (theVmo->Type[ct])
            {
              if (!MD_ISDEF_IN_RT_CT(SYM_MAT_DESC(theVmo->ms),rt,ct))
              {
                UserWrite("mat desc does not include column types of specified types\n");
                return (NOT_ACTIVE);
              }
              if (theVmo->vs)
                if (VD_NCMPS_IN_TYPE(SYM_VEC_DESC(theVmo->vs),ct)!=MD_ROWS_IN_RT_CT(SYM_MAT_DESC(theVmo->ms),rt,ct))
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
  struct VecMatPlotObj2D *theVmo;

  theVmo = &(thePlotObj->theVmo);

  /* print content */
  if (theVmo->Marker == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"marker","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"marker","NO");

  if (theVmo->Type[NODEVECTOR] == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"node-vecs","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"node-vecs","NO");

  if (theVmo->Type[EDGEVECTOR] == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"edge-vecs","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"edge-vecs","NO");

  if (theVmo->Type[ELEMVECTOR] == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"elem-vecs","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"elem-vecs","NO");

  if (theVmo->Idx == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"index","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"index","NO");

  if (theVmo->Connections == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"connections","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"connections","NO");

  if (theVmo->Extra == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"extra","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"extra","NO");

  switch (theVmo->Order)
  {
  case 1 : UserWriteF(DISPLAY_PO_FORMAT_SS,"order","1: vector order"); break;
  case 2 : UserWriteF(DISPLAY_PO_FORMAT_SS,"order","2: vector order, cuts black"); break;
  case 3 : UserWriteF(DISPLAY_PO_FORMAT_SS,"order","3: blockvector order"); break;
  }

  if (theVmo->Dependency == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"dependency","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"dependency","NO");

  if (theVmo->ConnectVectors == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"connect","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"connect","NO");

  if (theVmo->Boundary == YES)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"boundary","YES");
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"boundary","NO");

  if (theVmo->vs!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"vec data",SYM_NAME(theVmo->vs));
  else
    UserWriteF(DISPLAY_PO_FORMAT_SS,"vec data","NO");

  if (theVmo->ms!=NULL)
    UserWriteF(DISPLAY_PO_FORMAT_SS,"mat data",SYM_NAME(theVmo->ms));
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

static INT InitScalarFieldPlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  struct ElemScalarPlotObj2D *theEspo;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEspo = &(thePlotObj->theEspo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
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
    sprintf(buffer,"number of contours is greater than the limit (%d)",PO_MAXCONTOURS);
    PrintErrorMessage('E',"InitScalarFieldPlotObject_2D",buffer);
    ret = NOT_ACTIVE;
  }

  /* get plot procedure */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %s",buffer)!=1)
        break;
      if (strlen(buffer)>=NAMESIZE) break;
      theEspo->EvalFct = GetElementValueEvalProc(buffer);
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
  char buffer[128];

  theEspo = &(thePlotObj->theEspo);

  /* print content */
  if (theEspo->EvalFct != NULL)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theEspo->EvalFct));
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc","---");
  UserWrite(buffer);

  if (theEspo->PlotGrid == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"Grid","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"Grid","NO");
  UserWrite(buffer);

  sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"Range",(float)theEspo->min,(float)theEspo->max);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SI,"Depth",(int)theEspo->depth);
  UserWrite(buffer);
  switch (theEspo->mode)
  {
  case PO_COLOR :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"PlotMode","COLOR");
    UserWrite(buffer);
    break;
  case PO_CONTOURS_EQ :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"PlotMode","CONTOURS_EQ");
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SI,"NbOfCont",(int)theEspo->numOfContours);
    UserWrite(buffer);
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

static INT InitVectorFieldPlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  struct ElemVectorPlotObj2D *theEvpo;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEvpo = &(thePlotObj->theEvpo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
  V2_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  ret = ACTIVE;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theEvpo->max                    = 1.0;
    theEvpo->CutVectors     = YES;
    theEvpo->RasterSize     = PO_RADIUS(thePlotObj)/10.0;
    theEvpo->CutLenFactor   = 1.0;
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
      theEvpo->EvalFct = GetElementVectorEvalProc(buffer);
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
  char buffer[128];

  theEvpo = &(thePlotObj->theEvpo);

  /* print content */
  sprintf(buffer,DISPLAY_PO_FORMAT_SF,"maxValue",(float)theEvpo->max);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SF,"RasterSize",(float)theEvpo->RasterSize);
  UserWrite(buffer);
  if (theEvpo->CutVectors == YES)
  {
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"CutVectors","YES");
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SF,"CutLenFactor",(float)theEvpo->CutLenFactor);
    UserWrite(buffer);
  }
  else
  {
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"CutVectors","NO");
    UserWrite(buffer);
  }

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

static INT InitLinePlotObject_2D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  struct LinePlotObj2D *theLpo;
  INT i, ret;
  int iValue;
  float fValue[2];
  COORD dist;
  char buffer[128];

  theLpo = &(thePlotObj->theLpo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
  PO_MIDPOINT(thePlotObj)[0] = PO_MIDPOINT(thePlotObj)[1] = 0.5;
  PO_RADIUS(thePlotObj) = 0.70711;
  ret = ACTIVE;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theLpo->min                                                             = 0.0;
    theLpo->max                                                             = 1.0;
    theLpo->left[0] = theLpo->left[1]               = 0.0;
    theLpo->right[0] = theLpo->right[1]     = 1.0;
    theLpo->color                                                   = 0.0;
    theLpo->aspectratio                                             = 1.0;
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
      theLpo->EvalFct = GetElementValueEvalProc(buffer);
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
  char buffer[128];

  theLpo = &(thePlotObj->theLpo);

  /* print content */
  if (theLpo->EvalFct != NULL)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theLpo->EvalFct));
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc","---");
  UserWrite(buffer);

  sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"Range",(float)theLpo->min,(float)theLpo->max);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"left",(float)theLpo->left[0],(float)theLpo->left[1]);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"right",(float)theLpo->right[0],(float)theLpo->right[1]);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SF,"color",(float)theLpo->color);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SF,"asp. ratio",(float)theLpo->aspectratio);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SI,"Depth",(int)theLpo->depth);
  UserWrite(buffer);

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
  BVP_DESC theBVPDesc;
  struct DomainPlotObj3D *theDpo;
  INT i;
  int iValue;

  theDpo = &(thePlotObj->theDpo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
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
  char buffer[128];

  theDpo = &(thePlotObj->theDpo);

  /* print content */
  sprintf(buffer,"NbOfSteps  =%d\n",(int)theDpo->NbOfSteps);
  UserWrite(buffer);

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

static INT InitGridObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  struct GridPlotObj3D *theGpo;
  CUT *theCut;
  INT i;
  int iValue;
  char buffer[64];

  theGpo = &(thePlotObj->theGpo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
  V3_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  theCut = &(theGpo->theCut);

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theGpo->ElemColored     = NO;
    theGpo->WhichElem               = PO_ALL;
    CUT_STATUS(theCut)              = NOT_INIT;
  }

  /* color mode */
  for (i=1; i<argc; i++)
    if (argv[i][0]=='c')
    {
      if (sscanf(argv[i],"c %d",&iValue)!=1)
        break;
      if (iValue==1)
        theGpo->ElemColored = YES;
      else if (iValue==0)
        theGpo->ElemColored = NO;
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

  /* (re)define cut plane */
  if (SetCutPlane(theCut,argc,argv)) return (1);

  if (CUT_STATUS(theCut)==NOT_ACTIVE) return (NOT_ACTIVE);
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
  CUT *theCut;
  char buffer[128];

  theGpo = &(thePlotObj->theGpo);
  theCut = &(theGpo->theCut);

  /* print content */
  if (theGpo->ElemColored == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"COLORED","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"COLORED","NO");
  UserWrite(buffer);
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
  UserWrite(buffer);
  UserWrite("\n");

  /* print content of cut plane */
  if (DisplayCutPlane(theCut)) return (1);

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

static INT InitScalarFieldPlotObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  CUT *theCut;
  struct ElemScalarPlotObj3D *theEspo;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEspo = &(thePlotObj->theEspo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
  V3_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  theCut = &(theEspo->theCut);
  ret = ACTIVE;

  /* defaults */
  if (PO_STATUS(thePlotObj)==NOT_INIT)
  {
    theEspo->min                    = 0.0;
    theEspo->max                    = 1.0;
    theEspo->mode                   = PO_COLOR;
    theEspo->numOfContours  = 10;
    CUT_STATUS(theCut)              = NOT_INIT;
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
      theEspo->EvalFct = GetElementValueEvalProc(buffer);
      break;
    }
  if (theEspo->EvalFct == NULL)
  {
    UserWrite("cannot find plot procedure\n");
    ret = NOT_ACTIVE;
  }

  /* (re)define cut plane */
  if (SetCutPlane(theCut,argc,argv)) return (1);
  if (CUT_STATUS(theCut)==NOT_ACTIVE)
    ret = NOT_ACTIVE;

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
  CUT *theCut;
  struct ElemScalarPlotObj3D *theEspo;
  char buffer[128];

  theEspo = &(thePlotObj->theEspo);
  theCut = &(theEspo->theCut);

  /* print content */
  if (theEspo->EvalFct != NULL)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theEspo->EvalFct));
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc","---");
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"Range",(float)theEspo->min,(float)theEspo->max);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SI,"Depth",(int)theEspo->depth);
  UserWrite(buffer);
  switch (theEspo->mode)
  {
  case PO_COLOR :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"PlotMode","COLOR");
    UserWrite(buffer);
    break;
  case PO_CONTOURS_EQ :
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"PlotMode","CONTOURS_EQ");
    UserWrite(buffer);
    sprintf(buffer,DISPLAY_PO_FORMAT_SI,"NbOfCont",(int)theEspo->numOfContours);
    UserWrite(buffer);
  }
  UserWrite("\n");

  /* print content of cut plane */
  if (DisplayCutPlane(theCut)) return (1);

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

static INT InitVectorFieldPlotObject_3D (PLOTOBJ *thePlotObj, INT argc, char **argv)
{
  BVP_DESC theBVPDesc;
  struct ElemVectorPlotObj3D *theEvpo;
  CUT *theCut;
  char buffer[64];
  INT i, ret;
  int iValue;
  float fValue;

  theEvpo = &(thePlotObj->theEvpo);
  if (BVP_GetBVPDesc(MG_BVP(PO_MG(thePlotObj)),&theBVPDesc)) return (NOT_INIT);
  V3_COPY(BVPD_MIDPOINT(theBVPDesc),PO_MIDPOINT(thePlotObj))
  PO_RADIUS(thePlotObj) = BVPD_RADIUS(theBVPDesc);
  theCut = &(theEvpo->theCut);
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
      theEvpo->EvalFct = GetElementVectorEvalProc(buffer);
      break;
    }
  if (theEvpo->EvalFct == NULL)
  {
    UserWrite("cannot find plot procedure\n");
    ret = NOT_ACTIVE;
  }

  /* (re)define cut plane */
  if (SetCutPlane(theCut,argc,argv)) return (1);
  if (CUT_STATUS(theCut)==NOT_ACTIVE)
    ret = NOT_ACTIVE;

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
  CUT *theCut;
  char buffer[128];

  theEvpo = &(thePlotObj->theEvpo);
  theCut = &(theEvpo->theCut);

  /* print content */
  if (theEvpo->EvalFct != NULL)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc",ENVITEM_NAME(theEvpo->EvalFct));
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"EvalProc","---");
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SFF,"Range",0.0,(float)theEvpo->max);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SF,"RasterSize",(float)theEvpo->RasterSize);
  UserWrite(buffer);
  sprintf(buffer,DISPLAY_PO_FORMAT_SF,"CutLenFactor",(float)theEvpo->CutLenFactor);
  UserWrite(buffer);

  if (theEvpo->CutVector == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"CutVector","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"CutVector","NO");
  UserWrite(buffer);
  if (theEvpo->ProjectVector == YES)
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"ProjectVector","YES");
  else
    sprintf(buffer,DISPLAY_PO_FORMAT_SS,"ProjectVector","NO");
  UserWrite(buffer);
  UserWrite("\n");

  /* print content of cut plane */
  if (DisplayCutPlane(theCut)) return (1);

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

  /* set data and procedures of PLOTOBJTYPE */
        #ifdef __TWODIM__

  if ((thePOT=GetPlotObjType("Matrix"))  == NULL) return (1);
  thePOT->Dimension                               = TYPE_2D;
  thePOT->SetPlotObjProc                  = InitMatrixPlotObject;
  thePOT->DispPlotObjProc                 = DisplayMatrixPlotObject;

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
  if ((thePOT=GetPlotObjType("Matrix"))  == NULL) return (1);
  thePOT->Dimension                               = TYPE_3D;
  thePOT->SetPlotObjProc                  = InitMatrixPlotObject;
  thePOT->DispPlotObjProc                 = DisplayMatrixPlotObject;

  if ((thePOT=GetPlotObjType("EScalar")) == NULL) return (1);
  thePOT->Dimension                               = TYPE_3D;
  thePOT->SetPlotObjProc                  = InitScalarFieldPlotObject_3D;
  thePOT->DispPlotObjProc                 = DisplayScalarFieldPlotObject_3D;

  if ((thePOT=GetPlotObjType("EVector")) == NULL) return (1);
  thePOT->Dimension                               = TYPE_3D;
  thePOT->SetPlotObjProc                  = InitVectorFieldPlotObject_3D;
  thePOT->DispPlotObjProc                 = DisplayVectorFieldPlotObject_3D;

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
