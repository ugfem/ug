// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  remote.c														*/
/*																			*/
/* Purpose:   remote output interface, based on socket communication        */
/*																			*/
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   960820 kb  begin                                              */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


#ifdef RIF_SOCKETS   /* remote if must be enabled in ug.conf */

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/


#include "sockcomm.h"


/* standard C includes */
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>


/* interface includes */
#include "compiler.h"
#include "devices.h"
#include "initdev.h"
#include "debug.h"
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

#define MAXBUFSIZE 1024


#ifdef DebugSockets
#define BUFSTART        ptr=buf; BUFINT(MAGIC_COOKIE);
#else
#define BUFSTART    ptr=buf
#endif

#define BUFINT(i)   *(INT *)ptr=(INT)(i); ptr+=sizeof(INT)
#define BUFLONG(i)   *(long *)ptr=(long)(i); ptr+=sizeof(long)
#define BUFPOINT(p) BUFINT((p).x); BUFINT((p).y)
#define BUFTEXT(t)  { int l=strlen(t); BUFINT(l); memcpy(ptr,(t),l); ptr+=l; }
#define BUFEND          SocketWrite(theSocket, buf, ptr-buf);



/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);



/* current context */
static OUTPUTDEVICE *RemoteOutputDevice; /* outputdevice that has been initi */

static char buf[MAXBUFSIZE], *ptr;


static char *prog_name;                                                 /* my own name					*/
static unsigned int display_width;                              /* size of screen if needed     */
static unsigned int display_height;

static int theSocket;



/****************************************************************************/
/*																			*/
/* routines
   /*																			*/
/****************************************************************************/


static void IFMove (SHORT_POINT point)
{
  BUFSTART;
  BUFINT(DC_Move);
  BUFPOINT(point);
  BUFEND;
}

static void IFDraw (SHORT_POINT point)
{
  BUFSTART;
  BUFINT(DC_Draw);
  BUFPOINT(point);
  BUFEND;
}

static void IFPolyline (SHORT_POINT *points, INT n)
{
  BUFSTART;
  BUFINT(DC_Polyline);
  BUFINT(n);
  BUFEND;

  SocketWrite(theSocket, (char *)points, sizeof(SHORT_POINT)*n);
}

static void IFInversePolyline (SHORT_POINT *points, INT n)
{
  BUFSTART;
  BUFINT(DC_InversePolyline);
  BUFINT(n);
  BUFEND;

  SocketWrite(theSocket, (char *)points, sizeof(SHORT_POINT)*n);
}


static void IFPolygon (SHORT_POINT *points, INT n)
{
  BUFSTART;
  BUFINT(DC_Polygon);
  BUFINT(n);
  BUFEND;

  SocketWrite(theSocket, (char *)points, sizeof(SHORT_POINT)*n);
}


static void IFInversePolygon (SHORT_POINT *points, INT n)
{
  BUFSTART;
  BUFINT(DC_InversePolygon);
  BUFINT(n);
  BUFEND;

  SocketWrite(theSocket, (char *)points, sizeof(SHORT_POINT)*n);
}


static void IFErasePolygon (SHORT_POINT *points, INT n)
{
  BUFSTART;
  BUFINT(DC_ErasePolygon);
  BUFINT(n);
  BUFEND;

  SocketWrite(theSocket, (char *)points, sizeof(SHORT_POINT)*n);
}


static void IFPolymark (short n, SHORT_POINT *points)
{
  BUFSTART;
  BUFINT(DC_Polymark);
  BUFINT(n);
  BUFEND;

  SocketWrite(theSocket, (char *)points, sizeof(SHORT_POINT)*n);
}

static void IFInvPolymark (short n, SHORT_POINT *points)
{
  BUFSTART;
  BUFINT(DC_InvPolymark);
  BUFINT(n);
  BUFEND;

  SocketWrite(theSocket, (char *)points, sizeof(SHORT_POINT)*n);
}

static void IFText (const char *s, INT mode)
{
  BUFSTART;
  BUFINT(DC_Text);
  BUFINT(mode);
  BUFTEXT(s);        /* text am schluss, wegen alignment */
  BUFEND;
}


static void IFCenteredText (SHORT_POINT point, const char *s, INT mode)
{
  BUFSTART;
  BUFINT(DC_CenteredText);
  BUFPOINT(point);
  BUFINT(mode);
  BUFTEXT(s);        /* text am schluss, wegen alignment */
  BUFEND;
}

static void IFClearViewPort (void)
{
  BUFSTART;
  BUFINT(DC_ClearViewPort);
  BUFEND;
}


static void IFSetLineWidth (short w)
{
  BUFSTART;
  BUFINT(DC_SetLineWidth);
  BUFINT(w);
  BUFEND;
}

static void IFSetTextSize (short s)
{
  BUFSTART;
  BUFINT(DC_SetTextSize);
  BUFINT(s);
  BUFEND;
}

static void IFSetMarkerSize (short s)
{
  BUFSTART;
  BUFINT(DC_SetMarkerSize);
  BUFINT(s);
  BUFEND;
}

static void IFSetMarker (short s)
{
  BUFSTART;
  BUFINT(DC_SetMarker);
  BUFINT(s);
  BUFEND;
}

static void IFSetColor (long index)
{
  BUFSTART;
  BUFINT(DC_SetColor);
  BUFLONG(index);
  BUFEND;
}


static void IFSetPaletteEntry (long index, short r, short g, short b)
{
  BUFSTART;
  BUFINT(DC_SetPaletteEntry);
  BUFLONG(index);
  BUFINT(r);
  BUFINT(g);
  BUFINT(b);
  BUFEND;
}

static void IFSetNewPalette (long start, long count, short *r, short *g, short *b)
{
  BUFSTART;
  BUFINT(DC_SetNewPalette);
  /* TODO noch nicht fertig! */
  BUFEND;
}

static void IFGetPaletteEntry (long index, short *r, short *g, short *b)
{
  BUFSTART;
  BUFINT(DC_GetPaletteEntry);
  /* TODO noch nicht fertig! */
  BUFEND;
}

static void IFFlush (void)
{
  BUFSTART;
  BUFINT(DC_Flush);
  BUFEND;
}


/****************************************************************************/
/*
   InitXPort - implement basic drawing functions by X11

   SYNOPSIS:
   void InitXPort (OUTPUTDEVICE *thePort);

   PARAMETERS:
   .  thePort - PORT structure to initialize

   DESCRIPTION:
   This function  implements basic drawing functions by X11.

   RETURN VALUE:
   void
 */
/****************************************************************************/

void InitRemotePort (OUTPUTDEVICE *thePort)
{
  int l;

  BUFSTART;
  BUFINT(DC_InitRemotePort);
  BUFEND;

  /* init thePort with data from remote port */
  /*
          SocketRead(theSocket, (char *)thePort, sizeof(OUTPUTDEVICE));
   */
  thePort->black = SocketReadLong(theSocket);
  thePort->white = SocketReadLong(theSocket);
  thePort->red = SocketReadLong(theSocket);
  thePort->green = SocketReadLong(theSocket);
  thePort->blue = SocketReadLong(theSocket);
  thePort->cyan = SocketReadLong(theSocket);
  thePort->orange = SocketReadLong(theSocket);
  thePort->yellow = SocketReadLong(theSocket);
  thePort->darkyellow = SocketReadLong(theSocket);
  thePort->magenta = SocketReadLong(theSocket);
  thePort->hasPalette = SocketReadINT(theSocket);
  thePort->spectrumStart = SocketReadLong(theSocket);
  thePort->spectrumEnd = SocketReadLong(theSocket);
  thePort->signx = SocketReadINT(theSocket);
  thePort->signy = SocketReadINT(theSocket);

  thePort->PixelRatio = 1;


  /* now overwrite local values */

  /* init pointers to basic drawing functions */
  thePort->Move                   = IFMove;
  thePort->Draw                   = IFDraw;
  thePort->Polyline               = IFPolyline;
  thePort->InversePolyline= IFInversePolyline;
  thePort->Polygon                = IFPolygon;
  thePort->InversePolygon = IFInversePolygon;
  thePort->ErasePolygon   = IFErasePolygon;
  thePort->Polymark               = IFPolymark;
  thePort->InvPolymark    = IFInvPolymark;
  thePort->Text                   = IFText;
  thePort->CenteredText   = IFCenteredText;
  thePort->ClearViewPort  = IFClearViewPort;

  /* init pointers to set functions */
  thePort->SetLineWidth           = IFSetLineWidth;
  thePort->SetTextSize            = IFSetTextSize;
  thePort->SetMarker                      = IFSetMarker;
  thePort->SetMarkerSize          = IFSetMarkerSize;
  thePort->SetColor                       = IFSetColor;
  thePort->SetPaletteEntry        = IFSetPaletteEntry;
  thePort->SetNewPalette          = IFSetNewPalette;

  /* init pointers to miscellaneous functions */
  thePort->GetPaletteEntry        = IFGetPaletteEntry;
  thePort->Flush                          = IFFlush;
}


/*==========================================================================*/
/*																			*/
/* Outputdevice functions													*/
/*																			*/
/*==========================================================================*/


/****************************************************************************/
/*																			*/
/* Function:  OpenDocumentWindow											*/
/*																			*/
/* Purpose:   open a remote window and associate it with theView                        */
/*																			*/
/* Input:	  char *Windowtitle: window title								*/
/*			  UGWINDOW *theUgWindow: view for that window					*/
/*																			*/
/* Output:	  void *: pointer the window struct or							*/
/*					  NULL if an error occured								*/
/*																			*/
/****************************************************************************/


WINDOWID Remote_OpenOutput (const char *title, INT x, INT y, INT width, INT height, INT *Global_LL, INT *Global_UR, INT *Local_LL, INT *Local_UR, INT *error)
{
  WINDOWID win;

  BUFSTART;
  BUFINT(DC_OpenOutput);
  BUFINT(x);
  BUFINT(y);
  BUFINT(width);
  BUFINT(height);
  BUFTEXT(title);        /* text am schluss, wegen alignment */
  BUFEND;

  /* ... now the OpenOutput function is called on remote side */


  SocketReadINTN(theSocket, Global_LL, 2);
  SocketReadINTN(theSocket, Global_UR, 2);
  SocketReadINTN(theSocket, Local_LL, 2);
  SocketReadINTN(theSocket, Local_UR, 2);
  *error     = SocketReadINT(theSocket);
  win        = (WINDOWID)SocketReadINT(theSocket);

  /* return window ptr */
  return(win);
}


/****************************************************************************/
/*																			*/
/* Function:  CloseDocumentWindow											*/
/*																			*/
/* Purpose:   close the remote window associated with theView				*/
/*																			*/
/* Input:	  VIEW *theView: Port of that View will be inited				*/
/*																			*/
/* Output:	  INT: 0 if all was done well									*/
/*				   1 if an error ocurred									*/
/*																			*/
/****************************************************************************/

INT Remote_CloseOutput (WINDOWID win)
{
  BUFSTART;
  BUFINT(DC_CloseOutput);
  BUFINT(win);
  BUFEND;

  return(0);
}


/****************************************************************************/
/*																			*/
/* Functions: ActivateOutput												*/
/*																			*/
/* Purpose:   activate the window win										*/
/*																			*/
/* Input:	  WINDOWID win													*/
/*																			*/
/* Output:	  0 is OK														*/
/*																			*/
/****************************************************************************/

INT Remote_ActivateOutput (WINDOWID win)
{
  BUFSTART;
  BUFINT(DC_ActivateOutput);
  BUFINT(win);
  BUFEND;

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:   UpdateOutput                                                                                             */
/*																			*/
/* Purpose:    Draws all controls and highlights active tool				*/
/*																			*/
/* Input:	   GraphWindow *gwin											*/
/*																			*/
/* Output:	   0: OK														*/
/*			   1: error, could not complete                                                                 */
/*																			*/
/****************************************************************************/

INT Remote_UpdateOutput (WINDOWID win, INT tool)
{
  BUFSTART;
  BUFINT(DC_UpdateOutput);
  BUFINT(win);
  BUFINT(tool);
  BUFEND;

  return(0);
}

/****************************************************************************/
/*
   InitRemoteOutputDevice - Install output device 'remote'

   SYNOPSIS:
   OUTPUTDEVICE *InitOutputDevice (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function installs output device 'x11'.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      pointer to
   .n      NULL if an error occured.
 */
/****************************************************************************/

OUTPUTDEVICE *InitRemoteOutputDevice (void)
{
  /* create output device */
  if ((RemoteOutputDevice=CreateOutputDevice("screen"))==NULL)
    return(NULL);

  InitRemotePort(RemoteOutputDevice);
  RemoteOutputDevice->v.locked = 1;

  /* init output device 'screen' */
  RemoteOutputDevice->OpenOutput  = Remote_OpenOutput;
  RemoteOutputDevice->CloseOutput  = Remote_CloseOutput;
  RemoteOutputDevice->ActivateOutput  = Remote_ActivateOutput;
  RemoteOutputDevice->UpdateOutput  = Remote_UpdateOutput;


  printf("output device 'screen' for remote display created\n");

  return(RemoteOutputDevice);
}


/****************************************************************************/



/****************************************************************************/
/*
   GetScreenSize - Return sreen size

   SYNOPSIS:
   INT GetScreenSize (INT size[2]);

   PARAMETERS:
   .  size[2] - pointer to the size of screen

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

INT GetScreenSize (INT size[2])
{
  size[0] = display_width;
  size[1] = display_height;

  return (0);
}


/****************************************************************************/
/*
   GetNextUGEvent - Process an event from the system and pass it to ug

   SYNOPSIS:
   INT GetNextUGEvent (EVENT *theEvent, INT EventMask);

   PARAMETERS:
   .  theEvent - pointer to ug event

   .  EventMask -

   DESCRIPTION:
   This function processes an event from the system and passes it to ug if necessary.

   RETURN VALUE:
   INT
   .n     0 if no error occurred (ug or system)
   .n     1 if an error occurred (ug or system).
 */
/****************************************************************************/

INT GetNextUGEvent (EVENT *theEvent, INT Eventmask)
{
  INT len;

  BUFSTART;
  BUFINT(DC_GetNextUGEvent);
  BUFEND;
  /*
          SocketWriteData(theSocket, theEvent, sizeof(EVENT));
   */
  SocketWriteINT(theSocket, Eventmask);

  len = SocketReadINT(theSocket);
  if (len==0)
    return(0);

  /* an event has occured */
  SocketRead(theSocket, (char *)theEvent, len);
  return(0);
}


/****************************************************************************/
/*
   InitScreen - Init rest of GUI and return pointer to screen outputdevice

   SYNOPSIS:
   OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error);

   PARAMETERS:
   .  argcp - pointer to argument counter
   .  argv  - argument vector
   .  error - errorcode

   DESCRIPTION:
   This function inits rest of GUI and return ptr to screen outputdevice.

   RETURN VALUE:
   OUTPUTDEVICE *
   .n      POINTER if all is o.k.
   .n      NULL if an error occurred.

 */
/****************************************************************************/

OUTPUTDEVICE *InitScreen (int *argcp, char **argv, INT *error)
{
  struct sockaddr_in serv_addr;
  int i;
  OUTPUTDEVICE   *d;
  char buf[128], *hoststr, *portstr;
  int port;
  unsigned long inaddr;
  struct hostent *hp;
  int nodelay_flag, optlen=sizeof(nodelay_flag);


  /* connect procedure */
  prog_name = argv[0];
  strcpy(buf,"SHELLWINNAME");
  argv[0] = buf;


  /* init serv_addr socket structure */
  bzero((char *)&serv_addr, sizeof(serv_addr));
  serv_addr.sin_family  = AF_INET;

  if (*argcp<2)
  {
    fprintf(stderr,"%s: please specify host[:port] as first argument!\n",prog_name);
    exit(-1);
  }


  /*
          get first command line argument and interpret is as an
          address of the form:

          host[:port]

          where port is a port number between SERV_TCP_PORT_MIN and
          SERV_TCP_PORT_MAX, with SERV_TCP_PORT_DEFAULT being the default;
          host is the internet number or hostname of the
          machine where the daemon (ugd) is running.
   */
  hoststr = argv[1];
  portstr = strrchr(argv[1], ':');
  if (portstr==NULL)
  {
    /* no ':', no port number specified. use default. */
    port = SERV_TCP_PORT_DEFAULT;
  }
  else
  {
    /* port number specified. manipulate argument. */
    *portstr = 0;
    port = atoi(portstr+1);
    if (port<SERV_TCP_PORT_MIN || port>SERV_TCP_PORT_MAX)
    {
      fprintf(stderr,"%s: port number must be between %d and %d!\n",
              prog_name, SERV_TCP_PORT_MIN, SERV_TCP_PORT_MAX);
      exit(-1);
    }
  }

  /* lets try if hoststr is a dotted decimal number */
  if ((inaddr=inet_addr(hoststr)) != INADDR_NONE)
  {
    /* we got a xxx.yyy.zzz.zzz */
    bcopy((char *)&inaddr, (char *)&serv_addr.sin_addr, sizeof(inaddr));
  }
  else
  {
    /* no dotted decimal number, we try the hostname */
    if ((hp=gethostbyname(hoststr))==NULL)
    {
      fprintf(stderr,"%s: hostname error: %s (h_errno=%d)\n",
              prog_name, hoststr, h_errno);
      exit(-1);
    }

    bcopy(hp->h_addr, (char *)&serv_addr.sin_addr, hp->h_length);
  }

  serv_addr.sin_port        = htons(port);

  fprintf(stdout,"%s: connecting to %s, port %d ...\n",
          prog_name,
          InternetAddr((struct in_addr *)&serv_addr.sin_addr),
          port);


  /* now remove first argument from command line */
  for(i=1; i < (*argcp)-1; i++)
  {
    argv[i] = argv[i+1];
  }
  if (*argcp > 1) (*argcp)--;



  /* open internet TCP stream socket */
  if ((theSocket=socket(AF_INET, SOCK_STREAM, 0)) < 0)
  {
    fprintf(stderr,"%s: cannot open internet socket\n",prog_name);
    exit(-1);
  }
  getsockopt(theSocket, IPPROTO_TCP, TCP_NODELAY, (char *)&nodelay_flag,&optlen);
  printf("TCP_NODELAY=%d\n");


  /* connect to server */
  if (connect(theSocket, (struct sockaddr *)&serv_addr, sizeof(serv_addr))<0)
  {
    fprintf(stderr,"%s: cannot connect to ugd at %s\n",
            prog_name, hoststr);
    exit(-1);
  }



  /* now the socket connection is ready. */

  BUFSTART;
  BUFINT(DC_InitScreen);
  BUFINT((INT)*argcp);
  BUFEND;
  for(i=0; i<*argcp; i++)
    SocketWriteData(theSocket, argv[i], strlen(argv[i]));


  fprintf(stdout,"%s: connection to %s up and running.\n",
          prog_name, hoststr);


  d = InitRemoteOutputDevice();
  if (d==NULL) {*error=1; return(NULL);}

  *error = 0;

  return(d);
}

/****************************************************************************/

void ExitScreen (void)
{
  /* close socket connection */
  close(theSocket);
}


/****************************************************************************/
/*
   WriteString - write a string to a terminal window

   SYNOPSIS:
   void WriteString (const char *s);

   PARAMETERS:
   .  s -

   DESCRIPTION:
   This function writes a string to a terminal window.

   RETURN VALUE:
   void
 */
/****************************************************************************/

void WriteString (const char *s)
{
        #ifdef STDIF
  fputs(s,stdout);
        #ifdef Debug
  fflush(stdout);
        #endif
  /* printf("%s",s); */
        #else

  BUFSTART;
  BUFINT(DC_WriteString);
  BUFTEXT(s);
  BUFEND;

  return;
        #endif
}


/****************************************************************************/

void MousePosition (INT *point)
{
  BUFSTART;
  BUFINT(DC_MousePosition);
  BUFEND;
  SocketReadINTN(theSocket, point, 2);
}


/****************************************************************************/

INT MouseStillDown (void)
{
  BUFSTART;
  BUFINT(DC_MouseStillDown);
  BUFEND;
  return(SocketReadINT(theSocket));
}


/****************************************************************************/

void DrawInfoBox (WINDOWID win, const char *info)
{
  return;
}

INT WhichTool (WINDOWID win, const INT mouse[2], INT *tool)
{
  return (0);
}

#endif
