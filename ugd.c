// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugd.c                                                         */
/*                                                                          */
/* Purpose:   main function and routines for ug-daemon                      */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   960820 kb  begin                                              */
/*                                                                          */
/* Remarks:                                                                 */
/*            there are some TODOs:                                         */
/*                                                                          */
/*            little/big endian conversion isn't implemented yet.           */
/*            the socket interface is weak, rework abstraction.             */
/*                                                                          */
/*                                                                          */
/****************************************************************************/


#ifdef RIF_SOCKETS  /* remote if must be enabled in ug.conf */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>

/* low module */
#include "initlow.h"
#include "misc.h"
#include "general.h"
#include "defaults.h"

#ifndef __COMPILER__
#include "compiler.h"
#endif

#include "devices.h"
#include "sockcomm.h"

#include "debug.h"

/****************************************************************************/


#define MAXLINE 1024



/****************************************************************************/


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


static char buf[MAXLINE];

static OUTPUTDEVICE *theOutputDevice;


/****************************************************************************/
/*                                                                          */
/* subroutines                                                              */
/*                                                                          */
/****************************************************************************/


static INT InitUgDaemon (int argc, char **argv)
{
  INT err;

  /* init the low module */
  if ((err=InitLow())!=0)
  {
    fprintf(stderr, "ugd: ERROR in InitUgDaemon while InitLow (line %d): called routine line %d\n",
            (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ugd\n");

    return (1);
  }


  /* init the devices module */
  if ((err=InitDevices(&argc,argv))!=0)
  {
    fprintf(stderr, "ugd: ERROR in InitUgDaemon while InitDevices (line %d): called routine line %d\n",
            (int) HiWrd(err), (int) LoWrd(err));
    printf ("aborting ugd\n");

    return (1);
  }

  return (0);
}



/****************************************************************************/

static void ExecInitScreen (int sockfd)
{
  INT len;
  int i, n, argc, error;
  char **argv;

  n = read(sockfd, (char *)&argc, sizeof(int));
  argv = (char **) malloc(sizeof(char *) * argc);

  for(i=0; i<argc; i++)
  {
    len  = SocketReadINT(sockfd);
    argv[i] = (char *) malloc(len+1);
    n = SocketRead(sockfd, argv[i], len);
    argv[i][len] = 0;
  }

  /* TODO. speicherverwaltung */

  InitUgDaemon(argc,argv);

  theOutputDevice = GetDefaultOutputDevice();
}


static void ExecOpenOutput (int sockfd)
{
  INT len;
  INT x, y, width, height;
  INT Global_LL[2], Global_UR[2], Local_LL[2], Local_UR[2];
  INT error;
  WINDOWID winID;

  x = SocketReadINT(sockfd);
  y = SocketReadINT(sockfd);
  width = SocketReadINT(sockfd);
  height = SocketReadINT(sockfd);

  len = SocketReadINT(sockfd);
  SocketRead(sockfd, buf, len);
  buf[len] = 0;


  winID = (*theOutputDevice->OpenOutput)(buf, x, y, width, height,
                                         Global_LL, Global_UR, Local_LL, Local_UR, &error);

  SocketWriteINTN(sockfd, Global_LL, 2);
  SocketWriteINTN(sockfd, Global_UR, 2);
  SocketWriteINTN(sockfd, Local_LL, 2);
  SocketWriteINTN(sockfd, Local_UR, 2);
  SocketWriteINT(sockfd, error);
  SocketWriteINT(sockfd, winID);
}



#define SKIP_NO_EVENTS   5

static void ExecGetNextUGEvent (int sockfd)
{
  EVENT ev;
  INT emask;
  int i;

  emask = SocketReadINT(sockfd);

  /* skip up to SKIP_NO_EVENTS NO_EVENTs */
  i = SKIP_NO_EVENTS;
  do
  {
    i--;
    GetNextUGEvent(&ev, emask);
  } while (i>0 && EVENT_TYPE(ev)==NO_EVENT);

  switch (EVENT_TYPE(ev))
  {
  case NO_EVENT :
    /*
       printf("event NO_EVENT\n");
     */
    SocketWriteData(sockfd, (char *)&ev, sizeof(NO_UGEVENT));
    break;

  case TERM_GOAWAY :
  case TERM_CMDKEY :
  case DOC_GOAWAY :
  case DOC_ACTIVATE :
  case DOC_DRAG :
  case DOC_GROW :
  case DOC_CHANGETOOL :
  case DOC_CONTENTCLICK :
  case DOC_UPDATE :
    /*
       printf("event %d\n", EVENT_TYPE(ev));
     */
    SocketWriteData(sockfd, (char *)&ev, sizeof(EVENT));
    break;

  case TERM_STRING :
    /*
       printf("event STRING %s\n", ev.TermString.String);
     */
    SocketWriteData(sockfd, (char *)&ev, sizeof(TERM_STRING_EVENT));
    break;

  default :
    fprintf(stderr, "ugd: unknown UGEvent type=%d!\n", ev.Type);
  }
}



/****************************************************************************/


static void ug_frontend (int sockfd)
{
  int n, cmd=-1;
  INT len;
  WINDOWID win;


  for (;;)
  {
                #ifdef DebugSockets
    {
      int info;
      n = read(sockfd, (char *)&info, sizeof(int));

      if (info!=MAGIC_COOKIE)
      {
        fprintf(stderr, "ugd: invalid cookie for cmd=#%d in debugmode!\n",
                cmd);
        fflush(stdout);
      }
    }
                #endif


    n = read(sockfd, (char *)&cmd, sizeof(int));
    if (n==0)
    {
      return;
    }

    else if (n<0)
    {
      fprintf(stderr, "ugd: ERROR in readline.\n");
      exit(1);
    }

    /*
       printf("cmd %s/%d\n", cmd_text[cmd], cmd);
       fflush(stdout);
     */

    switch (cmd)
    {
    case DC_InitScreen :
      ExecInitScreen(sockfd);
      break;

    case DC_WriteString :
      len = SocketReadINT(sockfd);
      n = SocketRead(sockfd, buf, len);
      buf[len] = 0;
      WriteString(buf);
      break;

    case DC_MousePosition :
    {
      INT pos[2];
      MousePosition(pos);
      SocketWriteINTN(sockfd, pos, 2);
      break;
    }

    case DC_MouseStillDown :
      SocketWriteINT(sockfd, MouseStillDown());
      break;

    case DC_GetNextUGEvent :
      ExecGetNextUGEvent(sockfd);
      break;

    case DC_InitRemotePort :
      SocketWriteLong(sockfd, theOutputDevice->black);
      SocketWriteLong(sockfd, theOutputDevice->white);
      SocketWriteLong(sockfd, theOutputDevice->red);
      SocketWriteLong(sockfd, theOutputDevice->green);
      SocketWriteLong(sockfd, theOutputDevice->blue);
      SocketWriteLong(sockfd, theOutputDevice->cyan);
      SocketWriteLong(sockfd, theOutputDevice->orange);
      SocketWriteLong(sockfd, theOutputDevice->yellow);
      SocketWriteLong(sockfd, theOutputDevice->darkyellow);
      SocketWriteLong(sockfd, theOutputDevice->magenta);
      SocketWriteINT(sockfd, theOutputDevice->hasPalette);
      SocketWriteLong(sockfd, theOutputDevice->spectrumStart);
      SocketWriteLong(sockfd, theOutputDevice->spectrumEnd);
      SocketWriteINT(sockfd, theOutputDevice->signx);
      SocketWriteINT(sockfd, theOutputDevice->signy);

      /*
                                      SocketWrite(sockfd, (char *)theOutputDevice, sizeof(OUTPUTDEVICE));
       */
      break;

    case DC_OpenOutput :
      ExecOpenOutput(sockfd);
      break;

    case DC_CloseOutput :
      win = SocketReadINT(sockfd);
      (*theOutputDevice->CloseOutput)(win);
      break;

    case DC_ActivateOutput :
      win = SocketReadINT(sockfd);
      (*theOutputDevice->ActivateOutput)(win);
      break;

    case DC_UpdateOutput :
    {
      INT tool;
      win = SocketReadINT(sockfd);
      tool = SocketReadINT(sockfd);

      (*theOutputDevice->UpdateOutput)(win,tool);
      break;
    }

    case DC_Move :
    {
      SHORT_POINT p;
      p.x = SocketReadINT(sockfd);
      p.y = SocketReadINT(sockfd);
      (*theOutputDevice->Move)(p);
      break;
    }

    case DC_Draw :
    {
      SHORT_POINT p;
      p.x = SocketReadINT(sockfd);
      p.y = SocketReadINT(sockfd);
      (*theOutputDevice->Draw)(p);
      break;
    }


    case DC_Polyline :
    {
      INT np = SocketReadINT(sockfd);
      SocketRead(sockfd, buf, sizeof(SHORT_POINT)*np);
      (*theOutputDevice->Polyline)((SHORT_POINT *)buf, np);
      break;
    }

    case DC_InversePolyline :
    {
      INT np = SocketReadINT(sockfd);
      SocketRead(sockfd, buf, sizeof(SHORT_POINT)*np);
      (*theOutputDevice->InversePolyline)((SHORT_POINT *)buf, np);
      break;
    }

    case DC_Polygon :
    {
      INT np = SocketReadINT(sockfd);
      SocketRead(sockfd, buf, sizeof(SHORT_POINT)*np);
      (*theOutputDevice->Polygon)((SHORT_POINT *)buf, np);
      break;
    }

    case DC_InversePolygon :
    {
      INT np = SocketReadINT(sockfd);
      SocketRead(sockfd, buf, sizeof(SHORT_POINT)*np);
      (*theOutputDevice->InversePolygon)((SHORT_POINT *)buf, np);
      break;
    }

    case DC_ErasePolygon :
    {
      INT np = SocketReadINT(sockfd);
      SocketRead(sockfd, buf, sizeof(SHORT_POINT)*np);
      (*theOutputDevice->ErasePolygon)((SHORT_POINT *)buf, np);
      break;
    }

    case DC_Polymark :
    {
      INT np = SocketReadINT(sockfd);
      SocketRead(sockfd, buf, sizeof(SHORT_POINT)*np);
      (*theOutputDevice->Polymark)(np, (SHORT_POINT *)buf);
      break;
    }

    case DC_InvPolymark :
    {
      INT np = SocketReadINT(sockfd);
      SocketRead(sockfd, buf, sizeof(SHORT_POINT)*np);
      (*theOutputDevice->InvPolymark)(np, (SHORT_POINT *)buf);
      break;
    }

    case DC_Text :
    {
      INT mode = SocketReadINT(sockfd);
      len = SocketReadINT(sockfd);
      n = SocketRead(sockfd, buf, len);
      buf[len] = 0;

      (*theOutputDevice->Text)(buf, mode);
      break;
    }

    case DC_CenteredText :
    {
      INT mode;
      SHORT_POINT p;
      p.x = SocketReadINT(sockfd);
      p.y = SocketReadINT(sockfd);
      mode = SocketReadINT(sockfd);

      len = SocketReadINT(sockfd);
      n = SocketRead(sockfd, buf, len);
      buf[len] = 0;

      (*theOutputDevice->CenteredText)(p, buf, mode);
      break;
    }

    case DC_ClearViewPort :
    {
      (*theOutputDevice->ClearViewPort)();
      break;
    }

    case DC_SetLineWidth :
    {
      short s = SocketReadINT(sockfd);
      (*theOutputDevice->SetLineWidth)(s);
      break;
    }

    case DC_SetTextSize :
    {
      short s = SocketReadINT(sockfd);
      (*theOutputDevice->SetTextSize)(s);
      break;
    }

    case DC_SetMarkerSize :
    {
      short s = SocketReadINT(sockfd);
      (*theOutputDevice->SetMarkerSize)(s);
      break;
    }

    case DC_SetMarker :
    {
      short s = SocketReadINT(sockfd);
      (*theOutputDevice->SetMarker)(s);
      break;
    }

    case DC_SetColor :
    {
      long v = SocketReadLong(sockfd);
      (*theOutputDevice->SetColor)(v);
      break;
    }

    case DC_SetPaletteEntry :
    {
      long index;
      short r, g, b;

      index = SocketReadLong(sockfd);
      r = SocketReadINT(sockfd);
      g = SocketReadINT(sockfd);
      b = SocketReadINT(sockfd);
      (*theOutputDevice->SetPaletteEntry)(index, r, g, b);

    }

    /* TODO andere Palette-funktionen! */


    case DC_Flush :
      (*theOutputDevice->Flush)();
      break;


    default :
      fprintf(stderr, "ugd: unknown command %s (#%d)!\n", cmd_text[cmd], cmd);
      fflush(stdout);
      break;
    }
  }
}


/****************************************************************************/

INT main (int argc, char **argv)
{
  int sockfd, newsockfd, clilen, childpid;
  struct sockaddr_in cli_addr, serv_addr;
  int port;
  int nodelay_flag, optlen=sizeof(nodelay_flag);


  /* open internet TCP stream socket */
  if ((sockfd=socket(AF_INET, SOCK_STREAM, 0)) < 0)
  {
    fprintf(stderr, "ugd: can't open TCP socket.\n");
    exit(1);
  }
  getsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, (char *)&nodelay_flag,&optlen);
  printf("TCP_NODELAY=%d\n");


  /* optional first argument is local port number */
  if (argc>1)
  {
    port = atoi(argv[1]);
  }
  else
  {
    port = SERV_TCP_PORT_DEFAULT;
  }

  if (port<SERV_TCP_PORT_MIN || port>SERV_TCP_PORT_MAX)
  {
    fprintf(stderr,"%s: port number must be between %d and %d!\n",
            argv[0], SERV_TCP_PORT_MIN, SERV_TCP_PORT_MAX);
    exit(-1);
  }



  /* bind local address */
  bzero((char *)&serv_addr, sizeof(serv_addr));
  serv_addr.sin_family  = AF_INET;
  serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  serv_addr.sin_port        = htons(port);

  if (bind(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
  {
    perror("ugd");
    exit(1);
  }

  listen(sockfd, 5);

  fprintf(stdout, "%s: listening to port %d.\n", argv[0], port);


  /*
          now clients can connect, enter main loop.

          this is a concurrent server, i.e., every time
          a socket connection is established, the server
          forks a child process and connects it to this
          socket. the child process is the ug-frontend.
   */

  for (;;)
  {
    clilen = sizeof(cli_addr);
    newsockfd = accept(sockfd,
                       (struct sockaddr *) &cli_addr,
                       &clilen);



    if (newsockfd<0)
    {
      fprintf(stderr, "ugd: error during accept.\n");
      exit(1);
    }

    if ((childpid=fork()) < 0)
    {
      fprintf(stderr, "ugd: error during fork.\n");
      exit(1);
    }

    else if (childpid==0)
    {
      /* this is the child process */

      /* close original socket */
      close(sockfd);

      /* child process: do work and suicide afterwards. */
      ug_frontend(newsockfd);
      printf("ugd: client has closed connection.\n");
      exit(0);
    }
    else
    {
      /* this is the server (ugd) process */
      printf("ugd: starting ug-frontend with pid=%d for client %s\n",
             childpid,
             InternetAddr((struct in_addr *) &cli_addr.sin_addr));
    }

    close(newsockfd);
  }

  /* we will never reach this line */
  return(0);
}


/****************************************************************************/

#endif
