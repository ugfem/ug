// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      sockcomm.c                                                    */
/*                                                                          */
/* Purpose:   routines and examples for berkeley socket usage               */
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
/* Remarks:   some of these routines have been copied from                  */
/*            W.R. Stevens, Unix Network Programming, Prentice Hall         */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __COMPILER__
#include "compiler.h"
#endif

#include "sockcomm.h"





/****************************************************************************/


char *cmd_text[] =
{
  "nn",
  "InitScreen",
  "WriteString",
  "GetNextUGEvent",
  "MousePosition",
  "MouseStillDown",

  "InitRemotePort",

  "OpenOutput",
  "CloseOutput",
  "ActivateOutput",
  "UpdateOutput",

  "Move",
  "Draw",
  "Polyline",
  "InversePolyline",
  "Polygon",
  "InversePolygon",
  "ErasePolygon",
  "Polymark",
  "Text",
  "CenteredText",
  "ClearViewPort",
  "SetLineWidth",
  "SetTextSize",
  "SetMarkerSize",
  "SetMarker",
  "SetColor",
  "SetPaletteEntry",
  "SetNewPalette",
  "GetPaletteEntry",
  "Flush"
};





/****************************************************************************/
/*                                                                          */
/* subroutines                                                              */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/

/* read n bytes from stream socket */

int SocketRead (int fd, char *ptr, int nbytes)
{
  int nleft, nread;

  nleft = nbytes;
  while (nleft>0)
  {
    nread = read(fd, ptr, nleft);
    if (nread<0)
      return(nread);                     /* error, return <0 */
    else if (nread==0)
      break;                             /* EOF */

    nleft -= nread;
    ptr   += nread;
  }

  return(nbytes-nleft);         /* return >=0 */
}

/****************************************************************************/

/* write n bytes to a stream socket */

int SocketWrite (int fd, char *ptr, int nbytes)
{
  int nleft, nwritten;

  nleft = nbytes;
  while (nleft>0)
  {
    nwritten = write(fd, ptr, nleft);
    if (nwritten<=0)
      return(nwritten);                     /* error */

    nleft -= nwritten;
    ptr   += nwritten;
  }

  return(nbytes-nleft);
}


/****************************************************************************/

/* read a line from a descriptor */

int SocketReadString (int fd, char *ptr, int maxlen)
{
  int n, rc;
  char c;

  for(n=1; n<maxlen; n++)
  {
    if ((rc=read(fd, &c, 1)) ==1)
    {
      *ptr++ = c;
      if (c=='\n')
        break;
    }
    else if (rc==0)
    {
      if (n==1)
        return(0);                           /* EOF, no data read */
      else
        break;                               /* EOF, some data was read */
    }
    else
      return(-1);                    /* error */
  }

  *ptr = 0;
  return(n);
}


/****************************************************************************/

void SocketWriteCmd (int sockfd, int cmd)
{
  if (write(sockfd, (char *)&cmd, sizeof(int)) < 0)
  {
    fprintf(stderr, "ug: write error on socket.\n");
  }
}


void SocketWriteINT (int sockfd, INT val)
{
  if (write(sockfd, (char *)&val, sizeof(val)) < 0)
  {
    fprintf(stderr, "ug: write error on socket.\n");
  }
}


void SocketWriteINTN (int sockfd, INT *val, int nval)
{
  if (write(sockfd, (char *)val, nval*sizeof(INT)) < 0)
  {
    fprintf(stderr, "ug: write error on socket.\n");
  }
}



void SocketWriteLong (int sockfd, long val)
{
  if (write(sockfd, (char *)&val, sizeof(long)) < 0)
  {
    fprintf(stderr, "ug: write error on socket.\n");
  }
}


void SocketWriteString (int sockfd, char *str)
{
  int n = strlen(str);

  if (SocketWrite(sockfd, str, n) != n)
  {
    fprintf(stderr, "ug: write error on socket.\n");
  }
}


void SocketWriteData (int sockfd, const char *data, int len)
{
  if (write(sockfd, (char *)&len, sizeof(int)) < 0)
  {
    fprintf(stderr, "ug: write error on socket.\n");
  }

  if (SocketWrite(sockfd, (char *)data, len) < 0)
  {
    fprintf(stderr, "ug: write error on socket.\n");
  }
}




INT SocketReadINT (int sockfd)
{
  INT v, n;

  n = SocketRead(sockfd, (char *)&v, sizeof(INT));
  return(v);
}


void SocketReadINTN (int sockfd, INT *val, int nval)
{
  INT n;

  n = SocketRead(sockfd, (char *)val, nval*sizeof(INT));
}


long SocketReadLong (int sockfd)
{
  long v;
  INT n;

  n = SocketRead(sockfd, (char *)&v, sizeof(long));
  return(v);
}


/****************************************************************************/
