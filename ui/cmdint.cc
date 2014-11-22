// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  cmdint.c                                                                                                              */
/*																			*/
/* Purpose:   command interpreter for ug 2.0								*/
/*																			*/
/* Author:	  Peter Bastian, Nicolas Neuss									*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de				                */
/*																			*/
/* History:   18.02.92 begin, ug version 2.0								*/
/*			  10.05.92 begin, ug interpreter by Nicolas Neuss				*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment cmd
#endif

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <config.h>
#include <ctype.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>

/* low module */
#include "ugtypes.h"
#include "fileopen.h"
#include "heaps.h"
#include "ugenv.h"
#include "misc.h"
#include "defaults.h"
#include "general.h"

/* dev module */
#include "ugdevices.h"
#include "debug.h"

/* ui module */
#include "ugstruct.h"
#include "cmdint.hh"
#include "uginterface.hh"
#include "cmdline.hh"

/** \todo conflict with hierarchy */
#include "initug.h"

#ifdef __GUI__
#include <X11/Xlib.h>
#include "xmain.h"
#endif

#include "ppif_namespace.h"

/* glibc: no prototype in math.h unless _ISOC9X_SOURCE is defined */
long int lrint(double x);

USING_UG_NAMESPACES
  USING_PPIF_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define VERSION_STRING   "This is " UG_VERSION " from $Date$\n"

/* for interpreter */
#define DONE                            0
#define qErr                            -1

#define BUFSIZE                         32
#define MAXREPEAT                       16
#define MAXTOKENLENGTH          64
#define MAXSTRINGSIZE           256
#define PROGRAMBUFSIZE          8000
#define FILEBUFSIZE                     1000

#define MULTIPLE                        256
#define SKIPMODE                        512

#define STATUSCODE                      255
#define IF                                      1
#define ELSE                            2
#define REPEAT                          4

#define NUMBERID                        1
#define ALPHAID                         2
#define LSTRINGID                       3
#define EMPTYID                         4

#define Boolean int

#define ISALPHA(c)                      (isalpha((int)c) || (((char) c) == '_') || (((char) c) == ':'))
#define ISNUMBER(c)                      isdigit((int)c)

static char blanks[]                    = " \t\n";
/* static char normalBlanks[]		= " \t"; */
static char seperators[]                = ";{}";
/* static char tokenseperators[]	= " \n\t;{}()"; */
static char terminators[]               = ";}";         /* terminators for commands */

static INT scriptpaths_set=FALSE;
static INT dontexit=FALSE;      /* if TRUE set ':cmdstatus' rather than exiting	*/
static INT UseWithPerl=0;

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

struct stringOp {
  INT type;
  const char *sptr;
};

struct lstringOp {
  INT type;
  const char *sptr;
  int length;
};

struct realOp {
  INT type;
  DOUBLE value;
};

union Operand {
  struct lstringOp lo;
  struct stringOp so;
  struct realOp ro;
};

typedef union Operand OPERAND;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

INT NS_DIM_PREFIX cmdintbufsize=CMDINTBUFSIZE;


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static int doneFlag;
#ifdef ModelP
static int interactiveFlag=0;
#endif

/* variables for command interpreter */
static const char *cmdPtr,*cmdStart;
static long executePos=0;
static char *cmdBuffer;
static char *executeBuffer;
static INT programFlag=FALSE;
static long programbufsize=PROGRAMBUFSIZE;
static char *programBuffer;
static char fileBuffer[FILEBUFSIZE+1];
static char stringBuffer[MAXSTRINGSIZE];

/* for output during execute */
static int mutelevel=0;

#ifdef ModelP
static char *execCmdBuffer;
#endif

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static INT InterpretString(void);

/* are called by GetFactor, GetAnItem */
static INT GetEquation(OPERAND *result);
static INT GetCondition(DOUBLE *result);

/****************************************************************************/
/** \brief Skip "blanks" in a command string

   This function moves the command pointer `cmdPtr` until it encounters the first
   non-blank. Here `blank` means space, tab and `\n`-character
   as well as comments of the form '#....EOL'.

   \return The final non-blank character.
 */
/****************************************************************************/

static char SkipBlanks (void)
{
  char c;

  while ((c=*cmdPtr)!=(char) 0)
  {
    if (c=='#')
    {
      /* comment: skip end of line */
      while ((c=*(++cmdPtr))!='\n')
        if (c==(char) 0) return(c);
    }
    else
    {
      if (strchr(blanks,(int) c)!=NULL)
        cmdPtr++;
      else
        break;
    }
  }

  return(c);
}

/****************************************************************************/
/** \brief Skip all text up to the next seperator

   Skips the program text until it encounters the first seperator (i.e. one of ";{}").

   \return The final separator.
 */
/****************************************************************************/

static char SkipAll (void)
{
  char c;

  while ((c=*cmdPtr)!=(char) 0)
  {
    if (c=='#')
    {
      /* comment: skip end of line */
      while (c!='\n')
      {
        c=*(++cmdPtr);
        if (c==(char) 0) return(c);
      }
      continue;
    }

    if (c=='"')
    {
      /* string: skip to next " */
      while ((c=*(++cmdPtr))!='"')
        if (c==(char) 0) return(c);
      cmdPtr++;
      continue;
    }

    if (strchr(seperators,(int) c)==NULL)
      cmdPtr++;
    else
      break;
  }

  return(c);
}

/****************************************************************************/
/** \brief Reads an item (number, token) into a buffer

   \param itemType - address where type of item is written to
   \param buffer - buffer where item is written to

   This function skips blanks and reads the following number or token
   in the memory pointed to by 'buffer'. It understands also names
   like '..:variable' and computes indices of the form [`expression`].

   This function checks if the length of the item is larger than
   the constant 'MAXTOKENLENGTH'. If it encounters that error it
   prints an error message and terminates. So take care that your buffer
   has at least that size!

   \sa
   SkipBlanks, GetEquation, PrintErrorMessage

   \return Error code, if not equal zero.
 */
/****************************************************************************/

static INT GetAnItem (int *itemType,char *buffer)
{
  INT error;
  char c,localBuffer[MAXTOKENLENGTH];
  int k,l;
  OPERAND result;

  c=SkipBlanks();

  if (ISNUMBER(c)||((c=='.') && (cmdPtr[1]!='.')))
  {
    k=0;
    if (c!='.')
    {
      /* get number before decimal point */
      while (ISNUMBER(c))
      {
        if (k==MAXTOKENLENGTH-1)
        {
          PrintErrorMessage('E',"GetAnItem","token too long");
          return(8400);                                 /* token too long */
        }
        buffer[k++]=c;
        c=*(++cmdPtr);
      }

    }

    if (c=='.')
    {
      if (k==MAXTOKENLENGTH-1)
      {
        PrintErrorMessage('E',"GetAnItem","token too long");
        return(8400);                           /* token too long */
      }
      buffer[k++]=c;
      c=*(++cmdPtr);

      /* get number after decimal point */
      while (ISNUMBER(c))
      {
        if (k==MAXTOKENLENGTH-1)
        {
          PrintErrorMessage('E',"GetAnItem","token too long");
          return(8400);                                 /* token too long */
        }
        buffer[k++]=c;
        c=*(++cmdPtr);
      }
    }

    if ((c=='e')||(c=='E'))
    {
      if (k==MAXTOKENLENGTH-1)
      {
        PrintErrorMessage('E',"GetAnItem","token too long");
        return(8400);                           /* token too long */
      }
      buffer[k++]=c;
      c=*(++cmdPtr);

      if ((c=='-')||(c=='+'))
      {
        if (k==MAXTOKENLENGTH-1)
        {
          PrintErrorMessage('E',"GetAnItem","token too long");
          return(8400);                                 /* token too long */
        }
        buffer[k++]=c;
        c=*(++cmdPtr);
      }

      /* get exponent */
      while (ISNUMBER(c))
      {
        if (k==MAXTOKENLENGTH-1)
        {
          PrintErrorMessage('E',"GetAnItem","token too long");
          return(8400);                                 /* token too long */
        }
        buffer[k++]=c;
        c=*(++cmdPtr);
      }
    }

    buffer[k]=(char) 0;
    *itemType=NUMBERID;
    return(DONE);
  }

  if (ISALPHA(c)||((c=='.')&&(cmdPtr[1]=='.')))
  {
    k=0;
    do
    {
      if (k==MAXTOKENLENGTH-1)
      {
        PrintErrorMessage('E',"GetAnItem","token too long");
        return(8400);                           /* token too long */
      }
      buffer[k++]=c;
      c=*(++cmdPtr);
      if (c=='[')
      {
        cmdPtr++;
        if ((error=GetEquation(&result))!=DONE)
          return(error);

        switch (result.ro.type)
        {
        case NUMBERID :
          sprintf(localBuffer,"%-.14g",result.ro.value);
          l=strlen(localBuffer);
          if (k+l>=MAXTOKENLENGTH-1)
          {
            PrintErrorMessage('E',"GetAnItem","token too long");
            return(8400);                                       /* token too long */
          }
          strcpy(buffer+k,localBuffer);
          k+=l;
          break;

        case ALPHAID :
          l=strlen(result.so.sptr);
          if (k+l>=MAXTOKENLENGTH-1)
          {
            PrintErrorMessage('E',"GetAnItem","token too long");
            return(8400);                                       /* token too long */
          }
          strcpy(buffer+k,result.so.sptr);
          k+=l;
          break;

        case LSTRINGID :
          l=result.lo.length;
          if (k+l>=MAXTOKENLENGTH-1)
          {
            PrintErrorMessage('E',"GetAnItem","token too long");
            return(8400);                                       /* token too long */
          }
          strncpy(buffer+k,result.lo.sptr,(size_t) l);
          k+=l;
          break;
        }

        c=SkipBlanks();

        if (c!=']')
        {
          PrintErrorMessage('E',"GetAnItem","index does not terminate with ]");
          return(PARAMERRORCODE);
        }

        c=*(++cmdPtr);
      }
    }
    while (ISALPHA(c)||ISNUMBER(c)||(c=='.'));

    buffer[k]=(char) 0;
    *itemType=ALPHAID;
    return(DONE);
  }

  buffer[0]=(char) 0;
  *itemType=EMPTYID;
  return(DONE);
}

/****************************************************************************/
/** \brief Converts a string to a double

   \param strAdr - address of string or lstring
   \param maxLen - 0 if normal string, length of string if lstring
   \param endAdr - memory where the end address is stored
   \param type - address for type of string
   \param valuePtr - address where value is store

   This function converts the string at location 'strAdr' to a number if it is
   possible. It skips blanks at the beginning (and end), and tries
   to locate a floating point number at the beginning. If
   this is not possible '*type' is set to the constant 'ALPHAID' and
   '*valuePtr' to zero.

   An error is reported if the number string is longer than 'MAXTOKENLENGTH'.

   SEE ALSO:
   SkipBlanks, PrintErrorMessage

   RETURN VALUE:
   INT
   .n Error code, if not equal zero
 */
/****************************************************************************/

static INT ConvertStringToDouble(const char *strAdr, int maxLen, const char **endAdr, INT *type, DOUBLE *valuePtr)
{
  char c;
  int endPos,anfPos,k;
  DOUBLE sign;

  endPos=maxLen;
  if (maxLen==0)
    endPos=strlen(strAdr);

  /* skip trailing blanks */
  while (endPos>=0)
    if (strchr(blanks,(int) (c=strAdr[--endPos]))==NULL)
      break;

  k=0;
  /* skip blanks and minus signs */
  sign=1;
  while (k<=endPos)
  {
    if ((c=strAdr[k])!=' ')
    {
      if (c!='-')
        break;
      sign=-sign;
    }
    k++;
  }
  anfPos=k;

  if (ISNUMBER(c)||(c=='.'))
  {
    if (c!='.')
    {
      /* get number before decimal point */
      while (ISNUMBER(c))
      {
        if (k>endPos)
          break;
        else
          c=strAdr[k++];
      }
    }

    if (c=='.')
    {
      if (k<=endPos)
        c=strAdr[k++];

      /* get number after decimal point */
      while (ISNUMBER(c))
      {
        if (k>endPos)
          break;
        else
          c=strAdr[k++];
      }
    }

    if ((c=='e')||(c=='E'))
    {
      if (k<=endPos)
        c=strAdr[k++];

      if ((c=='-')||(c=='+'))
        if (k<=endPos)
          c=strAdr[k++];

      /* get exponent */
      while (ISNUMBER(c))
      {
        if (k>endPos)
          break;
        else
          c=strAdr[k++];
      }
    }
  }

  if (ISNUMBER(c))
  {
    *type=NUMBERID;
    if (valuePtr!=NULL)
    {
      if (endPos-anfPos+2<=MAXTOKENLENGTH)
      {
        strncpy(stringBuffer,strAdr+anfPos,endPos-anfPos+1);
        stringBuffer[endPos-anfPos+1]=(char) 0;
        *valuePtr=sign*atof(stringBuffer);
      }
      else
      {
        *valuePtr=0;
        PrintErrorMessage('E',"ConvertStringToDouble","number too long");
        return(8405);                           /* number too long */
      }
    }
  }
  else
  {
    *type=ALPHAID;
    if (valuePtr!=NULL)
      *valuePtr=0;
  }

  if (endAdr!=NULL)
    *endAdr=strAdr;

  return(DONE);
}

/****************************************************************************/
/** \brief Evaluates an OPERAND

   \param t - location where the value is put
   \param term - the operand to be evaluated

   This routine sets *t to the arithmetic value of the given operand.

   An error is reported if the operand is a string which does not contain a number.

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT GetValueOfOperand (DOUBLE *t, OPERAND *term)
{
  INT itemType,error;

  switch (term->ro.type)
  {
  case NUMBERID :
    *t=term->ro.value;
    error=FALSE;
    break;

  case ALPHAID :
    if ((error=ConvertStringToDouble(term->so.sptr,0,NULL,&itemType,t))!=DONE)
      return(error);
    error=(itemType!=NUMBERID);
    break;

  case LSTRINGID :
    if ((error=ConvertStringToDouble(term->lo.sptr,term->lo.length,NULL,&itemType,t))!=DONE)
      return(error);
    error=(itemType!=NUMBERID);
    break;

  default :
    error=TRUE;
  }

  if (error)
  {
    PrintErrorMessage('E',"GetValueOfOperand","wrong item type");
    return(8606);
  }
  else
    return(DONE);
}

/****************************************************************************/
/** \brief Compares two strings (normal or lstring)

   \param Diff - at this location the difference `val(string1)-val(string2)`
   is stored if both strings are numbers (if not the result of a string comparison).
   \param str1 - address of first string
   \param str2 - address of second string
   \param maxLen1 - should be 0 if str1 is a normal string and its length if it is an lstring
   \param maxLen2 - should be 0 if str2 is a normal string and its length if it is an lstring

   This routine compares two strings. If both can be evaluated to floating point
   numbers, their values are compared. The result is placed at position '*Diff'.

   SEE ALSO:
   ConvertStringToDouble

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT StringCompare(DOUBLE *Diff, const char *str1, const char *str2, int maxLen1, int maxLen2)
{
  DOUBLE value1,value2;
  int maxLen;
  INT type1,type2,error;

  *Diff=0;

  if ((error=ConvertStringToDouble(str1,maxLen1,NULL,&type1,&value1))!=DONE)
    return(error);

  if ((error=ConvertStringToDouble(str2,maxLen2,NULL,&type2,&value2))!=DONE)
    return(error);


  if ((type1==NUMBERID) && (type2==NUMBERID))
    *Diff=value1-value2;
  else
  {
    maxLen=(maxLen1>maxLen2) ? maxLen1 : maxLen2;
    if (maxLen<=0)
      *Diff=(DOUBLE) strcmp(str1,str2);
    else
      *Diff=(DOUBLE) strncmp(str1,str2,maxLen);
  }

  return(DONE);
}

/****************************************************************************/
/** \brief Reads a token from cmdPtr into a buffer

   \param buffer - the buffer for the token

   This command reads a variable name or a command into a buffer.
   and returns an error if it encounters something different. The buffer
   should have at least length 'MAXTOKENLENGTH'.

   SEE ALSO:
   GetAnItem

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT GetToken (char *buffer)
{
  int error,itemType;

  if ((error = GetAnItem(&itemType,buffer))!=DONE)
    return(error);

  if (itemType!=ALPHAID)
    return(1);
  else
    return(DONE);
}

/****************************************************************************/
/** \brief Evaluates a given term a, -a, @a, "...", `simple_function`(a) or (`expression`)

   \param result - pointer where the result is to be stored

   This function reads in the term addressed by the global pointer 'cmdPtr'.
   It handles (multiple) minus signs,
   substitution of variable values (including modification of variable names via the "@"-operator),
   strings and simple functions like `sqrt`, `exp`, `log`, `sin`, `cos`, `floor`.
   Last not least, it contains a call to the function 'EvaluateExpression'
   if an opening bracket "(" is encountered.

   SEE ALSO:
   SkipBlanks, EvaluateExpression, GetToken, GetStringVar

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/


static INT GetFactor (OPERAND *result)
{
  char c,c1;
  char buffer[MAXTOKENLENGTH];
  const char *stringAdr;
  const char *savedCmdPtr;
  int k;
  INT error;
  int itemType;
  int signflag;
  DOUBLE sign,t;
  OPERAND newResult;
  const char *lastname;
  ENVDIR *theDir;

  /* set error value */
  result->ro.type=EMPTYID;

  c=SkipBlanks();

  /* get signs */
  sign=1;
  signflag=FALSE;

  while (c=='-')
  {
    signflag=TRUE;
    sign=-sign;
    cmdPtr++;
    c=SkipBlanks();
  }

  /* get expression */
  switch (c)
  {
  case '"' :
  case '\'' :
    stringAdr=++cmdPtr;
    k=0;
    while ((c1=*cmdPtr)!=c)
    {
      if (c1==(char) 0)
      {
        PrintErrorMessage('E',"GetFactor","eof while reading string");
        return(8600);                           /* eof while reading string */
      }
      else
      {
        cmdPtr++;
        k++;
      }
    }
    cmdPtr++;                   /* skip terminating " */

    newResult.lo.type=LSTRINGID;
    newResult.lo.sptr=stringAdr;
    newResult.lo.length=k;
    break;

  case '@' :
    cmdPtr++;
    if ((error=GetToken(buffer))!=DONE)
    {
      PrintErrorMessage('E',"InterpretString","syntax error");
      return(error);                    /* syntax error */
    }

    if ((stringAdr=GetStringVar(buffer))!=NULL)
    {
      /* evaluate expression given in the variable */
      savedCmdPtr=cmdPtr;
      cmdPtr=stringAdr;
      if ((error=GetEquation(&newResult))!=DONE)
      {
        cmdPtr=savedCmdPtr;
        return(error);
      }

      if ((c=SkipBlanks())!=(char) 0)
      {
        PrintErrorMessage('E',"GetFactor","syntax error");
        return(8403);                           /* syntax error */
      }

      cmdPtr=savedCmdPtr;
    }
    break;

  case '(' :
    cmdPtr++;
    if ((error=GetEquation(&newResult))!=DONE)
      return(error);

    c=SkipBlanks();
    if (c!=')')
    {
      PrintErrorMessage('E',"GetFactor","syntax error");
      return(8403);                     /* syntax error */
    }

    cmdPtr++;                   /* skip terminating ) */
    break;

  default :
    if ((error=GetAnItem(&itemType,buffer))!=DONE)
      return(error);

    switch (itemType)
    {
    case NUMBERID :
      newResult.ro.type=NUMBERID;
      newResult.ro.value=atof(buffer);
      break;

    case ALPHAID :
      if (strcmp(buffer,"ugCmd")==0)
      {
        break;                                  /* will be implemented rsn */
      }
      /* is it a built in function? */
      if (    (strcmp(buffer,"exp")==0)
              ||      (strcmp(buffer,"log")==0)
              ||      (strcmp(buffer,"fabs")==0)
              ||      (strcmp(buffer,"floor")==0)
              ||      (strcmp(buffer,"sin")==0)
              ||      (strcmp(buffer,"cos")==0)
              ||      (strcmp(buffer,"sqrt")==0)
              )
      {
        if ((error=GetCondition(&t))!=DONE)
          return(error);

        newResult.ro.type=NUMBERID;

        if (strcmp(buffer,"exp")==0)
          newResult.ro.value=exp(t);

        if (strcmp(buffer,"log")==0)
          newResult.ro.value=log(t);

        if (strcmp(buffer,"fabs")==0)
          newResult.ro.value=fabs(t);

        if (strcmp(buffer,"floor")==0)
          newResult.ro.value=floor(t);

        if (strcmp(buffer,"sin")==0)
          newResult.ro.value=sin(t);

        if (strcmp(buffer,"cos")==0)
          newResult.ro.value=cos(t);

        if (strcmp(buffer,"sqrt")==0)
          newResult.ro.value=sqrt(t);
      }
      else if (strcmp(buffer,"def")==0)
      {
        /* check wether a string variable or a struct is defined */

        INT deref=FALSE;
        char *p;

        c=SkipBlanks();
        if (c!='(')
        {
          PrintErrorMessage('E',"def","'(' missing");
          return(8604);                                         /* ( missing */
        }
        cmdPtr++;

        c=SkipBlanks();
        if (c=='@')
        {
          deref = TRUE;
          cmdPtr++;
        }
        GetToken(buffer);
        p = buffer;
        if (deref)
          if ((p=GetStringVar(buffer))==NULL)
          {
            PrintErrorMessageF('E',"def","string var '%s' not found",buffer);
            return(__LINE__);                                                   /* ( missing */
          }

        /* check for string variable */
        newResult.ro.type  = NUMBERID;
        newResult.ro.value = (GetStringVar(p)!=NULL);

        if (newResult.ro.value==0)
          /* check for struct */
          if (((theDir=FindStructDir(buffer,&lastname))!=NULL)
              && (FindStructure(theDir,lastname)!=NULL))
            newResult.ro.value = 2;


        c=SkipBlanks();
        if (c!=')')
        {
          PrintErrorMessage('E',"def","')' missing");
          return(8604);                                         /* ) missing */
        }
        cmdPtr++;
      }
      else
      {
        /* get ug variable */

        if ((stringAdr=GetStringVar(buffer))==NULL)
        {
          char buf[64];
          sprintf(buf,"variable %s not found",buffer);
          PrintErrorMessage('E',"GetFactor",buf);
          return(8601);                                         /* variable not found */
        }

        newResult.so.type=ALPHAID;
        newResult.so.sptr=stringAdr;
        break;
      }

    default :
      break;
    }
    break;

  }             /* endswitch c */

  /* copy newResult to result and convert string to number if there is a minus sign in front */
  switch (result->ro.type=newResult.ro.type)
  {
  case NUMBERID :
    result->ro.value=sign*newResult.ro.value;
    break;

  case ALPHAID :
  case LSTRINGID :
    if (signflag)
    {
      result->ro.type=NUMBERID;
      result->ro.value=sign*atof(newResult.so.sptr);                    /* here might arise problems: 123".... */
    }
    else
    {
      result->lo.sptr=newResult.lo.sptr;
      result->lo.length=newResult.lo.length;
    }
    break;
  case EMPTYID :
    if (signflag)
    {
      PrintErrorMessage('E',"GetFactor","syntax error");
      return(8602);                     /* syntax error */
    }
    break;
  }

  return(DONE);
}

/****************************************************************************/
/** \brief  Evaluates a given term of the form a*b/c*d...

   \param result - memory address for the result

   This function evaluates a product (consisting of factors) in the program text.

   SEE ALSO:
   SkipBlanks, GetFactor

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/


static INT GetProduct(OPERAND *result)
{
  char c;
  INT error;
  OPERAND newResult;

  result->ro.type=EMPTYID;

  if ((error=GetFactor(&newResult))!=DONE)
    return(error);

  if (newResult.ro.type==EMPTYID)
    return(DONE);

  c=SkipBlanks();

  if ((c=='*') || (c=='/') || (c=='%'))
  {
    result->ro.type=NUMBERID;
    switch (newResult.ro.type)
    {
    case NUMBERID :
      result->ro.value=newResult.ro.value;
      break;

    case ALPHAID :
    case LSTRINGID :
      result->ro.value=atof(newResult.so.sptr);
      break;
    }
  }
  else
  {
    switch (result->ro.type=newResult.ro.type)
    {
    case NUMBERID :
      result->ro.value=newResult.ro.value;
      break;

    case ALPHAID :
      result->so.sptr=newResult.so.sptr;
      break;

    case LSTRINGID :
      result->lo.sptr=newResult.lo.sptr;
      result->lo.length=newResult.lo.length;
      break;
    }
    return(DONE);
  }

  do {

    cmdPtr++;

    if ((error=GetFactor(&newResult))!=DONE)
      return(error);

    if (newResult.ro.type==EMPTYID)
    {
      PrintErrorMessage('E',"GetProduct","incomplete operation");
      return(1);
    }

    switch (c)
    {
    case '*' :
      switch (newResult.ro.type)
      {
      case NUMBERID :
        result->ro.value*=newResult.ro.value;
        break;

      case ALPHAID :
      case LSTRINGID :
        result->ro.value*=atof(newResult.so.sptr);
        break;
      }
      break;

    case '/' :
      switch (newResult.ro.type)
      {
      case NUMBERID :
        result->ro.value/=newResult.ro.value;
        break;

      case ALPHAID :
      case LSTRINGID :
        result->ro.value/=atof(newResult.so.sptr);
        break;
      }
      break;

    case '%' :
      switch (newResult.ro.type)
      {
      case NUMBERID :
        result->ro.value = (DOUBLE)(lrint(result->ro.value)%lrint(newResult.ro.value));
        break;

      case ALPHAID :
      case LSTRINGID :
        result->ro.value = (DOUBLE)(lrint(result->ro.value)%lrint(atof(newResult.so.sptr)));
        break;
      }
      break;

    }                   /* switch */

    c=SkipBlanks();

  } while ((c=='*') || (c=='/') || (c=='%'));

  return(DONE);

}       /* GetProduct */


/****************************************************************************/
/** \brief Evaluates a given term of the form a+b-c+d...

   \param result - memory address for the result

   This function evaluates a sum (consisting of products) in the program text.

   SEE ALSO:
   SkipBlanks, GetProduct

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/


static INT GetSum(OPERAND *result)
{
  char c;
  INT error;
  OPERAND newResult;

  result->ro.type=NUMBERID;

  if ((error=GetProduct(&newResult))!=DONE)
    return(error);

  if (newResult.ro.type==EMPTYID)
    return(DONE);

  c=SkipBlanks();

  if ((c=='+') || (c=='-'))
  {
    result->ro.type=NUMBERID;
    switch (newResult.ro.type)
    {
    case NUMBERID :
      result->ro.value=newResult.ro.value;
      break;

    case ALPHAID :
    case LSTRINGID :
      result->ro.value=atof(newResult.so.sptr);
      break;
    }
  }
  else
  {
    switch (result->ro.type=newResult.ro.type)
    {
    case NUMBERID :
      result->ro.value=newResult.ro.value;
      break;

    case ALPHAID :
      result->so.sptr=newResult.so.sptr;
      break;

    case LSTRINGID :
      result->lo.sptr=newResult.lo.sptr;
      result->lo.length=newResult.lo.length;
      break;
    }
    return(DONE);
  }

  do {

    cmdPtr++;

    if ((error=GetProduct(&newResult))!=DONE)
      return(error);

    if (newResult.ro.type==EMPTYID)
    {
      PrintErrorMessage('E',"GetProduct","incomplete operation");
      return(1);
    }

    switch (c)
    {
    case '+' :
      switch (newResult.ro.type)
      {
      case NUMBERID :
        result->ro.value+=newResult.ro.value;
        break;

      case ALPHAID :
      case LSTRINGID :
        result->ro.value+=atof(newResult.so.sptr);
        break;
      }
      break;

    case '-' :
      switch (newResult.ro.type)
      {
      case NUMBERID :
        result->ro.value-=newResult.ro.value;
        break;

      case ALPHAID :
      case LSTRINGID :
        result->ro.value-=atof(newResult.so.sptr);
        break;
      }
      break;

    }                   /* switch */

    c=SkipBlanks();

  } while ((c=='+') || (c=='-'));

  return(DONE);

}       /* GetSum */


/****************************************************************************/
/** \brief Evaluates a term of the form a or a `rel` b with `rel` in ==, !=, <, >

   \param result - memory address for the result

   This function evaluates a relation (between sums) in the program text.
   Also string comparisons are handled.

   SEE ALSO:
   SkipBlanks, GetSum

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT GetEquation (OPERAND *result)
{
  char c;
  INT itemType,error;
  DOUBLE t;
  OPERAND term1,term2;

  result->ro.type=NUMBERID;

  if ((error=GetSum(&term1))!=DONE)
    return(error);

  if (term1.ro.type==EMPTYID)
    return(DONE);

  switch (c=SkipBlanks())
  {
  case '=' :
  case '!' :
    if ((*(++cmdPtr))!='=')
    {
      PrintErrorMessage('E',"GetEquation","syntax error");
      return(8603);                     /* syntax error */
    }
    cmdPtr++;
    if ((error=GetSum(&term2))!=DONE)
      return(error);
    break;

  case '<' :
  case '>' :
    if ((*(++cmdPtr))=='=')
    {
      if (c=='<')
        c='k';
      else
        c='g';
      cmdPtr++;
    }
    if ((error=GetSum(&term2))!=DONE)
      return(error);
    break;

  default :
    /* no equation type construct -> return term1 */
    switch (result->ro.type=term1.ro.type)
    {
    case NUMBERID :
      result->ro.value=term1.ro.value;
      break;

    case ALPHAID :
      result->so.sptr=term1.so.sptr;
      break;

    case LSTRINGID :
      result->lo.sptr=term1.lo.sptr;
      result->lo.length=term1.lo.length;
      break;
    }
    return(DONE);

  }             /* switch */

  if (term2.ro.type==EMPTYID)
  {
    PrintErrorMessage('E',"GetProduct","incomplete operation");
    return(1);
  }

  switch (term1.ro.type)
  {
  case NUMBERID :
    if ((error=GetValueOfOperand(&t,&term2))!=DONE)
      return(error);
    t=term1.ro.value-t;
    break;

  case ALPHAID :
    switch (term2.ro.type)
    {
    case NUMBERID :
      if ((error=ConvertStringToDouble(term1.so.sptr,0,NULL,&itemType,&t))!=DONE)
        return(error);
      if (itemType!=NUMBERID)
      {
        PrintErrorMessage('E',"GetEquation","wrong item id");
        return(8606);
      }

      t=t-term2.ro.value;
      break;

    case ALPHAID :
      if ((error=StringCompare(&t,term1.so.sptr,term2.so.sptr,0,0))!=DONE)
        return(error);
      break;

    case LSTRINGID :
      if ((error=StringCompare(&t,term1.so.sptr,term2.lo.sptr,strlen(term1.so.sptr),(size_t) term2.lo.length))!=DONE)
        return(error);
      break;
    }
    break;

  case LSTRINGID :
    switch (term2.ro.type)
    {
    case NUMBERID :
      if ((error=ConvertStringToDouble(term1.lo.sptr,term1.lo.length,NULL,&itemType,&t))!=DONE)
        return(error);
      if (itemType!=NUMBERID)
      {
        PrintErrorMessage('E',"GetEquation","wrong item id (number expected)");
        return(8606);
      }
      t=t-term2.ro.value;
      break;

    case ALPHAID :
      if ((error=StringCompare(&t,term1.lo.sptr,term2.so.sptr,(size_t) term1.lo.length,strlen(term2.so.sptr)))!=DONE)
        return(error);
      break;

    case LSTRINGID :
      if ((error=StringCompare(&t,term1.lo.sptr,term2.lo.sptr,(size_t) term1.lo.length,(size_t) term2.lo.length))!=DONE)
        return(error);
      break;
    }
    break;
  }

  result->ro.type=NUMBERID;
  switch (c)
  {
  case '=' :
    if (t==0)
      result->ro.value=1;
    else
      result->ro.value=0;
    break;

  case '!' :
    if (t!=0)
      result->ro.value=1;
    else
      result->ro.value=0;
    break;

  case '<' :
    if (t<0)
      result->ro.value=1;
    else
      result->ro.value=0;
    break;

  case '>' :
    if (t>0)
      result->ro.value=1;
    else
      result->ro.value=0;
    break;

  case 'k' :
    if (t<=0)
      result->ro.value=1;
    else
      result->ro.value=0;
    break;

  case 'g' :
    if (t>=0)
      result->ro.value=1;
    else
      result->ro.value=0;
    break;
  }

  return(DONE);

}

/****************************************************************************/
/** \brief Evaluates an condition enclosed in ()

   \param result - memory address for the result

   This function checks for enclosing brackets "(...)" and calls
   'GetEquation' for the term in between.
   The result of the "equation" is converted to a double and stored in '*result'.

   SEE ALSO:
   SkipBlanks, GetEquation, GetValueOfOperand

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT GetCondition (DOUBLE *result)
{
  char c;
  INT error;
  OPERAND newResult;

  c=SkipBlanks();

  if (c!='(')
  {
    PrintErrorMessage('E',"GetCondition","'(' missing");
    return(8604);               /* ( missing */
  }
  cmdPtr++;
  if ((error=GetEquation(&newResult))!=DONE)
    return(error);

  c=SkipBlanks();

  if (c!=')')
  {
    PrintErrorMessage('E',"GetCondition","')' missing");
    return(8604);               /* ) missing */
  }

  cmdPtr++;

  if ((error=GetValueOfOperand(result,&newResult))!=DONE)
    return(error);

  return(DONE);
}

/****************************************************************************/
/** \brief Evaluates an expression

   \param result - memory address for the result

   This function is a call of 'SkipBlanks' and 'GetEquation'.

   SEE ALSO:
   SkipBlanks, GetEquation

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT EvaluateExpression(OPERAND *result)
{
  char c;
  INT error;

  c=SkipBlanks();

  if ((error=GetEquation(result))!=DONE)
    return(error);

  return(DONE);
}

/****************************************************************************/
/** \brief Reads a full path name for the execute command

   \param buffer - memory address for the name

   This function reads a long path name including "/~."-signs from the program text
   into a buffer. This buffer should have at least length 'MAXTOKENLENGTH'.

   SEE ALSO:
   SkipBlanks

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT GetFullPathName (char *buffer)
{
  char c;
  char svarname[NAMESIZE],*sptr;
  int k,error;

  c=SkipBlanks();

  if (c=='@')
  {
    cmdPtr++;
    if ((error=GetToken(svarname))!=DONE)
    {
      PrintErrorMessage('E',"InterpretString","syntax error");
      return(error);                    /* syntax error */
    }

    if ((sptr=GetStringVar(svarname))==NULL)
    {
      char buf[64];
      sprintf(buf,"variable %s not found",svarname);
      PrintErrorMessage('E',"InterpretString",buf);
      return(8522);                     /* variable not found */
    }
    strcpy(buffer,sptr);
    return (DONE);
  }
  else if (ISALPHA(c)||(c=='/')||(c=='~')||((c=='.')&&(cmdPtr[1]=='.'))||((c=='.')&&(cmdPtr[1]=='/')))
  {
    k=0;
    do
    {
      if (k==MAXTOKENLENGTH-1)
      {
        PrintErrorMessage('E',"GetFullPathName","token too long");
        return(8402);                           /* token too long */
      }
      buffer[k++]=c;
      c=*(++cmdPtr);
    }
    while (ISALPHA(c)||ISNUMBER(c)||(c=='.')||(c=='/'));

    buffer[k]=(char) 0;

    return(DONE);
  }
  else
  {
    PrintErrorMessage('E',"GetFullPathName","invalid path");
    return(8403);               /* invalid path */
  }
}

FILE* NS_DIM_PREFIX FOpenScript (const char *script, const char *mode)
{
  if (scriptpaths_set)
    return FileOpenUsingSearchPaths(script,mode,"scriptpaths");
  else
    return fileopen(script,mode);
}

/****************************************************************************/
/** \brief Interprets a sequence of commands terminated by (char) 0.

   \param cmds - start of program string

   This command is a call to the routine 'InterpretString'. This subdivision has the purpose
   of making recursive calls possible by buffering the global variable 'cmdPtr'.

   It also implements a very crude programming facility via the commands
   'program' and 'endprogram'.

   SEE ALSO:
   InterpretString, GetMuteLevel

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

INT NS_DIM_PREFIX InterpretCommand (const char *cmds)
{
  int pLength;
  INT error;
  const char *oldCmdPtr,*oldCmdStart;

  mutelevel = GetMuteLevel();

  if ((strcmp(cmds,"program")==0)||(strcmp(cmds,"program\n")==0))
  {
    programFlag=TRUE;
    programBuffer[0]=0;
    return(DONE);
  }

  if ((strcmp(cmds,"endprogram")==0)||(strcmp(cmds,"endprogram\n")==0))
  {
    programFlag=FALSE;
    cmds=programBuffer;
  }

  if (programFlag==TRUE)
  {
    pLength=strlen(programBuffer);
    if (pLength+strlen(cmds)+1>programbufsize-1)
    {
      programBuffer[0]=(char) 0;
      programFlag=FALSE;

      PrintErrorMessage('E',"InterpretCommand","unexpected end");
      return(8512);
    }

    programBuffer[pLength]=(char) 13;
    programBuffer[pLength+1]=(char) 0;
    strcat(programBuffer,cmds);

    return(DONE);
  }

  /* save cmdPtr/cmdStart */
  oldCmdPtr=cmdPtr;
  oldCmdStart=cmdStart;

  /* Initialize cmdPtr */
  cmdPtr=cmdStart=cmds;

  /* call InterpretString */
  if ((error = InterpretString())!=DONE)
  {
    SetMuteLevel(0);
    return(error);
  }

  /* restore cmdPtr/cmdStart and return */
  cmdPtr=oldCmdPtr;
  cmdStart=oldCmdStart;
  return(DONE);

}       /* InterpretCommand */

/****************************************************************************/
/** \brief Interprets a string at address '*cmdPtr'.

   This function is the heart of the interpreter. It moves the program counter
   given by the global variable 'cmdPtr' around in the program text, thereby
   executing the program.

   It handles the `structure` commands 'if', 'then', 'else', 'repeat', 'break', 'continue'
   as well as a very simple input via 'input' and output via 'print'. A substitution
   operator '@' allows a kind of subroutines. Also the important
   'execute' (short 'ex') of a script is handled here.

   The interpreter uses only text variables and diretories of those, called `structures`.
   Of course, these variables may contain numbers, and with these the interpreter
   can perform simple computations ('+', '-', '*', '/'),
   comparisons ('==', '!=', '<', '>'), and some more complicated mathematical functions
   like 'sqrt', 'exp', 'log', 'sin', 'cos', 'floor'. Indirect referencing
   is possible by using the above '@'-operator or indexing via '[]'.

   If some line can not be handled by the interpreter, it is handed over to `ug`
   to be processed there. Variable passing is possible by using the
   '@'-operator or via the variables which can be accessed by `ug`.

   SEE ALSO:
   InterpretCommand, EvaluateCondition, EvaluateExpression, GetToken, GetStringVar, ExecCommand
   and of course the text on the interpreter.

   RETURN VALUE:
   INT
   .n Error code, if not zero.
 */
/****************************************************************************/

static INT InterpretString (void)
{
  Boolean termFlag;
  char c,c1;
  int k,value;
  INT error;
  int status,status0;
  long executePos0;

  DOUBLE res;

  int StatusStack[BUFSIZE];
  int RepeatStatusPos[MAXREPEAT];
  const char *RepeatPtr[MAXREPEAT];
  char buffer[MAXTOKENLENGTH];
  char filename[MAXTOKENLENGTH];
  char valueStr[32];
  char *sptr;

  int StatusPos,RepeatPos;
  const char *cmdStart;
  const char *cmdEnd;

  COMMAND *commandItem;

  OPERAND result;

  FILE *filePtr;


  status=MULTIPLE;
  StatusPos=0;
  RepeatPos=0;
  termFlag=FALSE;

  do
  {
    mutelevel = GetMuteLevel();

    if (UserInterrupt("InterpretString"))
      return(8512);                     /* UserInterrupt */

    /* if in terminating mode terminate until in multiple mode */
    if (termFlag)
    {
      if (status&MULTIPLE)
      {
        termFlag=FALSE;                                 /* leave terminating mode */
        continue;
      }

      if ((status&STATUSCODE)==IF)
      {
        cmdStart=cmdPtr;                                        /* search for else */

        if (GetToken(buffer)==DONE)
          if (strcmp(buffer,"else")==0)
          {
            status=status-IF+ELSE;
            status0=StatusStack[StatusPos-1];
            if ((status0&SKIPMODE)==0)
              status=(status^SKIPMODE);                                                         /* if-else active: switch to opposite */
            termFlag=FALSE;                                             /* leave terminating mode */
            continue;
          }

        cmdPtr=cmdStart;                                /* set cmdPtr back */
        status=StatusStack[--StatusPos];
        continue;
      }

      if ((status&STATUSCODE)==ELSE)
      {
        status=StatusStack[--StatusPos];
        continue;
      }

      if ((status&STATUSCODE)==REPEAT)
      {
        if ((status&SKIPMODE)==0)
        {
          if (RepeatPos==0)
          {
            PrintErrorMessage('E',"InterpretString","???");
            return(8506);                                       /* should not happen! */
          }
          else
          {
            cmdPtr=RepeatPtr[RepeatPos-1];                                              /* status is already ok */
            termFlag=FALSE;
            continue;
          }
        }
        else
        {
          status=StatusStack[--StatusPos];                                      /* leave inactive repeat loop */
          continue;
        }
      }

      PrintErrorMessage('E',"InterpretString","???");
      return(8524);                     /* should not happen: in which status are you? */
    }

    c=SkipBlanks();                     /* get next significant character */

    if (c==(char) 0)
    {
      if (StatusPos==0)
        break;
      else
      {
        PrintErrorMessage('E',"InterpretString","unexpected end");
        return(8501);                           /* unexpected end */
      }
    }

    if (c=='{')
    {
      if (StatusPos==BUFSIZE)
      {
        PrintErrorMessage('E',"InterpretString","program stack overflow");
        return(8502);                           /* program stack overflow */
      }
      else
      {
        StatusStack[StatusPos++]=status;
        status=((status&SKIPMODE)|MULTIPLE);
        cmdPtr++;
        continue;
      }
    }

    if (c==';')
    {
      if (status&MULTIPLE)
      {
        cmdPtr++;
        continue;
      }
      else
      {
        cmdPtr++;
        termFlag=TRUE;                                  /* in single mode ; has terminating effect */
        continue;
      }
    }

    if (c=='}')
    {
      if (status&MULTIPLE)
      {
        if (StatusPos==0)
        {
          PrintErrorMessage('E',"InterpretString","too many closing brackets");
          return(8507);                                 /* too many closing brackets */
        }
        else
        {
          cmdPtr++;
          status=StatusStack[--StatusPos];
          termFlag=TRUE;
          continue;
        }
      }
      else
      {
        PrintErrorMessage('E',"InterpretString","'}' in single mode (use ';' or '\n')");
        return(8523);                                   /* '}' in single mode (use ';' or '\n') */
      }
    }

    /* if skip mode, continue until a structure seperator ;{} is found */
    if (status&SKIPMODE)
    {
      SkipAll();
      continue;
    }

    /* check for a variable execution */
    if (c=='@')
    {
      cmdPtr++;

      if ((error=GetToken(buffer))!=DONE)
      {
        PrintErrorMessage('E',"InterpretString","syntax error");
        return(error);                          /* syntax error */
      }

      if ((cmdStart=GetStringVar(buffer))==NULL)
      {
        char buf[64];
        sprintf(buf,"variable %s not found",buffer);
        PrintErrorMessage('E',"InterpretString",buf);
        return(8522);                           /* variable not found */
      }

      /* execute string variable */
      if ((error=InterpretCommand(cmdStart))!=DONE)
        return(error);

      /* SkipAll();		skip all arguments */
      continue;
    }

    /* before checking for structure commands save command start for process command line */
    cmdStart=cmdPtr;

    if ((error=GetToken(buffer))!=DONE)
    {
      PrintErrorMessage('E',"InterpretString","syntax error");
      return(error);                    /* syntax error */
    }

    /* check for structure command */

    if (strcmp(buffer,"if")==0)
    {
      if (StatusPos==BUFSIZE)
      {
        PrintErrorMessage('E',"InterpretString","program stack overflow");
        return(8502);                           /* program stack overflow */
      }
      if ((error=GetCondition(&res))!=DONE)
        return(error);

      StatusStack[StatusPos++]=status;

      if (res!=0)
        status=IF;
      else
        status=(IF|SKIPMODE);

      continue;
    }

    if (strcmp(buffer,"repeat")==0)
    {
      if (status&SKIPMODE) continue;
      if ((StatusPos==BUFSIZE) || (RepeatPos==MAXREPEAT))
      {
        PrintErrorMessage('E',"InterpretString","stack overflow");
        return(8508);                           /* stack overflow */
      }
      StatusStack[StatusPos++]=status;
      status=REPEAT;
      RepeatStatusPos[RepeatPos]=StatusPos;
      RepeatPtr[RepeatPos++]=cmdPtr;
      continue;
    }

    if (strcmp(buffer,"break")==0)
    {
      /* All levels down to the repeat position must be skipped */
      status=(status|SKIPMODE);                                 /* skip actual level */
      if (RepeatPos>0)
        k=RepeatStatusPos[--RepeatPos];                                 /* skip from start of repeat position */
      else
        k=0;                            /* skip whole program */

      while (k<StatusPos)
      {
        StatusStack[k]=(StatusStack[k]|SKIPMODE);
        k++;
      }

      SkipAll();                        /* forget eventual arguments to break command */
      continue;
    }

    if (strcmp(buffer,"continue")==0)
    {
      /* search for break */
      if (RepeatPos==0)
      {
        PrintErrorMessage('E',"InterpretString","continue outside loop");
        return(1);
      }
      else
      {
        StatusPos=RepeatStatusPos[RepeatPos-1];
        status=REPEAT;
        cmdPtr=RepeatPtr[RepeatPos-1];
        continue;
      }
    }

    if (strcmp(buffer,"print")==0)
    {
      INT oldmute=GetMuteLevel();
      SetMuteLevel(0);

      c=(char) 0;
      do
      {
        /* skip , */
        c1=c;
        while ((c=SkipBlanks())==',')
        {c1=c; cmdPtr++;}

        if ((c==';')||(c==(char)0))
        {
          if (c1!=',')
            UserWrite("\n");
          break;
        }

        if ((error=EvaluateExpression(&result))!=DONE)
          return(error);

        switch (result.ro.type)
        {
        case NUMBERID :
          UserWriteF("%-.14g",result.ro.value);
          break;

        case ALPHAID :
          UserWrite(result.so.sptr);
          break;

        case LSTRINGID :
          if (result.lo.length>=MAXSTRINGSIZE)
          {
            PrintErrorMessage('E',"InterpretString","lstring too long");
            return(8521);
          }
          else
          {
            strncpy(stringBuffer,result.lo.sptr,(size_t) result.lo.length);
            stringBuffer[result.lo.length]=(char) 0;
            UserWrite(stringBuffer);
          }
          break;
        }

        c=SkipBlanks();
        if ((c!=',')&&(c!=';')&&(c!=(char) 0))
        {
          PrintErrorMessage('E',"InterpretString","error in print argument list");
          return(1);
        }
      }
      while (TRUE);

      SetMuteLevel(oldmute);
      continue;
    }

    if (strcmp(buffer,"input")==0)
    {
      UserRead(stringBuffer);

      error=GetToken(buffer);

      if (strcmp(buffer,"")==0)
        continue;

      if (error!=DONE)
      {
        PrintErrorMessage('E',"InterpretString","syntax error");
        return(8510);                           /* syntax error */
      }

            #ifdef ModelP
      sprintf(execCmdBuffer,"set %s = %s",buffer,stringBuffer);                   /* set cmd */
      error = ParExecCommand(execCmdBuffer);
                        #else
      error = SetStringVar(buffer,stringBuffer);
                        #endif

      if (error!=0)
        return (error);

      continue;
    }

    if (strcmp(buffer,"exit")==0)
    {
      GetToken(buffer);

      value = atoi(buffer);
      if ((value!=0) && (value!=1))
      {
        PrintErrorMessage('E',"InterpretString","syntax error (exit)");
        return(8510);                           /* syntax error */
      }

      dontexit = !value;

      continue;
    }

    if (strcmp(buffer,"mute")==0)
    {
      if ((error=EvaluateExpression(&result))!=DONE)
        return(error);

      if (result.ro.type!=NUMBERID)
        return (1234);                          /* syntax error: number for mutelevel expected */

      mutelevel = (INT) result.ro.value;

      SetStringValue(":oldmute",GetMuteLevel());
      SetMuteLevel(mutelevel);

      continue;
    }

    if ((strcmp(buffer,"execute")==0)||(strcmp(buffer,"ex")==0))
    {
      if ((error=GetFullPathName(buffer))!=DONE)
        return(error);

      filePtr = FOpenScript(buffer,"r");
      if (filePtr==NULL)
      {
        strcpy(filename,buffer);
        strcat(filename,".scr");
        filePtr = FOpenScript(filename,"r");
      }

      if (filePtr==NULL)
      {
        char buf[128];
        sprintf(buf,"could not open file %s",filename);
        PrintErrorMessage('E',"InterpretString",buf);
        return(8515);                           /* could not open file */
      }
      if ((error=setvbuf(filePtr,fileBuffer,_IOFBF,FILEBUFSIZE))!=0)
      {
        PrintErrorMessage('E',"InterpretString","could not allocate file buffer");
        return(8516);                           /* could not allocate file buffer */
      }

      rewind(filePtr);

      executePos0=executePos;
      while ((error=getc(filePtr))!=EOF)
      {
        if (executePos==cmdintbufsize-1)
        {
          fclose(filePtr);
          executePos=executePos0;                                       /* reset execute position */
          PrintErrorMessage('E',"InterpretString","buffer overflow");
          return(8517);                                 /* buffer overflow */
        }
        else if (error=='\r')
          /* convert PC-newlines into Mac ones */
          executeBuffer[executePos++]=(char) '\n';
        else
          executeBuffer[executePos++]=(char) error;
      }
      executeBuffer[executePos++]=(char) 0;

      fclose(filePtr);

      /* set :scriptname var */
      SetStringVar(":scriptname",buffer);

      /* execute buffer */
      error=InterpretCommand(&(executeBuffer[executePos0]));

      executePos=executePos0;                           /* reset execute position */
      if (error!=DONE)
        return(error);

      continue;
    }

    cmdEnd = cmdPtr;
    c=SkipBlanks();
    if (c!='=')
    {
      /* handle ug command */

      /* reset cmdPtr for correct expansion */
      cmdPtr=cmdStart;
      k=0;

      commandItem = GetCommand(buffer);

      /* search case INsensitive and permit abbreviations */
      if (commandItem==NULL)
        if ((commandItem=SearchUgCmd(buffer))!=NULL)
        {
          /* write correct commandname to cmdBuffer */
          strcpy(cmdBuffer,commandItem->v.name);
          k = strlen(cmdBuffer);

          /* set cmdPtr behind abbreviated cmd */
          cmdPtr = cmdEnd;
        }

      if (commandItem==NULL)
      {
        char buf[64];
        sprintf(buf,"command %s not found",buffer);
        PrintErrorMessage('E',"InterpretString",buf);
        return (8510);
      }

      /* expand arguments containig @ */
      do
      {
        if (strchr(terminators,(int) (c=*cmdPtr))!=NULL)
          c=(char) 0;

        if (c=='#')
        {
          /* comment: skip end of line */
          do
            cmdPtr++;
          while (((c=*cmdPtr)!='\n') && (c!=(char) 0));
        }

        if (c!=(char) 0)
          cmdPtr++;

        if (c=='@')
        {
          if ((error=GetToken(buffer))!=DONE)
          {
            PrintErrorMessage('E',"InterpretString","syntax error");
            return(error);                                      /* syntax error */
          }

          if ((sptr=GetStringVar(buffer))==NULL)
          {
            char buf[64];
            sprintf(buf,"variable %s not found",buffer);
            PrintErrorMessage('E',"InterpretString",buf);
            return(8522);                                       /* variable not found */
          }
          while ((c1=*(sptr++))!=(char) 0)
          {
            if (k==cmdintbufsize)
            {
              PrintErrorMessage('E',"InterpretString","command too long");
              return(8509);                                             /* command too long */
            }
            else
              cmdBuffer[k++]=c1;
          }
        }
        else
        {
          if (k==cmdintbufsize)
          {
            PrintErrorMessage('E',"InterpretString","command too long");
            return(8509);                                       /* command too long */
          }
          else
            cmdBuffer[k++]=c;
        }
      }
      while (c!=(char) 0);

      if (mutelevel>=0)
      {
        if (executePos!=0)
        {
          UserWrite(": ");
          UserWrite(cmdBuffer);
          UserWrite("\n");
        }
        else
        {
          WriteLogFile(": ");
          WriteLogFile(cmdBuffer);
          WriteLogFile("\n");
        }
      }
                        #ifdef ModelP
      error = ParExecCommand(cmdBuffer);
      /* if ((gresult==CMD_EXIT)||(gresult==CMD_DISASTROUS)) return(DONE); */
                        #else
      error = ExecCommand(cmdBuffer);
                        #endif
      if (dontexit && (error!=QUITCODE) && (error!=FATAL))
      {
        SetStringValue(":cmdstatus",error);
        error = DONE;
      }
      else if (error!=DONE)
      {
        if (GetMuteLevel()<0)
          SetMuteLevel(mutelevel=0);
        return(error);
      }
      continue;
    }
    else
    {
      /* There is only the possibility left that the statement
         is an assignment */
      cmdPtr++;
      if ((error=EvaluateExpression(&result))!=DONE)
        return(error);

      switch (result.ro.type)
      {
      case NUMBERID :
        sprintf(valueStr,"%-.14g",(double)result.ro.value);
                                        #ifdef ModelP
        PRINTDEBUG(ui,2,("%2d: set NUMBERID=%s\n",me,valueStr))
        sprintf(execCmdBuffer,"set %s %s",buffer,valueStr);
        error = ParExecCommand(execCmdBuffer);
                                        #else
        error = SetStringVar(buffer,valueStr);
                                        #endif
        break;

      case ALPHAID :
                                        #ifdef ModelP
        PRINTDEBUG(ui,2,("%2d: set ALPHAID=%s\n",me,result.so.sptr))
        sprintf(execCmdBuffer,"set %s %s",buffer,result.so.sptr);
        error = ParExecCommand(execCmdBuffer);
                                        #else
        error = SetStringVar(buffer,result.so.sptr);
                                        #endif
        break;

      case LSTRINGID :
                                        #ifdef ModelP
        sprintf(execCmdBuffer,"set %s ",buffer);
        strncat(execCmdBuffer,result.lo.sptr,result.lo.length);
        PRINTDEBUG(ui,2,("%2d: set LSTRINGID=%s\n",me,result.lo.sptr))
        error = ParExecCommand(execCmdBuffer);
                                        #else
        error = SetnStringVar(buffer,result.lo.sptr,result.lo.length);
                                        #endif
      }

      if (error)
        return (error);

      continue;
    }
  }
  while (TRUE);

  return(DONE);

}       /* InterpretString */



/****************************************************************************/
/** \brief Generate version string and print it

   This function generates the version string and prints it.

 */
/****************************************************************************/

static void PrintVersionString (void)
{
  char ver[128];
  int i,j,k;


  strcpy(ver,VERSION_STRING);
  for (i=0; i<100; i++)
  {
    if (ver[i] == '\0') break;
    if (ver[i] == '$') break;
  }
  k = 0;
  for (j=i+6; j<100; j++)
  {
    if (ver[j] == '$')
      k = 1;
    else
      ver[j-6-k] = ver[j];
    if (ver[j] == '\0') break;
  }

  /* print version */
  UserWrite(ver);
}


#ifdef ModelP

static INT InteractiveCommand(INT argc, char **argv)
{
  interactiveFlag=1;
  return OKCODE;
}

static INT NonInteractiveCommand(INT argc, char **argv)
{
  interactiveFlag=0;
  return OKCODE;
}

/****************************************************************************/
/** \brief Get next command to execute from master

   This function receives the next command which should be
   executed from master processor, unitl QUITCODE is returned.

 */
/****************************************************************************/

void NS_DIM_PREFIX ParCommandLoop (char *inpLine)
{
  INT error;
  char dummy[256];

  while (GetDoneFlag() == FALSE)
  {
    if (interactiveFlag) UserIn(dummy);
    error=ParExecCommand(inpLine);
    if (error==QUITCODE) SetDoneFlag();
  }

  return;
}
#endif

/****************************************************************************/
/** \brief Read commands from user interface

   \param argc - argument counter
   \param argv - argument vector

   This function reads commands from user interface and
   executes them until  QUITCODE is returned.

 */
/****************************************************************************/

void NS_DIM_PREFIX CommandLoop (int argc, char **argv)
{
  INT error;
  int i,kerr;
  char c,errLine[256],spcLine[256],buffer[256];
  char *inpLine;
  const char *strStart;
  int batch = FALSE;

  /* reset doneFlag */
  ResetDoneFlag();

  /* alloc input line buffer */
  if ((inpLine=(char *)malloc(cmdintbufsize))==NULL)
  {
    PrintErrorMessage('F',"CommandLoop()","could not allocate inpLine buffer");
    return;
  }
  inpLine[0] = (char) 0;

#ifdef ModelP
  if (me==master)
  {
    /* FOR MASTER PROCESSOR */
#endif

  for (i=1; i<argc; i++)
    if (argv[i][0]!='-')
      batch = TRUE;

  PrintVersionString();

  /* if (argc==-1): second start of CommandLoop */
  if (argc != -1)
  {
    /* execute init script */
    if (GetDefaultValue(DEFAULTSFILENAME,"initscript",buffer)==0)
    {
      char *cmds = strchr(buffer,'"');
      if (cmds!=NULL)
      {
        char *p = strchr(buffer+1,'"');
        if (p!=NULL) *p = '\0';
        strcpy(inpLine,cmds+1);
      }
      else
      {
        strcpy(inpLine,"execute ");
        strcat(inpLine,buffer);
      }
      error = InterpretCommand(inpLine);
      if (error==QUITCODE)
        SetDoneFlag();
    }
  }


  if (!batch)
  {
    while (GetDoneFlag() == FALSE)
    {
      if (UseWithPerl)
        WriteString("EOO\n");
      else
        WriteString(PROMPT);
#ifdef ModelP
      InterpretCommand("interactive");
#endif
      if (UserIn(inpLine)!=0)
      {
        PrintErrorMessage('E',"CommandLoop","process event error");
        continue;
      }
      if (GetDoneFlag() == TRUE) break;
#ifdef ModelP
      InterpretCommand("noninteractive");
#endif
      if ((error=InterpretCommand(inpLine))!=DONE)
      {
        if (error==QUITCODE)
          SetDoneFlag();
        else
        {
          UserWrite("Error position: ");

          /* search for beginning of line */
          strStart=cmdPtr;
          kerr=0;
          while (strStart>cmdStart)
          {
            if (*(--strStart)=='\n')
            {
              strStart++;
              break;
            }
            kerr++;
          }

          if (kerr<254)
          {
            int k=0;
            while (k<254)
            {
              c=*(strStart++);
              if ((c==(char) 0)||(c=='\n'))
                break;
              if (k<kerr)
              {
                if (c!='\t')
                  spcLine[k]=' ';
                else
                  spcLine[k]='\t';
              }
              errLine[k++]=c;
            }
            errLine[k++]='\n';
            errLine[k]=(char) 0;
            UserWrite(errLine);

            spcLine[kerr]=(char)0;
            UserWrite(spcLine);
            UserWrite("                ^\n");
          }
        }
      }
    }
  }
  else
  {
    i = 1;             /* first argument */
    while (i<argc && GetDoneFlag()==FALSE)
    {
      /* execute batch file */
      if (argv[i][0]!='-')
      {
        sprintf(inpLine,"execute %s\n",argv[i]);

        InterpretCommand(inpLine);                 /* execute command line argument */
        if (i + 1 < argc)
          if (strcmp(argv[i+1],"-noquit") == 0) {
            CommandLoop(-1,NULL);
            /* free command line buffer */
            free(inpLine);
            return;
          }
        InterpretCommand("quit\n");                /* end program */
        i++;
        continue;
      }
      /* set command from command line */
      if ((argv[i][0]=='-')&&(argv[i][1]=='S'))
      {
        if (i+1<argc)
        {
          sprintf(inpLine,"%s\n",(argv[i+1]));
          InterpretCommand(inpLine);
          i++;
        }
        else
          UserWrite("Error in command line option -S\n");
        i++;
        continue;
      }
      /* logon command from command line */
      if ((argv[i][0]=='-')&&(argv[i][1]=='L'))
      {
        if (i+1<argc)
        {
          sprintf(inpLine,"logon %s\n",(argv[i+1]));
          InterpretCommand(inpLine);
          i++;
        }
        else
          UserWrite("Error in command line option -L\n");
        i++;
        continue;
      }
      i++;
    }
  }

#ifdef ModelP
}
else
{
  /* FOR PROCESSORS WITH ME!=MASTER */

  ParCommandLoop(inpLine);

}
#endif

  /* free input line buffer */
  free(inpLine);

  /* reset doneFlag */
  ResetDoneFlag();

  /* call ExitUg() at the end of CommandLoop in order to avoid that
     the application programmer will forget to call it at the end of
     the application. */
#ifdef __GUI__
  {
    int gui_on = GUI_ON;
#ifdef ModelP
    Broadcast(&gui_on,sizeof(int));
#endif
    if (!gui_on)
#endif
  ExitUg();
#ifdef __GUI__
}
#endif
}


/****************************************************************************/
/*
   SetDoneFlag/ResetDoneFlag -

   SYNOPSIS:
   void SetDoneFlag (void); void ResetDoneFlag (void);

   PARAMETERS:
   \param none

   DESCRIPTION:
   Sets done flag to TRUE/FALSE

   RETURN VALUE:
   none
 */
/****************************************************************************/

void NS_DIM_PREFIX SetDoneFlag (void)
{
  doneFlag = TRUE;
}

void NS_DIM_PREFIX ResetDoneFlag (void)
{
  doneFlag = FALSE;
}

int NS_DIM_PREFIX GetDoneFlag (void)
{
  return(doneFlag);
}


/****************************************************************************/
/** \brief Initializes interpreter

   Allocates `execute` and `program` buffer, reads in `scriptpaths` from the defaults file,
   sets back the flag 'dontexit'.

   SEE ALSO:
   ReadSearchingPaths

   RETURN VALUE:
   INT
   .n Error code, if not zero.
   D*/
/****************************************************************************/


INT NS_DIM_PREFIX InitCommandInterpreter (INT argc, char **argv)
{
  char size[256];
  int i,ival;

  /* get cmdintbufsize from defaults file */
  if (GetDefaultValue(DEFAULTSFILENAME,"cmdintbufsize",size)==0)
  {
    sscanf(size," %d ",&ival);
    cmdintbufsize = ival;
  }

  /* alloc command buffer */
  if ((cmdBuffer=(char *)malloc(cmdintbufsize))==NULL)
  {
    PrintErrorMessage('F',"InitCommandInterpreter","could not allocate cmdBuffer buffer");
    return(__LINE__);
  }
  cmdBuffer[0] = (char) 0;

  /* alloc execute buffer */
  if ((executeBuffer=(char *)malloc(cmdintbufsize))==NULL)
  {
    PrintErrorMessage('F',"InitCommandInterpreter","could not allocate executeBuffer buffer");
    return(__LINE__);
  }
  executeBuffer[0] = (char) 0;

        #ifdef ModelP
  /* alloc buffer to broadcast command sequences */
  if ((execCmdBuffer=(char *)malloc(cmdintbufsize))==NULL)
  {
    PrintErrorMessage('F',"InitCommandInterpreter","could not allocate execCmdBuffer buffer");
    return(__LINE__);
  }
  execCmdBuffer[0] = (char) 0;
        #endif

  /* alloc program buffer */
  if ((programBuffer=(char *)malloc(programbufsize))==NULL)
  {
    PrintErrorMessage('F',"InitCommandInterpreter","could not allocate program buffer");
    return(__LINE__);
  }
  programBuffer[0] = (char) 0;

  /* read scriptpaths from defaults file (iff) */
  scriptpaths_set=FALSE;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"scriptpaths")==0)
    scriptpaths_set=TRUE;

  dontexit=FALSE;

  SetStringValue(":oldmute",GetMuteLevel());

  /* use with perl? */
  UseWithPerl=0;
  for (i=0; i<argc; i++)
    if(strcmp(argv[i],"-perl")==0)
      UseWithPerl=1;

#ifdef ModelP
  /* special commands to mark begin/end of interactive mode */
  if (CreateCommand("interactive", InteractiveCommand)==NULL) return(__LINE__);
  if (CreateCommand("noninteractive", NonInteractiveCommand)==NULL) return(__LINE__);
#endif

  /* return to application */
  return(0);
}


INT NS_DIM_PREFIX ExitCommandInterpreter()
{
  free(cmdBuffer);
  free(executeBuffer);
#ifdef ModelP
  free(execCmdBuffer);
#endif
  free(programBuffer);

  return 0;
}
