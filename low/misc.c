// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      misc.c                                                        */
/*                                                                          */
/* Purpose:   miscellaneous routines                                        */
/*                                                                          */
/* Author:      Klaus Johannsen                                             */
/*              Institut fuer Computeranwendungen                           */
/*              Universitaet Stuttgart                                      */
/*              Pfaffenwaldring 27                                          */
/*              70569 Stuttgart                                             */
/*              internet: ug@ica3.uni-stuttgart.de		                    */
/*                                                                          */
/* History:   08.12.94 begin, ug3-version                                   */
/*                                                                          */
/* Revision:  07.09.95                                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*              system include files                                        */
/*              application include files                                   */
/*                                                                          */
/****************************************************************************/

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#ifdef __NECSX4__
#include <sys/types.h>
#include <sys/syssx.h>
#endif

#include "compiler.h"
#include "general.h"
#include "misc.h"
#include "heaps.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*          in the corresponding include file!)                             */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

int UG_math_error = 0; /* This will be non zero after a math error occured  */
#ifndef ModelP
int me = 0;                     /* to have in the serial case this variable as a dummy */
int master = 0;         /* to have in the serial case this variable as a dummy */
int procs = 1;          /* to have in the serial case this variable as a dummy */
int _proclist_ = -1; /* to have in the serial case this variable as a dummy */
int _partition_ = 0; /* to have in the serial case this variable as a dummy */
#endif

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* forward declarations of macros                                           */
/*                                                                          */
/****************************************************************************/

#define MIN_DETERMINANT                 1e-8

/* for ExpandCShellVars */
#define CSHELL_VAR_BEGIN                "$("
#define CSHELL_VAR_BEGIN_LEN    2               /* strlen(CSHELL_VAR_BEGIN)			*/
#define CSHELL_VAR_END                  ")"
#define CSHELL_VAR_END_LEN              1               /* strlen(CSHELL_VAR_END)			*/

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****                                                                    ****/
/****        general routines                                            ****/
/****                                                                    ****/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   INT_2_bitpattern	- transform an INT into a bitpattern string

   SYNOPSIS:
   static void INT_2_bitpattern (INT n, char *text)

   PARAMETERS:
   .  n - integer to convert
   .  text - string of size >= 33 for conversion

   DESCRIPTION:
   This function transforms an INT into a bitpattern string consisting of 0s
   and 1s only.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void INT_2_bitpattern (INT n, char text[33])
{
  INT i;

  memset(text,'0',32*sizeof(char));

  for (i=0; i<32; i++)
    if ((1<<i)&n)
      text[31-i] = '1';
  text[32] = '\0';

  return;
}

/****************************************************************************/
/*																			*/
/* Function:  CenterInPattern												*/
/*																			*/
/* Purpose:   fill str up to PatLen with pattern and center text in it.		*/
/*			  Terminate the str with end									*/
/*																			*/
/* Input:	  s.a															*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   CenterInPattern	- compose a headline of chars with string centered

   SYNOPSIS:
   INT CenterInPattern (char *str, INT PatLen, const char *text, char p, const char *end)

   PARAMETERS:
   .  str - string to print to
   .  PatLen - width of headline
   .  text - text to center in headline
   .  p - char replicated for headline
   .  end - trailing string (optional)

   DESCRIPTION:
   This function transforms an INT into a bitpattern string consisting of 0s
   and 1s only.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

INT CenterInPattern (char *str, INT PatLen, const char *text, char p, const char *end)
{
  INT i,TextBegin,TextEnd,TextLen;

  TextLen   = strlen(text);
  TextBegin = (PatLen-TextLen)/2;
  TextEnd   = TextBegin+TextLen;

  if (TextLen>PatLen)
    return (CenterInPattern(str,PatLen," text too long ",p,end));

  for (i=0; i<TextBegin-1; i++)
    str[i] = p;
  str[i] = ' ';
  for (i=TextBegin; i<TextEnd; i++)
    str[i] = *(text++);
  str[i++] = ' ';
  for (; i<PatLen; i++)
    str[i] = p;
  str[PatLen] = '\0';
  if (end!=NULL)
    strcat(str,end);

  return (0);
}

/****************************************************************************/
/*D
   expandfmt - Expand (make explicit) charset-ranges in scanf

   SYNOPSIS:
   char *expandfmt (const char *fmt)

   PARAMETERS:
   .  fmt - pointer to char (const)

   DESCRIPTION:
   This function expands (make explicit) charset-ranges in scanf.
   For example '%<number>[...a-d...]' --> '%<number>[...abcd...]'.

   RETURN VALUE:
   char*
   .n        new pointer to char
   D*/
/****************************************************************************/

/* install a user math error handler */
int matherr(struct exception *x)
{
  /* your math error handling here */
  UG_math_error = 1;       /* lets be more specific later ... */
  return(0);   /* proceed with standard messages */
}
#ifdef __HP__
/* HP-UX 10.20 expects the name with an "_" as prefix for the math error handler */
int _matherr(struct exception *x)
{
  /* your math error handling here */
  UG_math_error = 1;       /* lets be more specific later ... */
  return(0);   /* proceed with standard messages */
}
#endif

#define FMTBUFFSIZE         1031

static char newfmt[FMTBUFFSIZE];

char *expandfmt (const char *fmt)
{
  const char *pos;
  char *newpos;
  char leftchar,rightchar;
  int newlen;

  /* calculate min size of newfmt */
  newlen = strlen(fmt);
  assert (newlen<FMTBUFFSIZE-1);

  pos    = fmt;
  newpos = newfmt;

  /* scan fmt for %<number>[ */
  while (*pos!='\0')
  {
    /* copy til '%' */
    while (*pos!='%' && *pos!='\0')
      *(newpos++) = *(pos++);

    if (*pos=='\0')
      break;

    *(newpos++) = *(pos++);

    /* copy til !isdigit */
    while (isdigit(*pos) && *pos!='\0')
      *(newpos++) = *(pos++);

    if (*pos=='\0')
      break;

    if (*pos!='[')
      continue;

    *(newpos++) = *(pos++);

    /* ']' following '[' is included in the charset */
    if ((*pos)==']')
      *(newpos++) = *(pos++);

    /* '^]' following '[' is included in the charset */
    else if ((*pos)=='^' && (*(pos+1))==']')
    {
      *(newpos++) = *(pos++);
      *(newpos++) = *(pos++);
    }

    while (*pos!=']' && *pos!='\0')
    {
      /* now we are inside the charset '[...]': */

      /* treat character ranges indicated by '-' */
      while (*pos!='-' && *pos!=']' && *pos!='\0')
        *(newpos++) = *(pos++);

      if (*pos=='\0')
        break;

      if ((*pos)==']')
        continue;

      /* gotya: is left char < right char? */
      leftchar  = *(pos-1);
      rightchar = *(pos+1);

      if (leftchar=='[' || rightchar==']')
      {
        *(newpos++) = *(pos++);
        continue;
      }

      if (leftchar>=rightchar)
      {
        *(newpos++) = *(pos++);
        continue;
      }

      if (leftchar+1==rightchar)
      {
        /* for example '...b-c...' */
        pos++;
        continue;
      }

      /* calc new size and expand range */
      newlen += rightchar-leftchar-2;
      assert (newlen<FMTBUFFSIZE-1);

      leftchar++;
      pos++;

      while (leftchar<rightchar)
      {
        if (leftchar=='^' || leftchar==']')
        {
          leftchar++;
          continue;
        }
        *(newpos++) = leftchar++;
      }
    }
  }

  *newpos = '\0';

  return (newfmt);
}

char *ExpandCShellVars (char *string)
{
  if (strstr(string,CSHELL_VAR_BEGIN)!=NULL)
  {
    /* shell var reference contained: copy string and expand vars */
    char *copy = StrDup(string);
    char *p0  = copy;                           /* current pos  */
    char *p1;                                           /* end of token */

    string[0] = '\0';

    while ((p1 = strstr(p0,CSHELL_VAR_BEGIN))!=NULL)
    {
      char *var;

      /* enclose current pos (p0) to CSHELL_VAR_BEGIN between p0,p1 */
      *p1 = '\0';

      /* copy verbatim */
      strcat(string,p0);

      /* advance p0 to begin of shell var name */
      p0 = p1+CSHELL_VAR_BEGIN_LEN;

      /* enclose current pos (p0) to CSHELL_VAR_END between p0,p1 */
      p1 = strstr(p0,CSHELL_VAR_END);
      if (p1==NULL)
      {
        free(copy);
        return NULL;
      }
      *p1 = '\0';

      /* copy shell variable iff */
      var = getenv(p0);
      if (var==NULL)
      {
        free(copy);
        return NULL;
      }
      strcat(string,var);

      /* advance p0 */
      p0 = p1+CSHELL_VAR_END_LEN;
    }
    /* copy remainder */
    strcat(string,p0);
    free(copy);
  }
  return string;
}

/****************************************************************************/
/*D
   StrTok - Copy a token out of a string (ANSI-C)

   SYNOPSIS:
   char *StrTok (char *s, const char *ct)

   PARAMETERS:
   .  s -  pointer to char
   .  ct - pointer to char (const)

   DESCRIPTION:
   This function copies a token out of a string.

   See also ANSI-C for description of this function.

   RETURN VALUE:
   char
   .n    pointer to char, modified string
   D*/
/****************************************************************************/

char *StrTok (char *s, const char *ct)
{
  static char *e;
  char *b;
  INT i, flag;

  if (s != NULL)
    b = s-1;
  else
    b = e+1;

  flag = 0;
  while (flag == 0)
  {
    b += 1;
    for (i=0; i<strlen(ct); i++)
      if ((*(ct+i)) == (*b))
        flag = 1;
  }
  e = b;
  flag = 0;
  while (flag == 0)
  {
    e += 1;
    for (i=0; i<strlen(ct); i++)
      if ((*(ct+i)) == (*e))
        flag = 1;
  }
  *e = '\0';

  return (b);
}

/****************************************************************************/
/*D
   strntok - Split a string into tokens each of maximal length 'n+1'

   SYNOPSIS:
   const char *strntok (const char *str, const char *sep, int n, char *token)

   PARAMETERS:
   .  str -   pointer to char (const)
   .  sep -   pointer to char (const)
   .  n -     integer, number of chars in token
   .  token - pointer to char

   DESCRIPTION:
   This function splits a string into tokens each of maximal length 'n+1'.
   A pointer to the next char following the token (its a sep) is returned.
   NB: possibly check whether the returned char is a sep.
   If not: the token was to long and only the first n chars where copied!

   RETURN VALUE:
   char
   .n     pointer to token
   .n     NULL if token larger than n.
   D*/
/****************************************************************************/

const char *strntok (const char *str, const char *sep, int n, char *token)
{
  int i;

  /* scan while current char is a seperator */
  while ((*str!='\0') && (strchr(sep,*str)!=NULL)) str++;

  /* copy into token */
  for (i=0; i<n; i++,str++)
    if ((*str!='\0') && (strchr(sep,*str)==NULL))
      token[i] = *str;
    else
      break;

  if (strchr(sep,*str)==NULL)
    return (NULL);                    /* ERROR: token too long! */

  /* 0-terminate string */
  token[i] = '\0';

  return (str);
}

/****************************************************************************/
/*D
   StrDup - duplicate string to memory allocated with malloc

   SYNOPSIS:
   char *StrDup (const char *s)

   PARAMETERS:
   .  s - string to duplicate

   DESCRIPTION:
   This function duplicates a string to memory allocated with malloc.

   RETURN VALUE:
   char*
   .n        pointer to new string
   D*/
/****************************************************************************/

char *StrDup (const char *s)
{
  char *p;

  p = (char*) malloc(strlen(s)+1);
  if (p!=NULL)
    strcpy(p,s);

  return(p);
}

/****************************************************************************/

static void Copy (char *a, const char *b, INT size)
{
  INT i;

  for (i=0; i<size; i++)
    a[i] = b[i];
}

/****************************************************************************/
/*D
   QSort - Sorting routine (standard function)

   SYNOPSIS:
   void QSort (void *base,INT n,INT size,int (*cmp)(const void *,const void *))

   PARAMETERS:
   .  base - pointer to void, field to be sorted
   .  n -    integer, length of string
   .  size - integer, number of characters to be sorted
   .  cmp -  pointer to function with two arguments

   DESCRIPTION:
   This function sorts the values returned by a function given as argument.

   See also standard description of 'QSort'.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/
void QSort (void *base, INT n, INT size, int (*cmp)(const void *, const void *))
{
  INT i, j;
  char *Base, v[4], t[4];
  INT cmp1, cmp2, compare,flag;

  if (n<2) return;
  flag=0;
  Base = (char*)base;
  Copy(v,Base+(n-1)*size,size);
  for (i=-1, j=n-1;;)
  {
    while (++i<n-1)
      if ((cmp1=cmp( (void *)v, (void *)(Base+i*size))) <= 0)
        break;
    while ( --j>0)
      if ((cmp2=cmp( (void *)v, (void *)(Base+j*size))) >= 0)
        break;
    if (i>=j)
      break;
    compare  = (cmp1<0);
    compare |= ((cmp2>0)<<1);
    switch (compare)
    {
    case (0) :
      QSort ( (void *)(Base+i*size), j-i+1, size, cmp);
      i--;
      while (++i<n-1)
        if (cmp( (void *)v, (void *)(Base+i*size)) < 0)
          break;
      flag=1;
      break;
    case (1) :
      Copy(t,Base+i*size,size);
      Copy(Base+i*size,Base+j*size,size);
      Copy(Base+j*size,t,size);
      i--;
      break;
    case (2) :
      Copy(t,Base+i*size,size);
      Copy(Base+i*size,Base+j*size,size);
      Copy(Base+j*size,t,size);
      j++;
      break;
    case (3) :
      Copy(t,Base+i*size,size);
      Copy(Base+i*size,Base+j*size,size);
      Copy(Base+j*size,t,size);
      break;
    }
    if (flag) break;
  }
  Copy(t,Base+i*size,size);
  Copy(Base+i*size,v,size);
  Copy(Base+(n-1)*size,t,size);
  QSort ( base, i, size, cmp);
  QSort ( (void *)(Base+(i+1)*size), n-i-1, size, cmp);
}

/****************************************************************************/
/*D
   SelectionSort - Sorting routine (standard)

   SYNOPSIS:
   void SelectionSort (void *base,INT n,INT size,int (*cmp)(const void *,const void *))

   PARAMETERS:
   .  base - pointer to void, field to be sorted
   .  n -    integer, length of string
   .  size - integer, number of characters to be sorted
   .  cmp -  pointer to function with two arguments

   DESCRIPTION:
   This function sorts the arguments of the function.

   See also standard description of 'SelectionSort'.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void SelectionSort (void *base, INT n, INT size, int (*cmp)(const void *, const void *))
{
  INT i,j,k1,k2,s;
  char Smallest[4];
  char *Base;

  if (n<2) return;
  Base = (char*)base;

  for (i=0; i<n; i++)
  {
    Copy(Smallest,Base+i*size,size);
    k1=i;
    for (s=0; s<n-i; s++)
    {
      k2=k1;
      for (j=i; j<n; j++)
      {
        if (j==k1) continue;
        if (cmp( (void *)Smallest, (void *)(Base+j*size)) > 0)
        {
          Copy(Smallest,Base+j*size,size);
          k1=j;
        }
      }
      if (k1==k2) break;
    }
    Copy(Smallest,Base+i*size,size);
    Copy(Base+i*size,Base+k1*size,size);
    Copy(Base+k1*size,Smallest,size);
  }
}


/****************************************************************************/
/*D
    ReadMemSizeFromString - Convert a (memory)size specification from String to MEM (long int)

   SYNOPSIS:
   INT ReadMemSizeFromString (const char *s, MEM *mem_size );

   PARAMETERS:
   .  s - input string
   .  mem_size - the specified mem size in byte

   DESCRIPTION:
   This function converts a (memory)size specification from String to type MEM (an integer type).
   The size specification contains an integer number followed by an optional unit specifier:
      G for gigabyte
      M for megabyte
      K for kilobyte
   (also the lower case chars are recognized).

   EXAMPLE:
      "10M" is converted to 10485760 (10 mega byte).

   RETURN VALUE:
   INT: 0 ok
        1 integer could not be read
        2 invalid unit specifier

   SEE ALSO:
   MEM, WriteMemSizeToString
   D*/
/****************************************************************************/

INT ReadMemSizeFromString (const char *s, MEM *mem_size )
{
  float mem;

  if (sscanf( s, "%e",&mem)!=1)
    return(1);

  switch( s[strlen(s)-1] )
  {
  case 'k' : case 'K' :             /* check for [kK]ilobyte-notation */
    *mem_size = floor(mem * KBYTE);
    return(0);
  case 'm' : case 'M' :             /* check for [mM]egabyte-notation */
    *mem_size = floor(mem * MBYTE);
    return(0);
  case 'g' : case 'G' :             /* check for [gG]igabyte-notation */
    *mem_size = floor(mem * GBYTE);
    return(0);
  case '0' : case '1' : case '2' : case '3' : case '4' : case '5' : case '6' : case '7' : case '8' : case '9' :     /* no recognized mem unit character recognized */
    *mem_size = floor(mem * 1);
    return(0);
  default :              /* unknown mem unit character */
    return(2);
  }
}

/****************************************************************************/
/*D
    WriteMemSizeToString - Convert a (memory)size MEM to string

   SYNOPSIS:
   INT WriteMemSizeToString (MEM mem_size, char *s)

   PARAMETERS:
   .  s - input string
   .  mem_size - the specified mem size in byte

   DESCRIPTION:
   This function writes a MEM size in MBytes to string in a format that is recognized by
   WriteMemSizeToString.

   RETURN VALUE:
   INT: 0 ok

   SEE ALSO:
   MEM, ReadMemSizeFromString
   D*/
/****************************************************************************/

INT WriteMemSizeToString (MEM mem_size, char *s)
{
  float mem = mem_size;

  /* write mem size in units of MByte */
  sprintf(s, "%g M",mem/MBYTE);
  return 0;
}


/* may be later on ugshell better ... */
#define UserWriteF printf

INT MemoryParameters (void)
{
  /* integer data types */
  char charType;
  short shortType;
  int intType;
  long longType;
  unsigned char ucharType;
  unsigned short ushortType;
  unsigned int uintType;
  unsigned long ulongType;

  /* floating point data types */
  float floatType;
  double doubleType;

  /* pointer data types */
  void *                          ptrType;

  int byteorder = 1;
  char hline[40]="---------------------------";

  UserWriteF("MEMORY specific parameters\n");

  UserWriteF("\n  Data type sizes are:\n");

  /* sizes of integer types */
  UserWriteF("    for integer types\n");
  UserWriteF("        type        | %.5s | %.5s | %.5s | %.5s\n","char","short","int","long");
  UserWriteF("    size   (signed) |   %2d |    %2d |  %2d |   %2d\n",
             sizeof(charType),sizeof(shortType),sizeof(intType),sizeof(longType));
  UserWriteF("    size (unsigned) |	%2d |    %2d |  %2d |   %2d\n",
             sizeof(ucharType),sizeof(ushortType),sizeof(uintType),sizeof(ulongType));

  /* sizes of floating point types */
  UserWriteF("    for floating point types\n");
  UserWriteF("        type | %6s | %6s\n","float","double");
  UserWriteF("        size |	   %2d |     %2d\n",
             sizeof(floatType),sizeof(doubleType));

  /* size of pointer types */
  UserWriteF("    for pointer types\n");
  UserWriteF("        ptr=%d\n",sizeof(ptrType));

  UserWriteF("\n  Ranges of data types are:\n");
  UserWriteF("    for integer types\n");
  UserWriteF("         type | %25s | %25s | %25s\n","<type>_MIN","<type>_MAX","U<type>_MAX");
  UserWriteF("         -------%.25s---%.25s---%.25s\n",hline,hline,hline);
  UserWriteF("         CHAR | %25d | %25d | %25d\n",CHAR_MIN,CHAR_MAX,UCHAR_MAX);
  UserWriteF("        SHORT | %25d | %25d | %25d\n",SHRT_MIN,SHRT_MAX,USHRT_MAX);
  UserWriteF("          INT | %25d | %25d | %25d\n",INT_MIN,INT_MAX,UINT_MAX);
  UserWriteF("         LONG | %25d | %25d | %25d\n",LONG_MIN,LONG_MAX,ULONG_MAX);

  UserWriteF("    for floating point types\n");
  UserWriteF("         type | %25s | %25s\n","<type>_MIN","<type>_MAX");
  UserWriteF("         -------%.25s---%.25s---%.25s\n",hline,hline,hline);
  UserWriteF("         FLT  | %.19E | %.19E\n",FLT_MIN,FLT_MAX);
  UserWriteF("         DBL  | %.19E| %.19E\n",DBL_MIN,DBL_MAX);

  UserWriteF("\n  Alignment and byteorder are:\n");
  UserWriteF("    alignment=%d byteorder=%s\n",ALIGNMENT,
             (((*((char *)&byteorder))==1) ? "BIGENDIAN" : "LITTLEENDIAN"));

  return(0);
}

#ifdef __NECSX4__
/* special high performance time system for NEC SX4 */
DOUBLE nec_clock( void )
{
  struct htms timebuf;
  DOUBLE dtime;

  if (syssx (HTIMES, (struct htms *)&timebuf) < 0)
    return -1.0;
  dtime = (timebuf.hutime / 1000000.0)+(timebuf.hstime / 1000000.0);
  return (dtime);
}
#endif
