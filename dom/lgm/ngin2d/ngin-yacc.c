// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*  A Bison parser, made from ngin.y
   by  GNU Bison version 1.25
 */

#define YYBISON 1  /* Identify Bison output.  */

#define ngparse ngparse
#define nglex nglex
#define ngerror ngerror
#define nglval nglval
#define ngchar ngchar
#define ngdebug ngdebug
#define ngnerrs ngnerrs
#define DOUBLE_VALUE    258
#define INT_VALUE       259
#define INODE   260
#define BNODE   261
#define LINE    262
#define ELEM    263
#define SIDE    264
#define TEND    265

#line 1 "ngin.y"

/****************************************************************************/
/*                                                                          */
/* File:      grid.y                                                        */
/*                                                                          */
/* Purpose:   parser for gridfiles                                          */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            Institut fuer Computeranwendungen                             */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: ug@ica3.uni-stuttgart.de                            */
/*                                                                          */
/* History:   19.2.98 begin,                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "ng.h"

#define alloca(p)               malloc(p)
#define SP_COPY(d,s)    {(d)->surf_id=(s)->surf_id; \
                         (d)->tri_id=(s)->tri_id; \
                         (d)->local[0]=(s)->local[0]; \
                         (d)->local[1]=(s)->local[1];}
#define LP_COPY(d,s)    {(d)->line_id=(s)->line_id; (d)->local=(s)->local;}

static LINE_POSITION LinePos;
static BND_NODE BndNode;
static INNER_NODE InnerNode;
static ELEM_SIDE ElemSide;
static NG_ELEMENT Elem;




#line 43 "ngin.y"
typedef union
{
  /* put RCS string here in order to get it into yacc-generated header file
     static char RCS_ID("$Header: /hosts/dom/cvs/df/gen/problems/dfcfg/dfcfg.y,v 0
     1998/02/20 16:58:46 birken Exp $",UG_RCS_STRING);
   */


  /* transfer lex->yacc */
  double dval;
  long ival;

  LINE_POSITION *lp;
  BND_NODE *bs;
  INNER_NODE *in;
  ELEM_SIDE *es;
  NG_ELEMENT *el;
} YYSTYPE;
#include <stdio.h>

#ifndef __cplusplus
#ifndef __STDC__
#define const
#endif
#endif



#define YYFINAL         50
#define YYFLAG          -32768
#define YYNTBASE        11

#define YYTRANSLATE(x) ((unsigned)(x) <= 265 ? ngtranslate[x] : 29)

static const char ngtranslate[] = {     0,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
                                        2,     2,     2,     2,     2,     1,     2,     3,     4,     5,
                                        6,     7,     8,     9,    10};

#if YYDEBUG != 0
static const short ngprhs[] = {     0,
                                    0,     4,     7,     9,    12,    13,    14,    20,    22,    25,
                                    29,    34,    40,    45,    51,    55,    57,    60,    65,    67,
                                    70,    71,    72,    80,    82,    85,    89,    91,    93};

static const short ngrhs[] = {    21,
                                  19,    12,     0,    21,    12,     0,    13,     0,    12,    13,
                                  0,     0,     0,     8,    14,    16,    15,    10,     0,    17,
                                  0,    17,    18,     0,    17,    18,    18,     0,    17,    18,
                                  18,    18,     0,    17,    18,    18,    18,    18,     0,    28,
                                  28,    28,    28,     0,    28,    28,    28,    28,    28,     0,
                                  9,    28,    28,     0,    20,     0,    19,    20,     0,     5,
                                  27,    27,    10,     0,    22,     0,    21,    22,     0,     0,
                                  0,     6,    27,    27,    23,    25,    24,    10,     0,    26,
                                  0,    25,    26,     0,     7,    28,    27,     0,     4,     0,
                                  3,     0,     4,     0};

#endif

#if YYDEBUG != 0
static const short ngrline[] = { 0,
                                 75,    77,    79,    81,    83,    87,    91,    92,    94,    95,
                                 96,    97,    99,   108,   118,   126,   128,   130,   138,   140,
                                 142,   149,   151,   152,   158,   164,   171,   173,   175};
#endif


#if YYDEBUG != 0 || defined (YYERROR_VERBOSE)

static const char * const ngtname[] = {   "$","error","$undefined.","DOUBLE_VALUE",
                                          "INT_VALUE","INODE","BNODE","LINE","ELEM","SIDE","TEND","Grid","ElemList","Elem",
                                          "@1","@2","ElemSpec","InnerElemSpec","ElemSide","InnerNodeList","InnerNode",
                                          "BndNodeList","BndNode","@3","@4","BndSpec","LinePosition","Coord","Id", NULL};
#endif

static const short ngr1[] = {     0,
                                  11,    11,    12,    12,    14,    15,    13,    16,    16,    16,
                                  16,    16,    17,    17,    18,    19,    19,    20,    21,    21,
                                  23,    24,    22,    25,    25,    26,    27,    27,    28};

static const short ngr2[] = {     0,
                                  3,     2,     1,     2,     0,     0,     5,     1,     2,     3,
                                  4,     5,     4,     5,     3,     1,     2,     4,     1,     2,
                                  0,     0,     7,     1,     2,     3,     1,     1,     1};

static const short ngdefact[] = {     0,
                                      0,     0,    19,    28,    27,     0,     0,     5,     2,     3,
                                      0,    16,    20,    21,     0,     0,     4,     1,    17,     0,
                                      0,    29,     6,     8,     0,     0,    22,    24,    18,     0,
                                      0,     9,     0,     0,     0,    25,     7,     0,    10,     0,
                                      26,    23,    15,    11,    13,    12,    14,     0,     0,     0};

static const short ngdefgoto[] = {    48,
                                      9,    10,    16,    30,    23,    24,    32,    11,    12,     2,
                                      3,    20,    35,    27,    28,     6,    25};

static const short ngpact[] = {     6,
                                    3,    13,-32768,-32768,-32768,     3,     3,-32768,     5,-32768,
                                    15,-32768,-32768,-32768,     3,    20,-32768,     5,-32768,    18,
                                    16,-32768,-32768,    21,    20,    20,    18,-32768,-32768,    17,
                                    20,    21,    20,     3,    19,-32768,-32768,    20,    21,    20,
                                    -32768,-32768,-32768,    21,    20,-32768,-32768,    31,    32,-32768};

static const short ngpgoto[] = {-32768,
                                22,    -4,-32768,-32768,-32768,-32768,   -28,-32768,    23,-32768,
                                33,-32768,-32768,-32768,     9,    -6,   -23};


#define YYLAST          36


static const short ngtable[] = {    14,
                                    15,    33,    34,    39,    17,     4,     5,    38,    21,    40,
                                    44,     1,     8,    17,    43,    46,    45,     7,     1,     7,
                                    8,    47,     8,    22,    26,    29,    37,    41,    42,    31,
                                    49,    50,    18,    19,    13,    36};

static const short ngcheck[] = {     6,
                                     7,    25,    26,    32,     9,     3,     4,    31,    15,    33,
                                     39,     6,     8,    18,    38,    44,    40,     5,     6,     5,
                                     8,    45,     8,     4,     7,    10,    10,    34,    10,     9,
                                     0,     0,    11,    11,     2,    27};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/share/bison.simple"

/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

#ifndef alloca
#ifdef __GNUC__
#define alloca __builtin_alloca
#else /* not GNU C.  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi)
#include <alloca.h>
#else /* not sparc */
#if defined (MSDOS) && !defined (__TURBOC__)
#include <malloc.h>
#else /* not MSDOS, or __TURBOC__ */
#if defined(_AIX)
#include <malloc.h>
 #pragma alloca
#else /* not MSDOS, __TURBOC__, or _AIX */
#ifdef __hpux
#ifdef __cplusplus
extern "C" {
  void *alloca (unsigned int);
};
#else /* not __cplusplus */
void *alloca ();
#endif /* not __cplusplus */
#endif /* __hpux */
#endif /* not _AIX */
#endif /* not MSDOS, or __TURBOC__ */
#endif /* not sparc.  */
#endif /* not GNU C.  */
#endif /* alloca not defined.  */

/* This is the parser code that is written into each bison parser
   when the %semantic_parser declaration is not specified in the grammar.
   It was written by Richard Stallman by simplifying the hairy parser
   used when %semantic_parser is specified.  */

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions, each action
   as one case of the switch.  */

#define ngerrok         (ngerrstatus = 0)
#define ngclearin       (ngchar = YYEMPTY)
#define YYEMPTY         -2
#define YYEOF           0
#define YYACCEPT        return (0)
#define YYABORT         return (1)
#define YYERROR         goto ngerrlab1
/* Like YYERROR except do call ngerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL          goto ngerrlab
#define YYRECOVERING()  (!!ngerrstatus)
#define YYBACKUP(token, value) \
  do                                                              \
    if (ngchar == YYEMPTY && nglen == 1)                          \
    { ngchar = (token), nglval = (value);                       \
      ngchar1 = YYTRANSLATE (ngchar);                           \
      YYPOPSTACK;                                               \
      goto ngbackup;                                            \
    }                                                           \
    else                                                          \
    { ngerror ("syntax error: cannot back up"); YYERROR; }      \
  while (0)

#define YYTERROR        1
#define YYERRCODE       256

#ifndef YYPURE
#define YYLEX           nglex()
#endif

#ifdef YYPURE
#ifdef YYLSP_NEEDED
#ifdef YYLEX_PARAM
#define YYLEX           nglex(&nglval, &nglloc, YYLEX_PARAM)
#else
#define YYLEX           nglex(&nglval, &nglloc)
#endif
#else /* not YYLSP_NEEDED */
#ifdef YYLEX_PARAM
#define YYLEX           nglex(&nglval, YYLEX_PARAM)
#else
#define YYLEX           nglex(&nglval)
#endif
#endif /* not YYLSP_NEEDED */
#endif

/* If nonreentrant, generate the variables here */

#ifndef YYPURE

int ngchar;                     /*  the lookahead symbol		*/
YYSTYPE nglval;                 /*  the semantic value of the		*/
                                /*  lookahead symbol			*/

#ifdef YYLSP_NEEDED
YYLTYPE nglloc;                 /*  location data for the lookahead	*/
                                /*  symbol				*/
#endif

int ngnerrs;                    /*  number of parse errors so far       */
#endif  /* not YYPURE */

#if YYDEBUG != 0
int ngdebug;                    /*  nonzero means print parse trace	*/
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif

/*  YYINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

#ifndef YYPARSE_RETURN_TYPE
#define YYPARSE_RETURN_TYPE int
#endif

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
YYPARSE_RETURN_TYPE ngparse (void);
#endif

#if __GNUC__ > 1                /* GNU C and GNU C++ define this.  */
#define __ng_memcpy(TO,FROM,COUNT)      __builtin_memcpy(TO,FROM,COUNT)
#else                           /* not GNU C or C++ */
#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__ng_memcpy (to, from, count)
char *to;
char *from;
int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__ng_memcpy (char *to, char *from, int count)
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif
#endif

#line 196 "/usr/share/bison.simple"

/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into ngparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
#ifdef __cplusplus
#define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else /* not __cplusplus */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
#endif /* not __cplusplus */
#else /* not YYPARSE_PARAM */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif /* not YYPARSE_PARAM */

YYPARSE_RETURN_TYPE
ngparse(YYPARSE_PARAM_ARG)
YYPARSE_PARAM_DECL
{
  register int ngstate;
  register int ngn;
  register short *ngssp;
  register YYSTYPE *ngvsp;
  int ngerrstatus;      /*  number of tokens to shift before error messages enabled */
  int ngchar1 = 0;              /*  lookahead token as an internal (translated) token number */

  short ngssa[YYINITDEPTH];     /*  the state stack			*/
  YYSTYPE ngvsa[YYINITDEPTH];   /*  the semantic value stack		*/

  short *ngss = ngssa;          /*  refer to the stacks thru separate pointers */
  YYSTYPE *ngvs = ngvsa;        /*  to allow ngoverflow to reallocate them elsewhere */

#ifdef YYLSP_NEEDED
  YYLTYPE nglsa[YYINITDEPTH];   /*  the location stack			*/
  YYLTYPE *ngls = nglsa;
  YYLTYPE *nglsp;

#define YYPOPSTACK   (ngvsp--, ngssp--, nglsp--)
#else
#define YYPOPSTACK   (ngvsp--, ngssp--)
#endif

  int ngstacksize = YYINITDEPTH;

#ifdef YYPURE
  int ngchar;
  YYSTYPE nglval;
  int ngnerrs;
#ifdef YYLSP_NEEDED
  YYLTYPE nglloc;
#endif
#endif

  YYSTYPE ngval;                /*  the variable used to return		*/
                                /*  semantic values from the action	*/
                                /*  routines				*/

  int nglen;

#if YYDEBUG != 0
  if (ngdebug)
    fprintf(stderr, "Starting parse\n");
#endif

  ngstate = 0;
  ngerrstatus = 0;
  ngnerrs = 0;
  ngchar = YYEMPTY;             /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  ngssp = ngss - 1;
  ngvsp = ngvs;
#ifdef YYLSP_NEEDED
  nglsp = ngls;
#endif

  /* Push a new state, which is found in  ngstate  .  */
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.  */
ngnewstate:

  *++ngssp = ngstate;

  if (ngssp >= ngss + ngstacksize - 1)
  {
    /* Give user a chance to reallocate the stack */
    /* Use copies of these so that the &'s don't force the real ones into memory. */
    YYSTYPE *ngvs1 = ngvs;
    short *ngss1 = ngss;
#ifdef YYLSP_NEEDED
    YYLTYPE *ngls1 = ngls;
#endif

    /* Get the current used size of the three stacks, in elements.  */
    int size = ngssp - ngss + 1;

#ifdef ngoverflow
    /* Each stack pointer address is followed by the size of
       the data in use in that stack, in bytes.  */
#ifdef YYLSP_NEEDED
    /* This used to be a conditional around just the two extra args,
       but that might be undefined if ngoverflow is a macro.  */
    ngoverflow("parser stack overflow",
               &ngss1, size * sizeof (*ngssp),
               &ngvs1, size * sizeof (*ngvsp),
               &ngls1, size * sizeof (*nglsp),
               &ngstacksize);
#else
    ngoverflow("parser stack overflow",
               &ngss1, size * sizeof (*ngssp),
               &ngvs1, size * sizeof (*ngvsp),
               &ngstacksize);
#endif

    ngss = ngss1; ngvs = ngvs1;
#ifdef YYLSP_NEEDED
    ngls = ngls1;
#endif
#else /* no ngoverflow */
      /* Extend the stack our own way.  */
    if (ngstacksize >= YYMAXDEPTH)
    {
      ngerror("parser stack overflow");
      return 2;
    }
    ngstacksize *= 2;
    if (ngstacksize > YYMAXDEPTH)
      ngstacksize = YYMAXDEPTH;
    ngss = (short *) alloca (ngstacksize * sizeof (*ngssp));
    __ng_memcpy ((char *)ngss, (char *)ngss1, size * sizeof (*ngssp));
    ngvs = (YYSTYPE *) alloca (ngstacksize * sizeof (*ngvsp));
    __ng_memcpy ((char *)ngvs, (char *)ngvs1, size * sizeof (*ngvsp));
#ifdef YYLSP_NEEDED
    ngls = (YYLTYPE *) alloca (ngstacksize * sizeof (*nglsp));
    __ng_memcpy ((char *)ngls, (char *)ngls1, size * sizeof (*nglsp));
#endif
#endif /* no ngoverflow */

    ngssp = ngss + size - 1;
    ngvsp = ngvs + size - 1;
#ifdef YYLSP_NEEDED
    nglsp = ngls + size - 1;
#endif

#if YYDEBUG != 0
    if (ngdebug)
      fprintf(stderr, "Stack size increased to %d\n", ngstacksize);
#endif

    if (ngssp >= ngss + ngstacksize - 1)
      YYABORT;
  }

#if YYDEBUG != 0
  if (ngdebug)
    fprintf(stderr, "Entering state %d\n", ngstate);
#endif

  goto ngbackup;
ngbackup:

  /* Do appropriate processing given the current state.  */
  /* Read a lookahead token if we need one and don't already have one.  */
  /* ngresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  ngn = ngpact[ngstate];
  if (ngn == YYFLAG)
    goto ngdefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* ngchar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (ngchar == YYEMPTY)
  {
#if YYDEBUG != 0
    if (ngdebug)
      fprintf(stderr, "Reading a token: ");
#endif
    ngchar = YYLEX;
  }

  /* Convert token to internal form (in ngchar1) for indexing tables with */

  if (ngchar <= 0)              /* This means end of input. */
  {
    ngchar1 = 0;
    ngchar = YYEOF;             /* Don't call YYLEX any more */

#if YYDEBUG != 0
    if (ngdebug)
      fprintf(stderr, "Now at end of input.\n");
#endif
  }
  else
  {
    ngchar1 = YYTRANSLATE(ngchar);

#if YYDEBUG != 0
    if (ngdebug)
    {
      fprintf (stderr, "Next token is %d (%s", ngchar, ngtname[ngchar1]);
      /* Give the individual parser a way to print the precise meaning
         of a token, for further debugging info.  */
#ifdef YYPRINT
      YYPRINT (stderr, ngchar, nglval);
#endif
      fprintf (stderr, ")\n");
    }
#endif
  }

  ngn += ngchar1;
  if (ngn < 0 || ngn > YYLAST || ngcheck[ngn] != ngchar1)
    goto ngdefault;

  ngn = ngtable[ngn];

  /* ngn is what to do for this token type in this state.
     Negative => reduce, -ngn is rule number.
     Positive => shift, ngn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (ngn < 0)
  {
    if (ngn == YYFLAG)
      goto ngerrlab;
    ngn = -ngn;
    goto ngreduce;
  }
  else if (ngn == 0)
    goto ngerrlab;

  if (ngn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YYDEBUG != 0
  if (ngdebug)
    fprintf(stderr, "Shifting token %d (%s), ", ngchar, ngtname[ngchar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (ngchar != YYEOF)
    ngchar = YYEMPTY;

  *++ngvsp = nglval;
#ifdef YYLSP_NEEDED
  *++nglsp = nglloc;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (ngerrstatus) ngerrstatus--;

  ngstate = ngn;
  goto ngnewstate;

  /* Do the default action for the current state.  */
ngdefault:

  ngn = ngdefact[ngstate];
  if (ngn == 0)
    goto ngerrlab;

  /* Do a reduction.  ngn is the number of a rule to reduce with.  */
ngreduce:
  nglen = ngr2[ngn];
  if (nglen > 0)
    ngval = ngvsp[1-nglen]; /* implement default value of the action */

#if YYDEBUG != 0
  if (ngdebug)
  {
    int i;

    fprintf (stderr, "Reducing via rule %d (line %d), ",
             ngn, ngrline[ngn]);

    /* Print the symbols being reduced, and their result.  */
    for (i = ngprhs[ngn]; ngrhs[i] > 0; i++)
      fprintf (stderr, "%s ", ngtname[ngrhs[i]]);
    fprintf (stderr, " -> %s\n", ngtname[ngr1[ngn]]);
  }
#endif


  switch (ngn) {

  case 5 :
#line 84 "ngin.y"
    {
      Elem.n_c=Elem.n_s=0;
      ;
      break;
    }
  case 6 :
#line 87 "ngin.y"
    {
      if (PutElement(&Elem)) YYABORT;
      ;
      break;
    }
  case 13 :
#line 100 "ngin.y"
    {
      Elem.subdom=(int)ngvsp[-3].ival;
      Elem.c_id[0]=(int)ngvsp[-2].ival;
      Elem.c_id[1]=(int)ngvsp[-1].ival;
      Elem.c_id[2]=(int)ngvsp[0].ival;
      Elem.n_c=3;
      ngval.el=&Elem;
      ;
      break;
    }
  case 14 :
#line 108 "ngin.y"
    {
      Elem.subdom=(int)ngvsp[-4].ival;
      Elem.c_id[0]=(int)ngvsp[-3].ival;
      Elem.c_id[1]=(int)ngvsp[-2].ival;
      Elem.c_id[2]=(int)ngvsp[-1].ival;
      Elem.c_id[3]=(int)ngvsp[0].ival;
      Elem.n_c=4;
      ngval.el=&Elem;
      ;
      break;
    }
  case 15 :
#line 119 "ngin.y"
    {
      Elem.side[Elem.n_s].c_id[0]=(int)ngvsp[-1].ival;
      Elem.side[Elem.n_s].c_id[1]=(int)ngvsp[0].ival;
      ngval.es=&(Elem.side[Elem.n_s]);
      Elem.n_s++;
      ;
      break;
    }
  case 18 :
#line 131 "ngin.y"
    {
      InnerNode.global[0]=ngvsp[-2].dval;
      InnerNode.global[1]=ngvsp[-1].dval;
      ngval.in=&InnerNode;
      PutInnerNode(&InnerNode);
      ;
      break;
    }
  case 21 :
#line 143 "ngin.y"
    {
      BndNode.n_lp=0;
      BndNode.global[0]=ngvsp[-1].dval;
      BndNode.global[1]=ngvsp[0].dval;
      ngval.bs=&BndNode;
      ;
      break;
    }
  case 22 :
#line 149 "ngin.y"
    {PutBndNode(&BndNode);;
     break;}
  case 24 :
#line 153 "ngin.y"
    {
      LP_COPY(&(BndNode.lp[BndNode.n_lp]),ngvsp[0].lp);
      BndNode.n_lp++;
      ngval.bs=&BndNode;
      ;
      break;
    }
  case 25 :
#line 158 "ngin.y"
    {
      LP_COPY(&(BndNode.lp[BndNode.n_lp]),ngvsp[0].lp);
      BndNode.n_lp++;
      ngval.bs=&BndNode;
      ;
      break;
    }
  case 26 :
#line 165 "ngin.y"
    {
      LinePos.line_id=(int)ngvsp[-1].ival;
      LinePos.local=(float)ngvsp[0].dval;
      ngval.lp=&LinePos;
      ;
      break;
    }
  case 27 :
#line 172 "ngin.y"
    {ngval.dval=(double)ngvsp[0].ival;;
     break;}
  case 29 :
#line 176 "ngin.y"
    {ngval.ival=ngvsp[0].ival;;
     break;}
  }
  /* the action file gets copied in in place of this dollarsign */
#line 498 "/usr/share/bison.simple"

  ngvsp -= nglen;
  ngssp -= nglen;
#ifdef YYLSP_NEEDED
  nglsp -= nglen;
#endif

#if YYDEBUG != 0
  if (ngdebug)
  {
    short *ssp1 = ngss - 1;
    fprintf (stderr, "state stack now");
    while (ssp1 != ngssp)
      fprintf (stderr, " %d", *++ssp1);
    fprintf (stderr, "\n");
  }
#endif

  *++ngvsp = ngval;

#ifdef YYLSP_NEEDED
  nglsp++;
  if (nglen == 0)
  {
    nglsp->first_line = nglloc.first_line;
    nglsp->first_column = nglloc.first_column;
    nglsp->last_line = (nglsp-1)->last_line;
    nglsp->last_column = (nglsp-1)->last_column;
    nglsp->text = 0;
  }
  else
  {
    nglsp->last_line = (nglsp+nglen-1)->last_line;
    nglsp->last_column = (nglsp+nglen-1)->last_column;
  }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  ngn = ngr1[ngn];

  ngstate = ngpgoto[ngn - YYNTBASE] + *ngssp;
  if (ngstate >= 0 && ngstate <= YYLAST && ngcheck[ngstate] == *ngssp)
    ngstate = ngtable[ngstate];
  else
    ngstate = ngdefgoto[ngn - YYNTBASE];

  goto ngnewstate;

ngerrlab:   /* here on detecting error */

  if (! ngerrstatus)
  /* If not already recovering from an error, report this error.  */
  {
    ++ngnerrs;

#ifdef YYERROR_VERBOSE
    ngn = ngpact[ngstate];

    if (ngn > YYFLAG && ngn < YYLAST)
    {
      int size = 0;
      char *msg;
      int x, count;

      count = 0;
      /* Start X at -ngn if nec to avoid negative indexes in ngcheck.  */
      for (x = (ngn < 0 ? -ngn : 0);
           x < (sizeof(ngtname) / sizeof(char *)); x++)
        if (ngcheck[x + ngn] == x)
          size += strlen(ngtname[x]) + 15, count++;
      msg = (char *) malloc(size + 15);
      if (msg != 0)
      {
        strcpy(msg, "parse error");

        if (count < 5)
        {
          count = 0;
          for (x = (ngn < 0 ? -ngn : 0);
               x < (sizeof(ngtname) / sizeof(char *)); x++)
            if (ngcheck[x + ngn] == x)
            {
              strcat(msg, count == 0 ? ", expecting `" : " or `");
              strcat(msg, ngtname[x]);
              strcat(msg, "'");
              count++;
            }
        }
        ngerror(msg);
        free(msg);
      }
      else
        ngerror ("parse error; also virtual memory exceeded");
    }
    else
#endif /* YYERROR_VERBOSE */
    ngerror("parse error");
  }

  goto ngerrlab1;
ngerrlab1:   /* here on error raised explicitly by an action */

  if (ngerrstatus == 3)
  {
    /* if just tried and failed to reuse lookahead token after an error, discard it.  */

    /* return failure if at end of input */
    if (ngchar == YYEOF)
      YYABORT;

#if YYDEBUG != 0
    if (ngdebug)
      fprintf(stderr, "Discarding token %d (%s).\n", ngchar, ngtname[ngchar1]);
#endif

    ngchar = YYEMPTY;
  }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  ngerrstatus = 3;              /* Each real token shifted decrements this */

  goto ngerrhandle;

ngerrdefault:  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  ngn = ngdefact[ngstate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (ngn) goto ngdefault;
#endif

ngerrpop:   /* pop the current state because it cannot handle the error token */

  if (ngssp == ngss) YYABORT;
  ngvsp--;
  ngstate = *--ngssp;
#ifdef YYLSP_NEEDED
  nglsp--;
#endif

#if YYDEBUG != 0
  if (ngdebug)
  {
    short *ssp1 = ngss - 1;
    fprintf (stderr, "Error: state stack now");
    while (ssp1 != ngssp)
      fprintf (stderr, " %d", *++ssp1);
    fprintf (stderr, "\n");
  }
#endif

ngerrhandle:

  ngn = ngpact[ngstate];
  if (ngn == YYFLAG)
    goto ngerrdefault;

  ngn += YYTERROR;
  if (ngn < 0 || ngn > YYLAST || ngcheck[ngn] != YYTERROR)
    goto ngerrdefault;

  ngn = ngtable[ngn];
  if (ngn < 0)
  {
    if (ngn == YYFLAG)
      goto ngerrpop;
    ngn = -ngn;
    goto ngreduce;
  }
  else if (ngn == 0)
    goto ngerrpop;

  if (ngn == YYFINAL)
    YYACCEPT;

#if YYDEBUG != 0
  if (ngdebug)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++ngvsp = nglval;
#ifdef YYLSP_NEEDED
  *++nglsp = nglloc;
#endif

  ngstate = ngn;
  goto ngnewstate;
}
#line 180 "ngin.y"



ngwrap (char *s)
{
  return (1);
}

ngerror (char *s)
{
  int line;
  char text[128];

  NP_Error(&line,text);
  NG_Print(ERROR_PREFIX "'%s', line %d\n",text,line);
  ngbreak();
}
