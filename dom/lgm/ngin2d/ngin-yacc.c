// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* A Bison parser, made by GNU Bison 1.875a.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

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
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* If NAME_PREFIX is specified substitute the variables and functions
   names.  */
#define yyparse ngparse
#define yylex   nglex
#define yyerror ngerror
#define yylval  nglval
#define yychar  ngchar
#define yydebug ngdebug
#define yynerrs ngnerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
/* Put the tokens into the symbol table, so that GDB and other debuggers
   know about them.  */
enum yytokentype {
  DOUBLE_VALUE = 258,
  INT_VALUE = 259,
  INODE = 260,
  BNODE = 261,
  LINE = 262,
  ELEM = 263,
  SIDE = 264,
  TEND = 265
};
#endif
#define DOUBLE_VALUE 258
#define INT_VALUE 259
#define INODE 260
#define BNODE 261
#define LINE 262
#define ELEM 263
#define SIDE 264
#define TEND 265




/* Copy the first part of user declarations.  */
#line 1 "ngin-yacc.y"

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

#include "ng2d.h"

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





/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 44 "ngin-yacc.y"
typedef union YYSTYPE {
  /* put RCS string here in order to get it into yacc-generated header file
     static char RCS_ID("$Header: /hosts/dom/cvs/df/gen/problems/dfcfg/dfcfg.y,v 1.1
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
/* Line 191 of yacc.c.  */
#line 165 "ngin-yacc.c"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 177 "ngin-yacc.c"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
/* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
  && (! defined (__cplusplus) \
  || (YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
  ((N) * (sizeof (short) + sizeof (YYSTYPE))                         \
   + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
  __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)              \
  do                                        \
  {                                       \
    register YYSIZE_T yyi;                \
    for (yyi = 0; yyi < (Count); yyi++)   \
      (To)[yyi] = (From)[yyi];            \
  }                                       \
  while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)                                        \
  do                                                                  \
  {                                                                 \
    YYSIZE_T yynewbytes;                                            \
    YYCOPY (&yyptr->Stack, Stack, yysize);                          \
    Stack = &yyptr->Stack;                                          \
    yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
    yyptr += yynewbytes / sizeof (*yyptr);                          \
  }                                                                 \
  while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
typedef signed char yysigned_char;
#else
typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  8
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   38

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  11
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  19
/* YYNRULES -- Number of rules. */
#define YYNRULES  30
/* YYNRULES -- Number of states. */
#define YYNSTATES  50

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   265

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
  0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
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
  2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
  5,     6,     7,     8,     9,    10
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
  0,     0,     3,     7,    10,    12,    15,    16,    17,    23,
  25,    28,    32,    37,    43,    48,    54,    58,    60,    63,
  68,    70,    73,    74,    75,    83,    85,    88,    92,    94,
  96
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
  12,     0,    -1,    22,    20,    13,    -1,    22,    13,    -1,
  14,    -1,    13,    14,    -1,    -1,    -1,     8,    15,    17,
  16,    10,    -1,    18,    -1,    18,    19,    -1,    18,    19,
  19,    -1,    18,    19,    19,    19,    -1,    18,    19,    19,
  19,    19,    -1,    29,    29,    29,    29,    -1,    29,    29,
  29,    29,    29,    -1,     9,    29,    29,    -1,    21,    -1,
  20,    21,    -1,     5,    28,    28,    10,    -1,    23,    -1,
  22,    23,    -1,    -1,    -1,     6,    28,    28,    24,    26,
  25,    10,    -1,    27,    -1,    26,    27,    -1,     7,    29,
  28,    -1,     4,    -1,     3,    -1,     4,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned char yyrline[] =
{
  0,    76,    76,    77,    80,    81,    84,    87,    84,    93,
  94,    95,    96,    97,   100,   108,   119,   127,   128,   131,
  139,   140,   143,   149,   143,   153,   158,   165,   172,   173,
  176
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "DOUBLE_VALUE", "INT_VALUE", "INODE",
  "BNODE", "LINE", "ELEM", "SIDE", "TEND", "$accept", "Grid", "ElemList",
  "Elem", "@1", "@2", "ElemSpec", "InnerElemSpec", "ElemSide",
  "InnerNodeList", "InnerNode", "BndNodeList", "BndNode", "@3", "@4",
  "BndSpec", "LinePosition", "Coord", "Id", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
  0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
  265
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
  0,    11,    12,    12,    13,    13,    15,    16,    14,    17,
  17,    17,    17,    17,    18,    18,    19,    20,    20,    21,
  22,    22,    24,    25,    23,    26,    26,    27,    28,    28,
  29
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
  0,     2,     3,     2,     1,     2,     0,     0,     5,     1,
  2,     3,     4,     5,     4,     5,     3,     1,     2,     4,
  1,     2,     0,     0,     7,     1,     2,     3,     1,     1,
  1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
  0,     0,     0,     0,    20,    29,    28,     0,     1,     0,
  6,     3,     4,     0,    17,    21,    22,     0,     0,     5,
  2,    18,     0,     0,    30,     7,     9,     0,     0,    23,
  25,    19,     0,     0,    10,     0,     0,     0,    26,     8,
  0,    11,     0,    27,    24,    16,    12,    14,    13,    15
};

/* YYDEFGOTO[NTERM-NUM]. */
static const yysigned_char yydefgoto[] =
{
  -1,     2,    11,    12,    18,    32,    25,    26,    34,    13,
  14,     3,     4,    22,    37,    29,    30,     7,    27
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -34
static const yysigned_char yypact[] =
{
  11,     3,    21,    14,   -34,   -34,   -34,     3,   -34,     3,
  -34,    16,   -34,     7,   -34,   -34,   -34,     3,    22,   -34,
  16,   -34,    18,    17,   -34,   -34,    19,    22,    22,    18,
  -34,   -34,    20,    22,    19,    22,     3,    23,   -34,   -34,
  22,    19,    22,   -34,   -34,   -34,    19,    22,   -34,   -34
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
  -34,   -34,    24,    -6,   -34,   -34,   -34,   -34,   -33,   -34,
  25,   -34,    28,   -34,   -34,   -34,     5,    -7,   -24
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
  16,    41,    17,    35,    36,    19,     5,     6,    46,    40,
  23,    42,     9,    48,    19,    10,    45,     1,    47,     9,
  1,     8,    10,    49,    10,    28,    24,    31,    33,    43,
  39,    15,     0,    44,    38,     0,     0,    20,    21
};

static const yysigned_char yycheck[] =
{
  7,    34,     9,    27,    28,    11,     3,     4,    41,    33,
  17,    35,     5,    46,    20,     8,    40,     6,    42,     5,
  6,     0,     8,    47,     8,     7,     4,    10,     9,    36,
  10,     3,    -1,    10,    29,    -1,    -1,    13,    13
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
  0,     6,    12,    22,    23,     3,     4,    28,     0,     5,
  8,    13,    14,    20,    21,    23,    28,    28,    15,    14,
  13,    21,    24,    28,     4,    17,    18,    29,     7,    26,
  27,    10,    16,     9,    19,    29,    29,    25,    27,    10,
  29,    19,    29,    28,    10,    29,    19,    29,    19,    29
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrlab1


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL          goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
  do                                                              \
    if (yychar == YYEMPTY && yylen == 1)                          \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      yytoken = YYTRANSLATE (yychar);                           \
      YYPOPSTACK;                                               \
      goto yybackup;                                            \
    }                                                           \
    else                                                          \
    {                                                           \
      yyerror ("syntax error: cannot back up");\
      YYERROR;                                                  \
    }                                                           \
  while (0)

#define YYTERROR        1
#define YYERRCODE       256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)         \
  Current.first_line   = Rhs[1].first_line;      \
  Current.first_column = Rhs[1].first_column;    \
  Current.last_line    = Rhs[N].last_line;       \
  Current.last_column  = Rhs[N].last_column;
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
  do {                                            \
    if (yydebug)                                  \
      YYFPRINTF Args;                             \
  } while (0)

# define YYDSYMPRINT(Args)                      \
  do {                                            \
    if (yydebug)                                  \
      yysymprint Args;                            \
  } while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)            \
  do {                                                            \
    if (yydebug)                                                  \
    {                                                           \
      YYFPRINTF (stderr, "%s ", Title);                         \
      yysymprint (stderr,                                       \
                  Token, Value);        \
      YYFPRINTF (stderr, "\n");                                 \
    }                                                           \
  } while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (cinluded).                                                   |
   `------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
short *bottom;
short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
  do {                                                            \
    if (yydebug)                                                  \
      yy_stack_print ((Bottom), (Top));                           \
  } while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
   `------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
int yyrule;
#endif
{
  int yyi;
  unsigned int yylineno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylineno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)          \
  do {                                    \
    if (yydebug)                          \
      yy_reduce_print (Rule);             \
  } while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
char *yydest;
const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
   `--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
FILE *yyoutput;
int yytype;
YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
  {
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
  {
  default :
    break;
  }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
   `-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
int yytype;
YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
  {

  default :
    break;
  }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
   `----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{

  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;             /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

  /*------------------------------------------------------------.
  | yynewstate -- Push a new state, which is found in yystate.  |
     `------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
   */
  yyssp++;

yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
  {
    /* Get the current used size of the three stacks, in elements.  */
    YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
    {
      /* Give user a chance to reallocate the stack. Use copies of
         these so that the &'s don't force the real ones into
         memory.  */
      YYSTYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;


      /* Each stack pointer address is followed by the size of the
         data in use in that stack, in bytes.  This used to be a
         conditional around just the two extra args, but that might
         be undefined if yyoverflow is a macro.  */
      yyoverflow ("parser stack overflow",
                  &yyss1, yysize * sizeof (*yyssp),
                  &yyvs1, yysize * sizeof (*yyvsp),

                  &yystacksize);

      yyss = yyss1;
      yyvs = yyvs1;
    }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
    goto yyoverflowlab;
# else
    /* Extend the stack our own way.  */
    if (YYMAXDEPTH <= yystacksize)
      goto yyoverflowlab;
    yystacksize *= 2;
    if (YYMAXDEPTH < yystacksize)
      yystacksize = YYMAXDEPTH;

    {
      short *yyss1 = yyss;
      union yyalloc *yyptr =
        (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
      if (! yyptr)
        goto yyoverflowlab;
      YYSTACK_RELOCATE (yyss);
      YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
      if (yyss1 != yyssa)
        YYSTACK_FREE (yyss1);
    }
# endif
#endif /* no yyoverflow */

    yyssp = yyss + yysize - 1;
    yyvsp = yyvs + yysize - 1;


    YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                (unsigned long int) yystacksize));

    if (yyss + yystacksize - 1 <= yyssp)
      YYABORT;
  }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

  /*-----------.
  | yybackup.  |
     `-----------*/
yybackup:

  /* Do appropriate processing given the current state.  */
  /* Read a lookahead token if we need one and don't already have one.  */
  /* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
  {
    YYDPRINTF ((stderr, "Reading a token: "));
    yychar = YYLEX;
  }

  if (yychar <= YYEOF)
  {
    yychar = yytoken = YYEOF;
    YYDPRINTF ((stderr, "Now at end of input.\n"));
  }
  else
  {
    yytoken = YYTRANSLATE (yychar);
    YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
  }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
  {
    if (yyn == 0 || yyn == YYTABLE_NINF)
      goto yyerrlab;
    yyn = -yyn;
    goto yyreduce;
  }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
     `-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
     `-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
  {
  case 6 :
#line 84 "ngin-yacc.y"
    {
      Elem.n_c=Elem.n_s=0;
      ;
    }
    break;

  case 7 :
#line 87 "ngin-yacc.y"
    {
      if (PutElement(&Elem)) YYABORT;
      ;
    }
    break;

  case 14 :
#line 100 "ngin-yacc.y"
    {
      Elem.subdom=(int)yyvsp[-3].ival;
      Elem.c_id[0]=(int)yyvsp[-2].ival;
      Elem.c_id[1]=(int)yyvsp[-1].ival;
      Elem.c_id[2]=(int)yyvsp[0].ival;
      Elem.n_c=3;
      yyval.el=&Elem;
      ;
    }
    break;

  case 15 :
#line 108 "ngin-yacc.y"
    {
      Elem.subdom=(int)yyvsp[-4].ival;
      Elem.c_id[0]=(int)yyvsp[-3].ival;
      Elem.c_id[1]=(int)yyvsp[-2].ival;
      Elem.c_id[2]=(int)yyvsp[-1].ival;
      Elem.c_id[3]=(int)yyvsp[0].ival;
      Elem.n_c=4;
      yyval.el=&Elem;
      ;
    }
    break;

  case 16 :
#line 119 "ngin-yacc.y"
    {
      Elem.side[Elem.n_s].c_id[0]=(int)yyvsp[-1].ival;
      Elem.side[Elem.n_s].c_id[1]=(int)yyvsp[0].ival;
      yyval.es=&(Elem.side[Elem.n_s]);
      Elem.n_s++;
      ;
    }
    break;

  case 19 :
#line 131 "ngin-yacc.y"
    {
      InnerNode.global[0]=yyvsp[-2].dval;
      InnerNode.global[1]=yyvsp[-1].dval;
      yyval.in=&InnerNode;
      PutInnerNode(&InnerNode);
      ;
    }
    break;

  case 22 :
#line 143 "ngin-yacc.y"
    {
      BndNode.n_lp=0;
      BndNode.global[0]=yyvsp[-1].dval;
      BndNode.global[1]=yyvsp[0].dval;
      yyval.bs=&BndNode;
      ;
    }
    break;

  case 23 :
#line 149 "ngin-yacc.y"
    {PutBndNode(&BndNode);;}
    break;

  case 25 :
#line 153 "ngin-yacc.y"
    {
      LP_COPY(&(BndNode.lp[BndNode.n_lp]),yyvsp[0].lp);
      BndNode.n_lp++;
      yyval.bs=&BndNode;
      ;
    }
    break;

  case 26 :
#line 158 "ngin-yacc.y"
    {
      LP_COPY(&(BndNode.lp[BndNode.n_lp]),yyvsp[0].lp);
      BndNode.n_lp++;
      yyval.bs=&BndNode;
      ;
    }
    break;

  case 27 :
#line 165 "ngin-yacc.y"
    {
      LinePos.line_id=(int)yyvsp[-1].ival;
      LinePos.local=(float)yyvsp[0].dval;
      yyval.lp=&LinePos;
      ;
    }
    break;

  case 28 :
#line 172 "ngin-yacc.y"
    {yyval.dval=(double)yyvsp[0].ival;;}
    break;

  case 30 :
#line 176 "ngin-yacc.y"
    {yyval.ival=yyvsp[0].ival;;}
    break;


  }

  /* Line 999 of yacc.c.  */
#line 1195 "ngin-yacc.c"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


  /*------------------------------------.
  | yyerrlab -- here on detecting error |
     `------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
  {
    ++yynerrs;
#if YYERROR_VERBOSE
    yyn = yypact[yystate];

    if (YYPACT_NINF < yyn && yyn < YYLAST)
    {
      YYSIZE_T yysize = 0;
      int yytype = YYTRANSLATE (yychar);
      char *yymsg;
      int yyx, yycount;

      yycount = 0;
      /* Start YYX at -YYN if negative to avoid negative indexes in
         YYCHECK.  */
      for (yyx = yyn < 0 ? -yyn : 0;
           yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
        if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
          yysize += yystrlen (yytname[yyx]) + 15, yycount++;
      yysize += yystrlen ("syntax error, unexpected ") + 1;
      yysize += yystrlen (yytname[yytype]);
      yymsg = (char *) YYSTACK_ALLOC (yysize);
      if (yymsg != 0)
      {
        char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
        yyp = yystpcpy (yyp, yytname[yytype]);

        if (yycount < 5)
        {
          yycount = 0;
          for (yyx = yyn < 0 ? -yyn : 0;
               yyx < (int) (sizeof (yytname) / sizeof (char *));
               yyx++)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
            {
              const char *yyq = ! yycount ? ", expecting " : " or ";
              yyp = yystpcpy (yyp, yyq);
              yyp = yystpcpy (yyp, yytname[yyx]);
              yycount++;
            }
        }
        yyerror (yymsg);
        YYSTACK_FREE (yymsg);
      }
      else
        yyerror ("syntax error; also virtual memory exhausted");
    }
    else
#endif /* YYERROR_VERBOSE */
    yyerror ("syntax error");
  }



  if (yyerrstatus == 3)
  {
    /* If just tried and failed to reuse lookahead token after an
       error, discard it.  */

    /* Return failure if at end of input.  */
    if (yychar == YYEOF)
    {
      /* Pop the error token.  */
      YYPOPSTACK;
      /* Pop the rest of the stack.  */
      while (yyss < yyssp)
      {
        YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
        yydestruct (yystos[*yyssp], yyvsp);
        YYPOPSTACK;
      }
      YYABORT;
    }

    YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
    yydestruct (yytoken, &yylval);
    yychar = YYEMPTY;

  }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


  /*----------------------------------------------------.
  | yyerrlab1 -- error raised explicitly by an action.  |
     `----------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
  {
    yyn = yypact[yystate];
    if (yyn != YYPACT_NINF)
    {
      yyn += YYTERROR;
      if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
      {
        yyn = yytable[yyn];
        if (0 < yyn)
          break;
      }
    }

    /* Pop the current state because it cannot handle the error token.  */
    if (yyssp == yyss)
      YYABORT;

    YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
    yydestruct (yystos[yystate], yyvsp);
    yyvsp--;
    yystate = *--yyssp;

    YY_STACK_PRINT (yyss, yyssp);
  }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


  /*-------------------------------------.
  | yyacceptlab -- YYACCEPT comes here.  |
     `-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

  /*-----------------------------------.
  | yyabortlab -- YYABORT comes here.  |
     `-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
  /*----------------------------------------------.
  | yyoverflowlab -- parser overflow comes here.  |
     `----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 180 "ngin-yacc.y"



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
