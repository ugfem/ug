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




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 51 "ngin-yacc.y"
typedef union YYSTYPE {
  /* put RCS string here in order to get it into yacc-generated header file
     static char RCS_ID("$Header: /hosts/dom/cvs/df/gen/problems/dfcfg/dfcfg.y,v 1.1
     1998/02/20 16:58:46 birken Exp $",UG_RCS_STRING);
   */


  /* transfer lex->yacc */
  double dval;
  long ival;

  NS_DIM_PREFIX LINE_POSITION *lp;
  NS_DIM_PREFIX BND_NODE *bs;
  NS_DIM_PREFIX INNER_NODE *in;
  NS_DIM_PREFIX ELEM_SIDE *es;
  NS_DIM_PREFIX NG_ELEMENT *el;
} YYSTYPE;
/* Line 1240 of yacc.c.  */
#line 75 "y.tab.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE nglval;
