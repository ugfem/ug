// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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
#define DOUBLE_VALUE    258
#define INT_VALUE       259
#define INODE   260
#define BNODE   261
#define LINE    262
#define ELEM    263
#define SIDE    264
#define TEND    265


extern YYSTYPE nglval;
