// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
typedef union
{
  /* put RCS string here in order to get it into yacc-generated header file
     static char RCS_ID("$Header: /hosts/dom/cvs/df/gen/problems/dfcfg/dfcfg.y,v 1.1
     1998/02/20 16:58:46 birken Exp $",UG_RCS_STRING);
   */


  /* transfer lex->yacc */
  double dval;
  long ival;

  LINE_POSITION *lp;
  SURFACE_POSITION *sp;
  BND_NODE *bs;
  INNER_NODE *in;
  ELEM_FACE *ef;
  NG_ELEMENT *el;
} YYSTYPE;
#define DOUBLE_VALUE    258
#define INT_VALUE       259
#define INODE   260
#define BNODE   261
#define SURFACE 262
#define LINE    263
#define ELEM    264
#define FACE    265
#define TEND    266


extern YYSTYPE nglval;
