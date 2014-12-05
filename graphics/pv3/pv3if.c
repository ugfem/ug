// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pv3if.c                                                       */
/*                                                                          */
/* Purpose:   Interface to pV3                                              */
/*                                                                          */
/* Author:    Michael Lampe                                                 */
/*            IWR Uni Heidelberg                                            */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            Email: Michael.Lampe@iwr.uni-heidelberg.de                    */
/*                                                                          */
/* History:   20030130 begin, ug3.8                                         */
/*                                                                          */
/* Remarks:   Don't despair. pV3 is written in Fortran 77 ;-)               */
/*                                                                          */
/*            TODO: 1) support for on-line visualization                    */
/*                  2) give pV3 info about partition boundaries             */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include "config.h"
#include <string.h>
#include <pV3.h>

#include "general.h"
#include "parallel.h"
#include "gm.h"
#include "shapes.h"
#include "commands.h"
#include "pv3if.h"

/****************************************************************************/
/*																			*/
/* defines in the following order										        */
/*																			*/
/*		compile time constants defining static data size (i.e. arrays)		*/
/*		other constants												                */
/*		macros																*/
/*																			*/
/****************************************************************************/

#define MAX_ITEM    50

#undef MIN
#define MIN(a, b)   ((a)<=(b) ? (a) : (b))

#define SURFACE_LOOP_BEGIN(mg, e) \
  { \
    int _i; \
    for (_i = 0; _i <= TOPLEVEL(mg); _i++) \
      for (e = FIRSTELEMENT(GRID_ON_LEVEL(mg, _i)); e != NULL; e = SUCCE(e)) { \
        if (!EstimateHere(e)) continue;

#define SURFACE_LOOP_END \
  } \
  }

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* initialization data */
static char title[]="UG Output";     /* title string                        */
static int cid;                      /* client id                           */
static char cname[]=" ";             /* client name                         */
static char dname[]=" ";             /* not needed                          */
static int iopt;                     /* unsteady control parameter          */
static int npgcut=0;                 /* # of programmer defined cuts        */
static char tpgcut[MAX_ITEM][32];    /* title for each cut                  */
static int nkeys=0;                  /* # of active keyboard keys           */
static int ikeys[MAX_ITEM];          /* X-keypress return code for each key */
static char tkeys[MAX_ITEM][32];     /* title for each key                  */
static int fkeys[MAX_ITEM];          /* type of function controlled         */
static float flims[MAX_ITEM][2];     /* function limits/scales              */
static int mirror=0;                 /* not needed                          */
static float repmat[16];             /* not needed                          */
static int maxblk=0;                 /* not needed                          */
static int istat;                    /* startup/terminate state             */

/* evaluation data */
static union {
  EVALUES *s;
  EVECTOR *v;
} eval[MAX_ITEM];
static char eval_name[MAX_ITEM][NAMESIZE];

/****************************************************************************/
/*                                                                          */
/*    utility functions                                                     */
/*                                                                          */
/****************************************************************************/

static void fstring(char *d, const char *s, size_t n)
{
  int i;

  strncpy(d, s, n);
  for(i = MIN(strlen(s), n); i < n; i++)
    d[i]=' ';
}

static void ClearVertexMarkers(MULTIGRID *mg)
{
  VERTEX *v;
  int i;

  for (i = 0; i <= TOPLEVEL(mg); i++)
    for (v = FIRSTVERTEX(GRID_ON_LEVEL(mg, i)); v != NULL; v = SUCCV(v))
      SETUSED(v,0);
}

/****************************************************************************/
/*                                                                          */
/*    callback for grid statistics                                          */
/*                                                                          */
/****************************************************************************/

void pVSTRUC(int *knode, int *kequiv, int *kcel1, int *kcel2, int *kcel3,
             int *kcel4, int *knptet, int *kptet, int *knblock, int *blocks,
             int *kphedra, int *ksurf, int *knsurf, int *hint)
{
  MULTIGRID *mg;
  ELEMENT *e;
  VERTEX *v;
  int i;

  *knode = *kcel1 = *kcel2 = *kcel3 = *kcel4 = *ksurf = 0;
  mg = GetCurrentMultigrid();
  ClearVertexMarkers(mg);
  SURFACE_LOOP_BEGIN(mg, e)
  switch (TAG(e))
  {
  case TETRAHEDRON :
    (*kcel1)++;
    break;
  case PYRAMID :
    (*kcel2)++;
    break;
  case PRISM :
    (*kcel3)++;
    break;
  case HEXAHEDRON :
    (*kcel4)++;
  }
  for (i = 0; i < CORNERS_OF_ELEM(e); i++) {
    v = MYVERTEX(CORNER(e, i));
    if (USED(v)) continue;
    SETUSED(v, 1);
    ID(v) = *knode;                           /* number vertices */
    (*knode)++;
  }
  /* check for domain boundary sides */
  if (OBJT(e) == BEOBJ) {
    for (i = 0; i < SIDES_OF_ELEM(e); i++)
      if (SIDE_ON_BND(e, i) && !InnerBoundary(e, i))
        (*ksurf)++;
  }
  SURFACE_LOOP_END
  *kequiv  = 0;
  *knptet  = 0;
  *kptet   = 0;
  *knblock = 0;
  *kphedra = 0;
  *knsurf  = 1;
}

/****************************************************************************/
/*                                                                          */
/*    callback for element vertices                                         */
/*                                                                          */
/****************************************************************************/

void pVCELL(int *cel1, int *cel2, int *cel3, int *cel4, int *nptet, int *ptet)
{
  MULTIGRID *mg;
  ELEMENT *e;

  mg = GetCurrentMultigrid();
  SURFACE_LOOP_BEGIN(mg, e)
  switch (TAG(e))
  {
  case TETRAHEDRON :
    cel1[0] = ID(MYVERTEX(CORNER(e, 3)))+1;
    cel1[1] = ID(MYVERTEX(CORNER(e, 0)))+1;
    cel1[2] = ID(MYVERTEX(CORNER(e, 1)))+1;
    cel1[3] = ID(MYVERTEX(CORNER(e, 2)))+1;
    cel1 += 4;
    break;
  case PYRAMID :
    cel2[0] = ID(MYVERTEX(CORNER(e, 0)))+1;
    cel2[1] = ID(MYVERTEX(CORNER(e, 1)))+1;
    cel2[2] = ID(MYVERTEX(CORNER(e, 2)))+1;
    cel2[3] = ID(MYVERTEX(CORNER(e, 3)))+1;
    cel2[4] = ID(MYVERTEX(CORNER(e, 4)))+1;
    cel2 += 5;
    break;
  case PRISM :
    cel3[0] = ID(MYVERTEX(CORNER(e, 1)))+1;
    cel3[1] = ID(MYVERTEX(CORNER(e, 4)))+1;
    cel3[2] = ID(MYVERTEX(CORNER(e, 5)))+1;
    cel3[3] = ID(MYVERTEX(CORNER(e, 2)))+1;
    cel3[4] = ID(MYVERTEX(CORNER(e, 3)))+1;
    cel3[5] = ID(MYVERTEX(CORNER(e, 0)))+1;
    cel3 += 6;
    break;
  case HEXAHEDRON :
    cel4[0] = ID(MYVERTEX(CORNER(e, 0)))+1;
    cel4[1] = ID(MYVERTEX(CORNER(e, 1)))+1;
    cel4[2] = ID(MYVERTEX(CORNER(e, 2)))+1;
    cel4[3] = ID(MYVERTEX(CORNER(e, 3)))+1;
    cel4[4] = ID(MYVERTEX(CORNER(e, 4)))+1;
    cel4[5] = ID(MYVERTEX(CORNER(e, 5)))+1;
    cel4[6] = ID(MYVERTEX(CORNER(e, 6)))+1;
    cel4[7] = ID(MYVERTEX(CORNER(e, 7)))+1;
    cel4 += 8;
    break;
  }
  SURFACE_LOOP_END
}

/****************************************************************************/
/*                                                                          */
/*    callback for vertex coordinates                                       */
/*                                                                          */
/****************************************************************************/

void pVGRID(float *xyz)
{
  MULTIGRID *mg;
  ELEMENT *e;
  VERTEX *v;
  int i;

  mg = GetCurrentMultigrid();
  ClearVertexMarkers(mg);
  SURFACE_LOOP_BEGIN(mg, e)
  for (i = 0; i < CORNERS_OF_ELEM(e); i++) {
    v = MYVERTEX(CORNER(e, i));
    if (USED(v)) continue;
    SETUSED(v, 1);
    *xyz++ = XC(v);
    *xyz++ = YC(v);
    *xyz++ = ZC(v);
  }
  SURFACE_LOOP_END
}

/****************************************************************************/
/*                                                                          */
/*    callback for surface data                                             */
/*                                                                          */
/****************************************************************************/

void pVSURFACE(int *nsurf, int *scon, int *scel, char *tsurf, int tsurf_len)
{
  MULTIGRID *mg;
  ELEMENT *e;
  VERTEX *v;
  int i, j, n;

  mg = GetCurrentMultigrid();

  /* domain surface */
  n = 0;
  SURFACE_LOOP_BEGIN(mg, e)
  if (OBJT(e) == BEOBJ) {
    for (i = 0; i < SIDES_OF_ELEM(e); i++)
      if (SIDE_ON_BND(e, i) && !InnerBoundary(e, i)) {
        *scon++ = 0;
        scel[3] = 0;
        for (j = 0; j < CORNERS_OF_SIDE(e, i); j++)
          scel[j] = ID(MYVERTEX(CORNER(e, CORNER_OF_SIDE(e, i, j))))+1;
        scel += 4;
        n++;
      }
  }
  SURFACE_LOOP_END
    nsurf[0] = n;
  nsurf[1] = 2;
  nsurf[2] = 1;
  fstring(tsurf, "Domain Surface", 20);
}

/****************************************************************************/
/*                                                                          */
/*    callback for scalar data                                              */
/*                                                                          */
/****************************************************************************/

void pVSCAL(int *key, float *s)
{
  MULTIGRID *mg;
  ELEMENT *e;
  VERTEX *v;
  PreprocessingProcPtr pre;
  ElementEvalProcPtr foo;
  double *cc[MAX_CORNERS_OF_ELEM], lc[3];
  int i, k;

  k = (*key)-1;
  mg = GetCurrentMultigrid();
  pre = eval[k].s->PreprocessProc;
  if (pre != NULL) pre(eval_name[k], mg);
  foo = eval[k].s->EvalProc;
  mg = GetCurrentMultigrid();
  ClearVertexMarkers(mg);
  SURFACE_LOOP_BEGIN(mg, e)
  for (i = 0; i < CORNERS_OF_ELEM(e); i++)
    cc[i] = CVECT(MYVERTEX(CORNER(e, i)));
  for (i = 0; i < CORNERS_OF_ELEM(e); i++) {
    v = MYVERTEX(CORNER(e, i));
    if (USED(v)) continue;
    SETUSED(v, 1);
    LocalCornerCoordinates(3, TAG(e), i, lc);
    *s++ = (float)foo(e, cc, lc);
  }
  SURFACE_LOOP_END
}

/****************************************************************************/
/*                                                                          */
/*    callback for vector data                                              */
/*                                                                          */
/****************************************************************************/

void pVVECT(int *key, float *V)
{
  MULTIGRID *mg;
  ELEMENT *e;
  VERTEX *v;
  PreprocessingProcPtr pre;
  ElementVectorProcPtr foo;
  double *cc[MAX_CORNERS_OF_ELEM], lc[3], vv[3];
  int i, k;

  k = (*key)-1;
  mg = GetCurrentMultigrid();
  pre = eval[k].v->PreprocessProc;
  if (pre != NULL) pre(eval_name[k], mg);
  foo = eval[k].v->EvalProc;
  mg = GetCurrentMultigrid();
  ClearVertexMarkers(mg);
  SURFACE_LOOP_BEGIN(mg, e)
  for (i = 0; i < CORNERS_OF_ELEM(e); i++)
    cc[i] = CVECT(MYVERTEX(CORNER(e, i)));
  for (i = 0; i < CORNERS_OF_ELEM(e); i++) {
    v = MYVERTEX(CORNER(e, i));
    if (USED(v)) continue;
    SETUSED(v, 1);
    LocalCornerCoordinates(3, TAG(e), i, lc);
    foo(e, cc, lc, vv);
    *V++ = (float)vv[0];
    *V++ = (float)vv[1];
    *V++ = (float)vv[2];
  }
  SURFACE_LOOP_END
}

/****************************************************************************/
/*                                                                          */
/*   pv3scal $k <key> $T <title> $ns <eval> $s <scal> $f <from> $t <to>      */
/*                                                                          */
/****************************************************************************/

static INT pv3ScalCommand(INT argc, char **argv)
{
  char k, T[32], ns[NAMELEN], s[NAMELEN];
  float f, t;

  if (strncmp(argv[1],"k",1) != 0) return PARAMERRORCODE;
  sscanf(argv[1], "k %c", &k);
  if (strncmp(argv[2],"T",1) != 0) return PARAMERRORCODE;
  sscanf(argv[2], "T %s", T);
  if (strncmp(argv[3],"ns",2) != 0) return PARAMERRORCODE;
  sscanf(argv[3], "ns %s", ns);
  if (strncmp(argv[4],"s",1) != 0) return PARAMERRORCODE;
  sscanf(argv[4], "s %s", s);
  if (strncmp(argv[5],"f",1) != 0) return PARAMERRORCODE;
  sscanf(argv[5], "f %f", &f);
  if (strncmp(argv[6],"t",1) != 0) return PARAMERRORCODE;
  sscanf(argv[6], "t %f", &t);
  ikeys[nkeys] = k;
  fstring(tkeys[nkeys], T, 32);
  fkeys[nkeys] = 1;
  flims[nkeys][0] = f;
  flims[nkeys][1] = t;
  eval[nkeys].s = GetElementValueEvalProc(ns);
  strcpy(eval_name[nkeys], s);
  nkeys++;
  return OKCODE;
}

/****************************************************************************/
/*                                                                          */
/*   pv3vect $k <key> $T <title> $nv <eval> $s <vect> $f <scale>            */
/*                                                                          */
/****************************************************************************/

static INT pv3VectCommand(INT argc, char **argv)
{
  char k, T[32], nv[NAMELEN], s[NAMELEN];
  float f;

  if (strncmp(argv[1],"k",1) != 0) return PARAMERRORCODE;
  sscanf(argv[1], "k %c", &k);
  if (strncmp(argv[2],"T",1) != 0) return PARAMERRORCODE;
  sscanf(argv[2], "T %s", T);
  if (strncmp(argv[3],"nv",2) != 0) return PARAMERRORCODE;
  sscanf(argv[3], "nv %s", nv);
  if (strncmp(argv[4],"s",1) != 0) return PARAMERRORCODE;
  sscanf(argv[4], "s %s", s);
  if (strncmp(argv[5],"f",1) != 0) return PARAMERRORCODE;
  sscanf(argv[5], "f %f", &f);
  ikeys[nkeys] = k;
  fstring(tkeys[nkeys], T, 32);
  fkeys[nkeys] = 2;
  flims[nkeys][0] = f;
  flims[nkeys][1] = f;
  eval[nkeys].v = GetElementVectorEvalProc(nv);
  strcpy(eval_name[nkeys], s);
  nkeys++;
  return OKCODE;
}

/****************************************************************************/
/*                                                                          */
/*   pv3plot                                                                */
/*                                                                          */
/****************************************************************************/

static INT pv3PlotCommand(INT argc, char **argv)
{
  iopt  = 0;       /* steady grid & data */
  istat = 3;       /* wait for and terminate with pV3 server */
  cid   = me+1;
  pV_INIT(title, &cid, cname, dname, &iopt, &npgcut, tpgcut[0],
          &nkeys, ikeys, tkeys[0], fkeys, flims[0], &mirror, repmat,
          &maxblk, &istat, strlen(title), strlen(cname), strlen(dname),
          32, 32);
  if (istat != 0) return CMDERRORCODE;
  return OKCODE;
}

/****************************************************************************/

INT InitPV3(void)
{
  if (CreateCommand("pv3scal", pv3ScalCommand) == NULL) return __LINE__;
  if (CreateCommand("pv3vect", pv3VectCommand) == NULL) return __LINE__;
  if (CreateCommand("pv3plot", pv3PlotCommand) == NULL) return __LINE__;
  return 0;
}
