// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      ug-famg.c														*/
/*																			*/
/* Purpose:   ug - famg interface											*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   Aug 1997 begin, ug version 3.7								*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <math.h>

/* ug library */
#include "gm.h"        /* for data structure               */
#include "evm.h"       /* for data structure               */
#include "ugdevices.h" /* for UserWrite, PrintErrorMessage */
#include "np.h"        /* for CreateNumProc                */
#include "debug.h"
#include "ugstruct.h"
#include "iter.h"

#include "famginterface.h"

/* test */
#include "wpm.h"
#include "wop.h"
#include "connectuggrape.h"
#include "uginterface.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/
#define DISPLAY_NP_FORMAT_SE                    "%-16.13s = %-7.4E\n"


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/
REP_ERR_FILE;

static struct CMG_Interface cmg_interface;
static struct CMG_Parameter cmg_parameter;

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/*  Class Definition                                                                    */
/*                                                                          */
/****************************************************************************/

typedef struct
{
  NP_ITER iter;

  INT heap;
  INT n1;
  INT n2;
  INT gamma;
  INT cgnodes;
  DOUBLE coarsening;
  DOUBLE strong;
  INT adaptive;
  INT maxit;
  DOUBLE alimit;
  DOUBLE rlimit;
  DOUBLE divlimit;
  DOUBLE reduction;
  INT famg_mark_key;
} NP_CMG;

/****************************************************************************/
/*                                                                          */
/* functions                                                                                            */
/*                                                                          */
/****************************************************************************/


static INT CMGInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_CMG *np;
  char *str;

  np = (NP_CMG *) theNP;

  if (ReadArgvINT("h",&(cmg_parameter.heap),argc,argv))
    cmg_parameter.heap = 1e+7;
  if (ReadArgvINT("n1",&(cmg_parameter.n1),argc,argv))
    cmg_parameter.n1 = 1;
  if (ReadArgvINT("n2",&(cmg_parameter.n2),argc,argv))
    cmg_parameter.n2 = 1;
  if (ReadArgvINT("g",&(cmg_parameter.gamma),argc,argv))
    cmg_parameter.gamma = 1;
  if (ReadArgvINT("cgn",&(cmg_parameter.cgnodes),argc,argv))
    cmg_parameter.cgnodes = 1;


  cmg_parameter.ilut = 1e+10;
  GetStringValueDouble(":cmg:ilut",&(cmg_parameter.ilut));

  cmg_parameter.cgilut = 0.0;
  GetStringValueDouble(":cmg:cgilut",&(cmg_parameter.cgilut));

  cmg_parameter.conloops = 0;
  GetStringValueInt(":cmg:conloops",&(cmg_parameter.conloops));

  cmg_parameter.mincoarse = 0.8;
  GetStringValueDouble(":cmg:mincoarse",&(cmg_parameter.mincoarse));

  cmg_parameter.type = 0;
  GetStringValueInt(":cmg:type",&(cmg_parameter.type));

  cmg_parameter.stv = 0;
  GetStringValueInt(":cmg:stv",&(cmg_parameter.stv));

  cmg_parameter.tol = 0.95;
  GetStringValueDouble(":cmg:tol",&(cmg_parameter.tol));

  cmg_parameter.sigma = 0.45;
  GetStringValueDouble(":cmg:sigma",&(cmg_parameter.sigma));

  cmg_parameter.omegar = 1.0;
  GetStringValueDouble(":cmg:omegar",&(cmg_parameter.omegar));

  cmg_parameter.omegal = 1.0;
  GetStringValueDouble(":cmg:omegal",&(cmg_parameter.omegal));

  cmg_parameter.error1 = 1e-6;
  GetStringValueDouble(":cmg:error1",&(cmg_parameter.error1));

  cmg_parameter.error2 = 1.0;
  GetStringValueDouble(":cmg:error2",&(cmg_parameter.error2));

  cmg_parameter.maxit = 100;
  GetStringValueInt(":cmg:maxit",&(cmg_parameter.maxit));

  cmg_parameter.alimit = 1e-14;
  GetStringValueDouble(":cmg:alimit",&(cmg_parameter.alimit));

  cmg_parameter.rlimit = 1e-10;
  GetStringValueDouble(":cmg:rlimit",&(cmg_parameter.rlimit));

  cmg_parameter.divlimit= 10.0;
  GetStringValueDouble(":cmg:divlimit",&(cmg_parameter.divlimit));

  cmg_parameter.reduction = 1.0;
  GetStringValueDouble(":cmg:reduction",&(cmg_parameter.reduction));

  strcpy(cmg_parameter.solver,"linit");
  str = GetStringVar(":cmg:solver");
  if(str != NULL) strcpy(cmg_parameter.solver,str);

  strcpy(cmg_parameter.presmoother,"fgs");
  str = GetStringVar(":cmg:presmoother");
  if(str != NULL) strcpy(cmg_parameter.presmoother,str);

  strcpy(cmg_parameter.postsmoother,"bgs");
  str = GetStringVar(":cmg:postsmoother");
  if(str != NULL) strcpy(cmg_parameter.postsmoother,str);

  strcpy(cmg_parameter.cgsmoother,"ilut");
  str = GetStringVar(":cmg:cgsmoother");
  if(str != NULL) strcpy(cmg_parameter.cgsmoother,str);


  np->heap = cmg_parameter.heap;
  np->n1 = cmg_parameter.n1;
  np->n2 = cmg_parameter.n2;
  np->gamma = cmg_parameter.gamma;
  np->cgnodes = cmg_parameter.cgnodes;
  np->maxit = cmg_parameter.maxit;
  np->alimit = cmg_parameter.alimit;
  np->rlimit = cmg_parameter.rlimit;
  np->divlimit = cmg_parameter.divlimit;
  np->reduction = cmg_parameter.reduction;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT CMGDisplay (NP_BASE *theNP)
{
  NP_CMG *np;

  np = (NP_CMG *) theNP;

  NPIterDisplay(&np->iter);

  UserWrite("configuration parameters:\n");
  UserWriteF(DISPLAY_NP_FORMAT_SI,"h",(int)np->heap);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n1",(int)np->n1);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n2",(int)np->n2);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"g",(int)np->gamma);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"cgn",(int)np->cgnodes);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"maxit",(int)np->maxit);
  UserWriteF(DISPLAY_NP_FORMAT_SE,"alimit",(double)np->alimit);
  UserWriteF(DISPLAY_NP_FORMAT_SE,"rlimit",(double)np->rlimit);
  UserWriteF(DISPLAY_NP_FORMAT_SE,"divlim",(double)np->divlimit);
  UserWriteF(DISPLAY_NP_FORMAT_SE,"red",(double)np->reduction);

  return (0);
}

static INT CMGPreProcess  (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           INT *baselevel, INT *result)
{
  MULTIGRID *mg;
  GRID *grid;
  NODE *node;
  VECTOR *vec,*w,*fw;
  MATRIX *m;
  SHORT xc,bc,mc,xmask,bmask;
  INT i,j,lev,found,ll, nnb, offset;
  DOUBLE d, sum;

  int n, nl, nv, *index, *start;
  double *entry, *vector[CMG_NVECTORS];
  void **extra;
  NP_CMG *np;

  np = (NP_CMG *) theNP;

  mg = NP_MG(theNP);
  MarkTmpMem(MGHEAP(mg),&np->famg_mark_key);   /* release in PostProcess */


  if (MD_IS_SCALAR(A) && VD_IS_SCALAR(x) && VD_IS_SCALAR(b))
  {
    xc    = VD_SCALCMP(x);
    mc    = MD_SCALCMP(A);
    bc    = VD_SCALCMP(b);
    xmask  = VD_SCALTYPEMASK(x);
    bmask  = VD_SCALTYPEMASK(b);
  }
  else
  {
    UserWrite("Not a scalar equation. \n");
    REP_ERR_RETURN(1);
  }


  /* count unknowns */
  n = 0;
  for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
  {
    for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
    {
      if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (FINE_GRID_DOF(vec)))
      {
        VINDEX(vec) = n;
        n++;
      }
    }
  }

  grid =  GRID_ON_LEVEL(mg,level);
  for (vec=FIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
  {
    if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (NEW_DEFECT(vec)))
    {
      VINDEX(vec) = n;
      n++;
    }
  }


  /* ug node information */
  extra = (void **) GetTmpMem(MGHEAP(mg),n*sizeof(void *),np->famg_mark_key);
  if (extra == NULL)
  {
    UserWrite("CMGCreateSystem: not enough memory. \n");
    REP_ERR_RETURN(1);
  }

  i = 0;
  for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
  {
    for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
    {
      if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (FINE_GRID_DOF(vec)))
      {
        extra[i] = (void *) MYVERTEX(VMYNODE(vec));
        i++;
      }
    }
  }

  grid =  GRID_ON_LEVEL(mg,level);
  for (vec=FIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
  {
    if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&xmask) && (NEW_DEFECT(vec)))
    {
      extra[i] = (void *) MYVERTEX(VMYNODE(vec));
      i++;
    }
  }

  /* allocate row/column start array */
  start = (int *) GetTmpMem(MGHEAP(mg),(n+1)*sizeof(int),np->famg_mark_key);
  if (start == NULL)
  {
    UserWrite("ug - cmg: not enough memory. \n");
    REP_ERR_RETURN(1);
  }

  /* allocate vectors */
  for(j = 0; j < CMG_NVECTORS; j++)
  {
    vector[j] = (DOUBLE *) GetTmpMem(MGHEAP(mg),n*sizeof(DOUBLE),np->famg_mark_key);
    if (vector[j] == NULL)
    {
      UserWrite("ug - cmg: not enough memory. \n");
      REP_ERR_RETURN(1);
    }
  }

  /* copy UG system matrix  into interface data structure */
  /* first step: count links */
  i = 0;
  start[0] = 0;
  for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
  {
    for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
    {
      if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&bmask) && (FINE_GRID_DOF(vec)))
      {
        nnb = 1;
        for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
        {
          w = MDEST(m);
          found = 0;
          if(FINE_GRID_DOF(w))
          {
            if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
            {
              nnb++;
            }
            found = 1;
          }
          if(!found)
          {
            node = VMYNODE(w);
            while(CORNERTYPE(node))
            {
              node = (NODE *)NFATHER(node);
              fw = NVECTOR(node);
              if(FINE_GRID_DOF(fw))
              {
                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                {
                  nnb++;
                }
                found = 1;
                break;
              }
            }
          }
          if(!found)
          {
            node = (NODE *)SONNODE(VMYNODE(w));
            ll = lev; ll++;
            while(node != NULL)
            {
              fw = NVECTOR(node);
              if((FINE_GRID_DOF(fw) && (ll < level))
                 || (NEW_DEFECT(fw) && (ll == level)))
              {
                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                {
                  nnb++;
                }
                found = 1;
                break;
              }
              node = (NODE *)SONNODE(node);
              ll++;
            }
          }
          if(!found)
          {
            UserWrite("error in CMGSolve \n");
          }
        }
        i++;
        start[i] = nnb+start[i-1];
      }
    }
  }
  grid =  GRID_ON_LEVEL(mg,level);
  for (vec=FIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
  {
    if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
    {
      nnb = 1;
      for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
      {
        w = MDEST(m);
        found = 0;
        if(NEW_DEFECT(w))
        {
          if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
          {
            nnb++;
          }
          found = 1;
        }
        if(!found)
        {
          node = VMYNODE(w);
          while(CORNERTYPE(node))
          {
            node = (NODE *)NFATHER(node);
            fw = NVECTOR(node);
            if(FINE_GRID_DOF(fw))
            {
              if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
              {
                nnb++;
              }
              found = 1;
              break;
            }
          }
        }
        if(!found)
        {
          UserWrite("error in CMGSolve \n");
        }
      }
      i++;
      start[i] = nnb+start[i-1];
    }
  }

  if(i != n)
  {
    UserWrite("error in CMGPreProcess. \n");
    REP_ERR_RETURN(1);
  }

  /* allocate index and matrix array */
  nl = start[n];
  index = (int *) GetTmpMem(MGHEAP(mg),nl*sizeof(int),np->famg_mark_key);
  if (index == NULL)
  {
    UserWrite("ug - cmg: not enough memory. \n");
    REP_ERR_RETURN(1);
  }
  entry = (double *) GetTmpMem(MGHEAP(mg),nl*sizeof(double),np->famg_mark_key);
  if (entry == NULL)
  {
    UserWrite("ug - cmg: not enough memory. \n");
    REP_ERR_RETURN(1);
  }

  /* second step: save matrix entries */
  i = 0; offset = 0;
  for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
  {
    for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
    {
      if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&bmask) && (FINE_GRID_DOF(vec)))
      {
        if(offset != start[i])
        {
          UserWrite("error in CMGPreProcess. \n");
          REP_ERR_RETURN(1);
        }
        entry[offset] = MVALUE(VSTART(vec),mc);
        index[offset] = i;
        offset++;
        for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
        {
          w = MDEST(m);
          found = 0;
          if(FINE_GRID_DOF(w))
          {
            if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
            {
              entry[offset] = MVALUE(m,mc);
              index[offset] = VINDEX(w);
              offset++;
            }
            found = 1;
          }
          if(!found)
          {
            node = VMYNODE(w);
            while(CORNERTYPE(node))
            {
              node = (NODE *)NFATHER(node);
              fw = NVECTOR(node);
              if(FINE_GRID_DOF(fw))
              {
                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                {
                  entry[offset] = MVALUE(m,mc);
                  index[offset] = VINDEX(fw);
                  offset++;
                }
                found = 1;
                break;
              }
            }
          }
          if(!found)
          {
            node = (NODE *)SONNODE(VMYNODE(w));
            ll = lev; ll++;
            while(node != NULL)
            {
              fw = NVECTOR(node);
              if((FINE_GRID_DOF(fw) && (ll < level))
                 || (NEW_DEFECT(fw) && (ll == level)))
              {
                if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
                {
                  entry[offset] = MVALUE(m,mc);
                  index[offset] = VINDEX(fw);
                  offset++;
                }
                found = 1;
                break;
              }
              node = (NODE *)SONNODE(node);
              ll++;
            }
          }
          if(!found)
          {
            UserWrite("error in CMGSolve \n");
          }
        }
        i++;
      }
    }
  }

  grid =  GRID_ON_LEVEL(mg,level);
  for (vec=FIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
  {
    if( (!VSKIPME(vec,0)) && (VDATATYPE(vec)&bmask) && (NEW_DEFECT(vec)))
    {
      if(offset != start[i])
      {
        UserWrite("error in CMGPreProcess. \n");
        REP_ERR_RETURN(1);
      }
      entry[offset] = MVALUE(VSTART(vec),mc);
      index[offset] = i;
      offset++;
      for (m=MNEXT(VSTART(vec)); m!=NULL; m=MNEXT(m))
      {
        w = MDEST(m);
        found = 0;
        if(NEW_DEFECT(w))
        {
          if ( (!VSKIPME(w,0)) &&  (VDATATYPE(w)&xmask))
          {
            entry[offset] = MVALUE(m,mc);
            index[offset] = VINDEX(w);
            offset++;
          }
          found = 1;
        }
        if(!found)
        {
          node = VMYNODE(w);
          while(CORNERTYPE(node))
          {
            node = (NODE *)NFATHER(node);
            fw = NVECTOR(node);
            if(FINE_GRID_DOF(fw))
            {
              if ( (!VSKIPME(fw,0)) &&  (VDATATYPE(fw)&xmask))
              {
                entry[offset] = MVALUE(m,mc);
                index[offset] = VINDEX(fw);
                offset++;
              }
              found = 1;
              break;
            }
          }
        }
        if(!found)
        {
          UserWrite("error in CMGSolve \n");
        }
      }
      i++;
    }
  }

  cmg_interface.n = n;
  cmg_interface.nl = nl;
  cmg_interface.nv = CMG_NVECTORS;
  cmg_interface.start = start;
  cmg_interface.index = index;
  cmg_interface.entry = entry;
  cmg_interface.extra = extra;
  for(j = 0; j < CMG_NVECTORS; j++)
  {
    cmg_interface.vector[j] = vector[j];
  }

  return(0);
}


static INT CMGSolve (NP_ITER *theNP, INT level,
                     VECDATA_DESC *c, VECDATA_DESC *b, MATDATA_DESC *A,
                     INT *result)
{
  MULTIGRID *mg;
  GRID *grid;
  NODE *snode;
  VERTEX *vertex;
  VECTOR *vec,*w,*svec,*fvec;
  MATRIX *mat;
  SHORT cc,bc,mc,cmask;
  INT i,lev,n;
  DOUBLE *unknown, *rhs, *tv, *tvT, *defect,norm,step,eps,x,y;
  PICTURE *thePic;

  mg = NP_MG(theNP);
  grid = GRID_ON_LEVEL(mg,level);

  if (MD_IS_SCALAR(A) && VD_IS_SCALAR(c) && VD_IS_SCALAR(b))
  {
    cc    = VD_SCALCMP(c);
    mc    = MD_SCALCMP(A);
    bc    = VD_SCALCMP(b);
    cmask  = VD_SCALTYPEMASK(c);
  }
  else
  {
    UserWrite("Not a scalar equation. \n");
    REP_ERR_RETURN(1);
  }

  unknown = cmg_interface.vector[CMG_UNKNOWN];
  rhs = cmg_interface.vector[CMG_RHS];
  defect = cmg_interface.vector[CMG_DEFECT];
  tv = cmg_interface.vector[CMG_TVA];
  tvT = cmg_interface.vector[CMG_TVB];
  n = cmg_interface.n;

  i = 0;
  for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
  {
    for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
    {
      if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&cmask) && (FINE_GRID_DOF(vec)))
      {
        vertex = MYVERTEX(VMYNODE(vec));
        rhs[i] = VVALUE(vec,bc);
        tv[i] = 1.0;
        tvT[i] = tv[i];
        i++;
      }
    }
  }

  grid =  GRID_ON_LEVEL(mg,level);
  for (vec=FIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
  {
    if( (!VSKIPME(vec,0)) &&  (VDATATYPE(vec)&cmask) && (NEW_DEFECT(vec)))
    {
      vertex = MYVERTEX(VMYNODE(vec));
      rhs[i] = VVALUE(vec,bc);
      x = XC(vertex);
      y = YC(vertex);
      tv[i] = 1.0;
      tvT[i] = tv[i];
      i++;
    }
  }
  if(i != n)
  {
    UserWrite("error in CMGPreProcess. \n");
    REP_ERR_RETURN(1);
  }
  UserWriteF("unknowns: %d \n",n);

  /* solve */
  CMGSolveSystem(&cmg_interface,&cmg_parameter);

  i = 0;
  for (lev=FULLREFINELEVEL(mg); lev<level; lev++)
  {
    for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!= NULL; vec=SUCCVC(vec))
    {
      if((VDATATYPE(vec)&cmask) && (FINE_GRID_DOF(vec)))
      {
        if(!VSKIPME(vec,0))
        {
          VVALUE(vec,cc) = unknown[i];
          VVALUE(vec,bc) = defect[i];
          i++;
        }
        else
        {
          VVALUE(vec,cc) = 0.0;
          VVALUE(vec,bc) = 0.0;
        }
        /* make solution consistent on higher levels */
        /* we are not responsible for consistency on lower levels */
        snode = (NODE *)SONNODE(VMYNODE(vec));
        while (snode != NULL)
        {
          svec = NVECTOR(snode);
          if (VDATATYPE(svec)&cmask)
          {
            VVALUE(svec,cc) = VVALUE(vec,cc);
            VVALUE(svec,bc) = VVALUE(vec,bc);
          }
          snode = (NODE *)SONNODE(snode);
        }
      }
    }
  }


  grid =  GRID_ON_LEVEL(mg,level);
  for (vec=FIRSTVECTOR(grid); vec!= NULL; vec=SUCCVC(vec))
  {
    if((VDATATYPE(vec)&cmask) && (NEW_DEFECT(vec)))
    {
      if(!VSKIPME(vec,0))
      {
        VVALUE(vec,cc) = unknown[i];
        VVALUE(vec,bc) = defect[i];
        i++;
      }
      else
      {
        VVALUE(vec,cc) = 0.0;
        VVALUE(vec,bc) = 0.0;
      }
    }
  }
  if(i != n)
  {
    UserWrite("error in CMGPreProcess. \n");
    REP_ERR_RETURN(1);
  }



  return(0);
}


static INT CMGPostProcess (NP_ITER *theNP, INT level,
                           VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                           INT *result)
{
  MULTIGRID *mg;
  NP_CMG *np;

  np = (NP_CMG *) theNP;

  mg = NP_MG(theNP);
  ReleaseTmpMem(MGHEAP(mg),np->famg_mark_key);   /* mark in PreProcess */

  return (0);
}

static INT CMGConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = CMGInit;
  theNP->Display = CMGDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = CMGPreProcess;
  np->Iter = CMGSolve;
  np->PostProcess = CMGPostProcess;

  return(0);
}


INT InitCMG ()
{
  if (CreateClass(ITER_CLASS_NAME ".cmg",sizeof(NP_CMG),CMGConstruct))
    REP_ERR_RETURN (__LINE__);

  return(0);

}
