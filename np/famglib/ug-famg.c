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

static struct FAMG_Interface famg_interface;
static struct FAMG_Parameter famg_parameter;

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
} NP_FAMG;

/****************************************************************************/
/*                                                                          */
/* functions                                                                                            */
/*                                                                          */
/****************************************************************************/


static INT FAMGInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_FAMG *np;
  char *str;

  np = (NP_FAMG *) theNP;

  if (ReadArgvINT("h",&(famg_parameter.heap),argc,argv))
    famg_parameter.heap = 1e+7;
  if (ReadArgvINT("n1",&(famg_parameter.n1),argc,argv))
    famg_parameter.n1 = 1;
  if (ReadArgvINT("n2",&(famg_parameter.n2),argc,argv))
    famg_parameter.n2 = 1;
  if (ReadArgvINT("g",&(famg_parameter.gamma),argc,argv))
    famg_parameter.gamma = 1;
  if (ReadArgvINT("cgn",&(famg_parameter.cgnodes),argc,argv))
    famg_parameter.cgnodes = 1;


  famg_parameter.ilut = 1e+10;
  GetStringValueDouble(":famg:ilut",&(famg_parameter.ilut));

  famg_parameter.cgilut = 0.0;
  GetStringValueDouble(":famg:cgilut",&(famg_parameter.cgilut));

  famg_parameter.conloops = 0;
  GetStringValueInt(":famg:conloops",&(famg_parameter.conloops));

  famg_parameter.mincoarse = 0.8;
  GetStringValueDouble(":famg:mincoarse",&(famg_parameter.mincoarse));

  famg_parameter.type = 0;
  GetStringValueInt(":famg:type",&(famg_parameter.type));

  famg_parameter.stv = 0;
  GetStringValueInt(":famg:stv",&(famg_parameter.stv));

  famg_parameter.tol = 0.95;
  GetStringValueDouble(":famg:tol",&(famg_parameter.tol));

  famg_parameter.sigma = 0.45;
  GetStringValueDouble(":famg:sigma",&(famg_parameter.sigma));

  famg_parameter.omegar = 1.0;
  GetStringValueDouble(":famg:omegar",&(famg_parameter.omegar));

  famg_parameter.omegal = 1.0;
  GetStringValueDouble(":famg:omegal",&(famg_parameter.omegal));

  famg_parameter.error1 = 1e-6;
  GetStringValueDouble(":famg:error1",&(famg_parameter.error1));

  famg_parameter.error2 = 1.0;
  GetStringValueDouble(":famg:error2",&(famg_parameter.error2));

  famg_parameter.maxit = 100;
  GetStringValueInt(":famg:maxit",&(famg_parameter.maxit));

  famg_parameter.alimit = 1e-14;
  GetStringValueDouble(":famg:alimit",&(famg_parameter.alimit));

  famg_parameter.rlimit = 1e-10;
  GetStringValueDouble(":famg:rlimit",&(famg_parameter.rlimit));

  famg_parameter.divlimit= 10.0;
  GetStringValueDouble(":famg:divlimit",&(famg_parameter.divlimit));

  famg_parameter.reduction = 1.0;
  GetStringValueDouble(":famg:reduction",&(famg_parameter.reduction));

  strcpy(famg_parameter.solver,"linit");
  str = GetStringVar(":famg:solver");
  if(str != NULL) strcpy(famg_parameter.solver,str);

  strcpy(famg_parameter.presmoother,"fgs");
  str = GetStringVar(":famg:presmoother");
  if(str != NULL) strcpy(famg_parameter.presmoother,str);

  strcpy(famg_parameter.postsmoother,"bgs");
  str = GetStringVar(":famg:postsmoother");
  if(str != NULL) strcpy(famg_parameter.postsmoother,str);

  strcpy(famg_parameter.cgsmoother,"ilut");
  str = GetStringVar(":famg:cgsmoother");
  if(str != NULL) strcpy(famg_parameter.cgsmoother,str);


  np->heap = famg_parameter.heap;
  np->n1 = famg_parameter.n1;
  np->n2 = famg_parameter.n2;
  np->gamma = famg_parameter.gamma;
  np->cgnodes = famg_parameter.cgnodes;
  np->maxit = famg_parameter.maxit;
  np->alimit = famg_parameter.alimit;
  np->rlimit = famg_parameter.rlimit;
  np->divlimit = famg_parameter.divlimit;
  np->reduction = famg_parameter.reduction;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT FAMGDisplay (NP_BASE *theNP)
{
  NP_FAMG *np;

  np = (NP_FAMG *) theNP;

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

static INT FAMGPreProcess  (NP_ITER *theNP, INT level,
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
  double *entry, *vector[FAMG_NVECTORS];
  void **extra;
  NP_FAMG *np;

  np = (NP_FAMG *) theNP;

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
    UserWrite("FAMGCreateSystem: not enough memory. \n");
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
    UserWrite("ug - famg: not enough memory. \n");
    REP_ERR_RETURN(1);
  }

  /* allocate vectors */
  for(j = 0; j < FAMG_NVECTORS; j++)
  {
    vector[j] = (DOUBLE *) GetTmpMem(MGHEAP(mg),n*sizeof(DOUBLE),np->famg_mark_key);
    if (vector[j] == NULL)
    {
      UserWrite("ug - famg: not enough memory. \n");
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
            UserWrite("error in FAMGSolve \n");
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
          UserWrite("error in FAMGSolve \n");
        }
      }
      i++;
      start[i] = nnb+start[i-1];
    }
  }

  if(i != n)
  {
    UserWrite("error in FAMGPreProcess. \n");
    REP_ERR_RETURN(1);
  }

  /* allocate index and matrix array */
  nl = start[n];
  index = (int *) GetTmpMem(MGHEAP(mg),nl*sizeof(int),np->famg_mark_key);
  if (index == NULL)
  {
    UserWrite("ug - famg: not enough memory. \n");
    REP_ERR_RETURN(1);
  }
  entry = (double *) GetTmpMem(MGHEAP(mg),nl*sizeof(double),np->famg_mark_key);
  if (entry == NULL)
  {
    UserWrite("ug - famg: not enough memory. \n");
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
          UserWrite("error in FAMGPreProcess. \n");
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
            UserWrite("error in FAMGSolve \n");
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
        UserWrite("error in FAMGPreProcess. \n");
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
          UserWrite("error in FAMGSolve \n");
        }
      }
      i++;
    }
  }

  famg_interface.n = n;
  famg_interface.nl = nl;
  famg_interface.nv = FAMG_NVECTORS;
  famg_interface.start = start;
  famg_interface.index = index;
  famg_interface.entry = entry;
  famg_interface.extra = extra;
  for(j = 0; j < FAMG_NVECTORS; j++)
  {
    famg_interface.vector[j] = vector[j];
  }

  return(0);
}


static INT FAMGSolve (NP_ITER *theNP, INT level,
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

  unknown = famg_interface.vector[FAMG_UNKNOWN];
  rhs = famg_interface.vector[FAMG_RHS];
  defect = famg_interface.vector[FAMG_DEFECT];
  tv = famg_interface.vector[FAMG_TVA];
  tvT = famg_interface.vector[FAMG_TVB];
  n = famg_interface.n;

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
    UserWrite("error in FAMGPreProcess. \n");
    REP_ERR_RETURN(1);
  }
  UserWriteF("unknowns: %d \n",n);

  /* solve */
  FAMGSolveSystem(&famg_interface,&famg_parameter);

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
    UserWrite("error in FAMGPreProcess. \n");
    REP_ERR_RETURN(1);
  }



  return(0);
}


static INT FAMGPostProcess (NP_ITER *theNP, INT level,
                            VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A,
                            INT *result)
{
  MULTIGRID *mg;
  NP_FAMG *np;

  np = (NP_FAMG *) theNP;

  mg = NP_MG(theNP);
  ReleaseTmpMem(MGHEAP(mg),np->famg_mark_key);   /* mark in PreProcess */

  return (0);
}

static INT FAMGConstruct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init = FAMGInit;
  theNP->Display = FAMGDisplay;
  theNP->Execute = NPIterExecute;

  np = (NP_ITER *) theNP;
  np->PreProcess = FAMGPreProcess;
  np->Iter = FAMGSolve;
  np->PostProcess = FAMGPostProcess;

  return(0);
}


INT InitFAMG ()
{
  if (CreateClass(ITER_CLASS_NAME ".famg",sizeof(NP_FAMG),FAMGConstruct))
    REP_ERR_RETURN (__LINE__);

  if (InitFAMGGraph())
    REP_ERR_RETURN (__LINE__);

  return(0);
}
