// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      formats.c                                                     */
/*                                                                          */
/* Purpose:   definition of user data and symbols                           */
/*                                                                          */
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de	                                                */
/*																			*/
/* History:   27.03.95 begin, ug version 3.0								*/
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "devices.h"
#include "enrol.h"
#include "compiler.h"
#include "misc.h"
#include "gm.h"
#include "ugenv.h"
#include "ugm.h"
#include "algebra.h"
#include "helpmsg.h"
#include "general.h"

#include "np.h"

#include "formats.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* limits for SYMBOL handling */
#define MAX_VEC                         30
#define MAX_MAT                         6
#define MAX_SUB                         3
#define SYMNAMESIZE                     16

#define MAX_PRINT_SYM           5

/* format for PrintVectorData and PrintMatrixData */
#define VFORMAT                         " %c=%11.4lE"
#define MFORMAT                         " %c%c=%11.4lE"

/* seperators */
#define NAMESEP                         ':'
#define BLANKS                          " \t"

/* macros for SUBVEC */
#define SUBV_NAME(s)            ((s)->Name)
#define SUBV_NCOMPS(s)          ((s)->Comp)
#define SUBV_NCOMP(s,tp)        ((s)->Comp[tp])
#define SUBV_COMP(s,tp)     ((s)->Comps[tp])
#define SUBV_COMPPTR(s)     ((s)->Comps)

/* macros for SUBMAT */
#define SUBM_NAME(s)            ((s)->Name)
#define SUBM_RCOMPS(s)          ((s)->RComp)
#define SUBM_CCOMPS(s)          ((s)->CComp)
#define SUBM_RCOMP(s,tp)        ((s)->RComp[tp])
#define SUBM_CCOMP(s,tp)        ((s)->CComp[tp])
#define SUBM_COMP(s,tp)     ((s)->Comps[tp])
#define SUBM_COMPPTR(s)     ((s)->Comps)

/* macros for VEC_FORMAT */
#define VF_COMPS(vf)            ((vf)->Comp)
#define VF_COMP(vf,tp)          ((vf)->Comp[tp])
#define VF_COMPNAMES(vf)        ((vf)->CompNames)
#define VF_COMPNAME(vf,i)       ((vf)->CompNames[i])
#define VF_FIRSTCOMP(vf,tp)     ((vf)->FirstComp[tp])
#define VF_SUB(vf,i)            ((vf)->SubVec+(i))
#define VF_NSUB(vf)                     ((vf)->nsub)

/* macros for MAT_FORMAT */
#define MF_RCOMPS(mf)           ((mf)->RComp)
#define MF_RCOMP(mf,tp)         ((mf)->RComp[tp])
#define MF_CCOMPS(mf)           ((mf)->CComp)
#define MF_CCOMP(mf,tp)         ((mf)->CComp[tp])
#define MF_COMPNAMES(mf)        ((mf)->CompNames)
#define MF_COMPNAME(mf,i)       ((mf)->CompNames[i])
#define MF_FIRSTCOMP(mf,tp)     ((mf)->FirstComp[tp])
#define MF_SUB(mf,i)            ((mf)->SubMat+(i))
#define MF_NSUB(mf)                     ((mf)->nsub)

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

typedef struct {

  char Name[SYMNAMESIZE];
  SHORT Comp[NVECTYPES];
  SHORT Comps[NVECTYPES];

} SUBVEC;

typedef struct {

  char Name[SYMNAMESIZE];
  SHORT RComp[NMATTYPES];
  SHORT CComp[NMATTYPES];
  SHORT Comps[NMATTYPES];

} SUBMAT;

typedef struct {

  ENVITEM v;

  SHORT Comp[NVECTYPES];
  char CompNames[MAX_VEC_COMP];

  SHORT nsub;
  SUBVEC SubVec[MAX_SUB];

} VEC_FORMAT;

typedef struct {

  ENVITEM v;

  SHORT RComp[NMATTYPES];
  SHORT CComp[NMATTYPES];
  char CompNames[2*MAX_MAT_COMP];

  SHORT nsub;
  SUBMAT SubMat[MAX_SUB];

} MAT_FORMAT;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* printing routine ptrs */
static ConversionProcPtr PrintVectorDataPtr[NVECTYPES];
static ConversionProcPtr PrintMatrixDataPtr[NMATTYPES];

/* print symbol counters and lists */
static INT NPrintVectorSymbols=0;
static INT NPrintMatrixSymbols=0;

/* printing formats */
static INT VectorPrintingFormat = 0x3;
static INT MatrixPrintingFormat = 0xF;
static INT PrintMatrixLine          = 0x3;

/* user data block id */
static BLOCK_ID nsrMGUDid;

/* environment dir and var ids */
static INT theNewFormatDirID;                   /* env type for NewFormat dir           */
static INT theVecVarID;                                 /* env type for VEC_FORMAT vars         */
static INT theMatVarID;                                 /* env type for MAT_FORMAT vars         */

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


static INT PrintTypeVectorData (INT type, void *data, const char *indent, char *s)
{
  INT i;

  /* TODO: fill in print */

  /* remove last \n */
  *s = '\0';

  return(0);
}

static INT PrintNodeVectorData (void *data, const char *indent, char *s)
{
  return (PrintTypeVectorData(NODEVECTOR,data,indent,s));
}

static INT PrintElemVectorData (void *data, const char *indent, char *s)
{
  return (PrintTypeVectorData(ELEMVECTOR,data,indent,s));
}

static INT PrintEdgeVectorData (void *data, const char *indent, char *s)
{
  return (PrintTypeVectorData(EDGEVECTOR,data,indent,s));
}

#ifdef __THREEDIM__
static INT PrintSideVectorData (void *data, const char *indent, char *s)
{
  return (PrintTypeVectorData(SIDEVECTOR,data,indent,s));
}
#endif

static INT PrintTypeMatrixData (INT type, void *data, const char *indent, char *s)
{
  INT i;

  /* TODO: fill in print */

  /* remove last \n */
  *s = '\0';

  return(0);
}

static INT PrintNodeNodeMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(NODEVECTOR,NODEVECTOR),data,indent,s));
}

static INT PrintNodeElemMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(NODEVECTOR,ELEMVECTOR),data,indent,s));
}

static INT PrintNodeEdgeMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(NODEVECTOR,EDGEVECTOR),data,indent,s));
}

#ifdef __THREEDIM__
static INT PrintNodeSideMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(NODEVECTOR,SIDEVECTOR),data,indent,s));
}
#endif

static INT PrintElemElemMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(ELEMVECTOR,ELEMVECTOR),data,indent,s));
}

static INT PrintElemEdgeMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(ELEMVECTOR,EDGEVECTOR),data,indent,s));
}

#ifdef __THREEDIM__
static INT PrintElemSideMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(ELEMVECTOR,SIDEVECTOR),data,indent,s));
}
#endif

static INT PrintEdgeEdgeMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(EDGEVECTOR,EDGEVECTOR),data,indent,s));
}

#ifdef __THREEDIM__
static INT PrintEdgeSideMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(EDGEVECTOR,SIDEVECTOR),data,indent,s));
}

static INT PrintSideSideMatrixData (void *data, const char *indent, char *s)
{
  return (PrintTypeMatrixData(MTP(SIDEVECTOR,SIDEVECTOR),data,indent,s));
}
#endif


static VEC_FORMAT *GetVectorTemplate (MULTIGRID *theMG, char *template)
{
  ENVITEM *item;

  if (ChangeEnvDir("/Multigrid") == NULL) return (NULL);
  item = (ENVITEM *)ChangeEnvDir(ENVITEM_NAME(theMG));
  if (item == NULL) return (NULL);
  if (template != NULL)
    for (item=ENVITEM_DOWN(item); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theVecVarID)
        if (strcmp(ENVITEM_NAME(item),template)==0)
          return ((VEC_FORMAT *)item);
  for (item=ENVITEM_DOWN(item); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theVecVarID)
      return ((VEC_FORMAT *)item);

  return (NULL);
}

VECDATA_DESC *CreateVecDescOfTemplate (MULTIGRID *theMG,
                                       char *name, char *template)
{
  VECDATA_DESC *vd;
  VEC_FORMAT *vf;
  SUBVEC *subv;
  INT i,type;
  char buffer[NAMESIZE];

  if (template != NULL)
    vf = GetVectorTemplate(theMG,template);
  else
    vf = GetVectorTemplate(theMG,name);
  if (vf == NULL) {
    PrintErrorMessage('E',"CreateVecDescOfTemplate",
                      "no vector template");
    return(NULL);
  }
  vd = CreateVecDesc(theMG,name,VF_COMPNAMES(vf),VF_COMPS(vf));
  if (vd == NULL) {
    PrintErrorMessage('E',"CreateVecDescOfTemplate",
                      "cannot create vector descriptor");
    return(NULL);
  }
  for (i=0; i<VF_NSUB(vf); i++) {
    subv = VF_SUB(vf,i);
    strcpy(buffer,SUBV_NAME(subv));
    strcat(buffer,name);
    if (CreateSubVecDesc(theMG,vd,buffer,
                         SUBV_NCOMPS(subv),SUBV_COMPPTR(subv))
        == NULL) {
      PrintErrorMessage('E',"CreateVecDescOfTemplate",
                        "cannot create subvector descriptor");
      return(NULL);
    }
  }

  return (vd);
}

INT CreateVecDescCmd (MULTIGRID *theMG, INT argc, char **argv)
{
  char *token,template[NAMESIZE];

  if (ReadArgvChar("t",template,argc,argv))
    template == NULL;
  token = strtok(argv[0],BLANKS);
  token = strtok(NULL,BLANKS);
  while (token!=NULL) {
    if (CreateVecDescOfTemplate(theMG,token,template) == NULL) {
      PrintErrorMessage('E'," CreateVecDescCmd",
                        "cannot create vector descriptor");
      return(1);
    }
    token = strtok(NULL,BLANKS);
  }

  return (0);
}

static MAT_FORMAT *GetMatrixTemplate (MULTIGRID *theMG, char *template)
{
  ENVITEM *item;

  if (ChangeEnvDir("/Multigrid") == NULL) return (NULL);
  item = (ENVITEM *)ChangeEnvDir(ENVITEM_NAME(theMG));
  if (item == NULL) return (NULL);
  if (template != NULL)
    for (item=ENVITEM_DOWN(item); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theMatVarID)
        if (strcmp(ENVITEM_NAME(item),template)==0)
          return ((MAT_FORMAT *)item);

  for (item=ENVITEM_DOWN(item); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theMatVarID)
      return ((MAT_FORMAT *)item);

  return (NULL);
}

MATDATA_DESC *CreateMatDescOfTemplate (MULTIGRID *theMG,
                                       char *name, char *template)
{
  MATDATA_DESC *md;
  MAT_FORMAT *mf;
  SUBMAT *subm;
  INT i,type;
  char buffer[NAMESIZE];

  if (template != NULL)
    mf = GetMatrixTemplate(theMG,template);
  else
    mf = GetMatrixTemplate(theMG,name);
  if (mf == NULL) {
    PrintErrorMessage('E'," CreateMatDescOfTemplate",
                      "no matrix template");
    return(NULL);
  }
  md = CreateMatDesc(theMG,name,MF_COMPNAMES(mf),
                     MF_RCOMPS(mf),MF_CCOMPS(mf));
  if (md == NULL) {
    PrintErrorMessage('E',"CreateMatDescOfTemplate",
                      "cannot create matrix descriptor");
    return(NULL);
  }
  for (i=0; i<MF_NSUB(mf); i++) {
    subm = MF_SUB(mf,i);
    strcpy(buffer,SUBM_NAME(subm));
    strcat(buffer,name);
    if (CreateSubMatDesc(theMG,md,buffer,SUBM_RCOMPS(subm),
                         SUBM_CCOMPS(subm),SUBM_COMPPTR(subm))
        == NULL) {
      PrintErrorMessage('E'," CreateMatDescOfTemplate",
                        "cannot create submatrix descriptor");
      return(NULL);
    }
  }

  return (md);
}

INT CreateMatDescCmd (MULTIGRID *theMG, INT argc, char **argv)
{
  char *token,template[NAMESIZE];

  if (ReadArgvChar("t",template,argc,argv))
    template == NULL;
  token = strtok(argv[0],BLANKS);
  token = strtok(NULL,BLANKS);
  while (token!=NULL) {
    if (CreateMatDescOfTemplate(theMG,token,template) == NULL) {
      PrintErrorMessage('E'," CreateMatDescCmd",
                        "cannot create matrix descriptor");
      return(1);
    }
    token = strtok(NULL,BLANKS);
  }

  return (0);
}

static VEC_FORMAT *CreateVecTemplate (char *name, SHORT *comp, INT n)
{
  VEC_FORMAT *vf;
  char buffer[NAMESIZE];
  INT type;

  if (ChangeEnvDir("/newformat")==NULL)
    return(NULL);

  if (name == NULL) sprintf(buffer,"vt%02d",n);
  else strcpy(buffer,name);
  vf = (VEC_FORMAT *) MakeEnvItem (buffer,theVecVarID,sizeof(VEC_FORMAT));
  if (vf==NULL) return (NULL);
  for (type=0; type<NVECTYPES; type++)
    VF_COMP(vf,type) = comp[type];
  VF_NSUB(vf) = 0;

  return (vf);
}

static MAT_FORMAT *CreateMatTemplate (char *name,
                                      SHORT *rcomp, SHORT *ccomp,INT n)
{
  MAT_FORMAT *mf;
  char buffer[NAMESIZE];
  INT type;

  if (ChangeEnvDir("/newformat")==NULL)
    return(NULL);
  if (name == NULL) sprintf(buffer,"mt%02d",n);
  else strcpy(buffer,name);
  mf = (MAT_FORMAT *) MakeEnvItem (buffer,theMatVarID,sizeof(MAT_FORMAT));
  if (mf==NULL) return (NULL);
  for (type=0; type<NMATTYPES; type++) {
    MF_RCOMP(mf,type) = rcomp[type];
    MF_CCOMP(mf,type) = ccomp[type];
  }
  MF_NSUB(mf) = 0;

  return (mf);
}

/****************************************************************************/
/*D
        newformat - init a format and allocate symbols

        DESCRIPTION:
        The 'newformat' command enrols a format for multigrid user data.
        It also creates templates for vector and matrix descriptors.

   .vb
        newformat <format_name> [$V <vec_size>: {<n_vec>|<template>*}]
                              [$comp <comp_names> {$sub <sub_name> <comps>}*]]
                            [$M <mat_size>: {<n_mat>|<template>*}
                              [$d <mtype> <depth>]]
                            [$I <mat_size>] [$N] [$e <size>] [$n <size>]
   .ve

   .   V - vector
   .   M - matrix
   .   I - interpolation matrix
   .   N - node element list
   .   e - element data
   .   n - node data

        EXAMPLE:
   .vb
    newformat scalar  $V n1: 4
                      $M n1xn1: 1;

        createvector sol rhs;
        creatematrix MAT;
   .ve
    SEE ALSO:
    'createvector', 'creatematrix'
   D*/
/****************************************************************************/

INT CreateFormatCmd (INT argc, char **argv)
{
  FORMAT *newFormat;
  VectorDescriptor vd[MAXVECTORS];
  MatrixDescriptor md[MAXMATRICES];
  ENVITEM *item;
  ENVDIR *dir;
  VEC_FORMAT *vf;
  MAT_FORMAT *mf;
  SUBVEC *subv;
  SUBMAT *subm;
  INT i,j,size,type,currtype,rtype,ctype,nvec,nmat,nsc[NMATTYPES],nvd,nmd;
  INT edata,ndata,nodeelementlist;
  SHORT offset[NMATOFFSETS],ConnDepth[NMATTYPES],ImatTypes[NVECTYPES];
  SHORT FirstVecComp[NVECTYPES],FirstMatComp[NMATTYPES];
  SHORT VComp[NVECTYPES],RComps[NMATTYPES],CComps[NMATTYPES];
  char formatname[NAMESIZE],*names,*token,tp,rt,ct,*p;
  char buffer[NAMESIZE];
  int n,nr,nc,depth;

  /* scan name of format */
  if (sscanf(argv[0],expandfmt(CONCAT3("newformat %",
                                       NAMELENSTR,"[ -~]")),formatname)!=1) {
    PrintErrorMessage('E',"newformat","no format name specified");
    return (1);
  }
  for (type=0; type<NVECTYPES; type++)
    ImatTypes[type] = FirstVecComp[type] = 0;
  for (type=0; type<NMATTYPES; type++)
    FirstMatComp[type] = ConnDepth[type] = 0;
  for (type=0; type<NVECTYPES; type++) ImatTypes[type] = 0;
  for (type=0; type<NMATTYPES; type++) ConnDepth[type] = 0;
  nvec = nmat = 0;
  edata = ndata = nodeelementlist = 0;
  /* install the /newformat directory */
  if (ChangeEnvDir("/")==NULL) {
    PrintErrorMessage('F',"InitFormats","could not changedir to root");
    return(__LINE__);
  }
  if (MakeEnvItem("newformat",theNewFormatDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitFormats",
                      "could not install '/newformat' dir");
    return(__LINE__);
  }
  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'V' :
      /* create a vector template */

      /* find name seperator */
      if ((names=strchr(argv[i],NAMESEP))==NULL)
      {
        PrintErrorMessage('E',"newformat","seperate names by a colon ':' from the description");
        return (1);
      }
      *(names++) = '\0';

      /* read types and sizes */
      token = strtok(argv[i]+1,BLANKS);
      for (type=0; type<NVECTYPES; type++) VComp[type] = 0;
      while (token!=NULL) {
        if (sscanf(token,"%c%d",&tp,&n)!=2) {
          PrintErrorMessage('E',"newformat",
                            "could not scan type and size");
          return (1);
        }
        switch (tp) {
        case 'n' : type = NODEVECTOR; break;
        case 'k' : type = EDGEVECTOR; break;
        case 'e' : type = ELEMVECTOR; break;
        case 's' : type = SIDEVECTOR; break;
        default :
          PrintErrorMessage('E',"newformat","specify n,k,e,s for the type (or change config to include type)");
          return (1);
        }
        if (VComp[type] !=0 ) {
          PrintErrorMessage('E',"newformat",
                            "double vector type specification");
          return (1);
        }
        VComp[type] = n;
        token = strtok(NULL,BLANKS);
      }
      /* read names of templates */
      if (sscanf(names,"%d",&n) == 1) {
        vf = CreateVecTemplate(NULL,VComp,nvec++);
        if (vf == NULL) {
          PrintErrorMessage('E',"newformat",
                            "could not allocate environment storage");
          return (2);
        }
      }
      else {
        n = 0;
        token = strtok(names,BLANKS);
        while (token!=NULL) {
          n++;
          if (n >= MAX_VEC) {
            PrintErrorMessage('E',"newformat",
                              "max number of main matrix symbols exceeded");
            return (1);
          }
          vf = CreateVecTemplate(buffer,VComp,nvec++);
          if (vf == NULL) {
            PrintErrorMessage('E',"newformat",
                              "could not allocate environment storage");
            return (2);
          }
          token = strtok(NULL,BLANKS);
        }
      }
      for (type=0; type<NVECTYPES; type++)
        FirstVecComp[type] += n * VF_COMP(vf,type);

      break;
      /* check next arg for compnames */
      if (i+1<argc)
        if (strncmp(argv[i+1],"comp",4)!=0)
        {
          for (j=0; j<MAX_VEC_COMP; j++) VF_COMPNAME(vf,j) = ' ';
          break;
        }
      i++;
      if (sscanf(argv[i],"comp %s",VF_COMPNAMES(vf))!=1)
      {
        PrintErrorMessage('E',"newformat","no vector comp names specified with comp option");
        return (1);
      }
      ConstructVecOffsets(VF_COMPS(vf),offset);
      if (strlen(VF_COMPNAMES(vf))!=offset[NVECTYPES])
      {
        PrintErrorMessage('E',"newformat","number of vector comp names != number of comps");
        return (1);
      }

      /* check next args for subv */
      while ((i+1<argc) && (strncmp(argv[i+1],"sub",3)==0))
      {
        i++;

        if (VF_NSUB(vf)>=MAX_SUB)
        {
          PrintErrorMessage('E',"newformat","max number of vector subs exceeded");
          return (1);
        }
        subv = VF_SUB(vf,VF_NSUB(vf));
        VF_NSUB(vf)++;

        /* subv name */
        token = strtok(argv[i]+3,BLANKS);
        if (token==NULL)
        {
          PrintErrorMessage('E',"newformat","specify name of subv");
          return (1);
        }
        strcpy(SUBV_NAME(subv),token);

        /* subv comps */
        for (type=0; type<NVECTYPES; type++) nsc[type] = 0;
        while ((token=strtok(NULL,BLANKS))!=NULL)
        {
          if (strlen(token)!=1)
          {
            PrintErrorMessage('E',"newformat","specify one char per subv comp");
            return (1);
          }
          if (strchr(VF_COMPNAMES(vf),*token)==NULL)
          {
            PrintErrorMessage('E',"newformat","wrong subv comp");
            return (1);
          }
          n = strchr(VF_COMPNAMES(vf),*token) - VF_COMPNAMES(vf);
          for (type=0; type<NVECTYPES; type++)
            if (n<offset[type+1]) break;

          if (nsc[type]>=MAX_VEC_COMP)
          {
            PrintErrorMessage('E',"newformat","max number of subv comps exceeded");
            return (1);
          }
          nsc[type]++;
          SUBV_COMP(subv,type) = n-offset[type];
        }
        for (type=0; type<NVECTYPES; type++) SUBV_NCOMP(subv,type) = nsc[type];
      }
      break;

    case 'M' :
      /* create a matrix */

      /* find name seperator */
      if ((names=strchr(argv[i],NAMESEP))==NULL)
      {
        PrintErrorMessage('E',"newformat","seperate names by a colon ':' from the description");
        return (1);
      }
      *(names++) = '\0';

      /* read types and sizes */
      token = strtok(argv[i]+1,BLANKS);
      for (type=0; type<NMATTYPES; type++)
        RComps[type] = CComps[type] = 0;
      while (token!=NULL)
      {
        if (sscanf(token,"%c%dx%c%d",&rt,&nr,&ct,&nc)!=4)
        {
          PrintErrorMessage('E',"newformat","could not scan type and size");
          return (1);
        }
        switch (rt)
        {
        case 'n' : rtype = NODEVECTOR; break;
        case 'k' : rtype = EDGEVECTOR; break;
        case 'e' : rtype = ELEMVECTOR; break;
        case 's' : rtype = SIDEVECTOR; break;
        default :
          PrintErrorMessage('E',"newformat","specify n,k,e,s for the row type");
          return (1);
        }
        switch (ct)
        {
        case 'n' : ctype = NODEVECTOR; break;
        case 'k' : ctype = EDGEVECTOR; break;
        case 'e' : ctype = ELEMVECTOR; break;
        case 's' : ctype = SIDEVECTOR; break;
        default :
          PrintErrorMessage('E',"newformat","specify n,k,e,s for the col type");
          return (1);
        }
        type = MTP(rtype,ctype);
        if (RComps[type] !=0 ) {
          PrintErrorMessage('E',"newformat",
                            "double matrix type specification");
          return (1);
        }
        RComps[mf,type] = nr;
        CComps[type] = nc;
        token = strtok(NULL,BLANKS);
      }

      /* read names of templates */
      if (sscanf(names,"%d",&n) == 1) {
        mf = CreateMatTemplate(NULL,RComps,CComps,nmat++);
        if (mf == NULL) {
          PrintErrorMessage('E',"newformat",
                            "could not allocate environment storage");
          return (2);
        }
      }
      else {
        n = 0;
        token = strtok(names,BLANKS);
        while (token!=NULL) {
          n++;
          if (n>=MAX_MAT) {
            PrintErrorMessage('E',"newformat",
                              "max number of main matrix symbols exceeded");
            return (1);
          }
          mf = CreateMatTemplate(buffer,RComps,CComps,nmat++);
          if (mf == NULL) {
            PrintErrorMessage('E',"newformat",
                              "could not allocate environment storage");
            return (2);
          }
          token = strtok(NULL,BLANKS);
        }
      }
      for (type=0; type<NMATTYPES; type++)
        FirstMatComp[type] += n*MF_RCOMP(mf,type)*MF_CCOMP(mf,type);

      break;
      /* check next arg for compnames */
      if (i+1<argc)
        if (strncmp(argv[i+1],"comp",4)!=0) {
          for (j=0; j<2*MAX_MAT_COMP; j++)
            MF_COMPNAME(mf,j) = ' ';
          break;
        }
      i++;
      if (sscanf(argv[i],"comp %s",MF_COMPNAMES(mf))!=1) {
        PrintErrorMessage('E',"newformat",
                          "no matrix comp names specified with comp option");
        return (1);
      }
      ConstructMatOffsets(MF_RCOMPS(mf),MF_CCOMPS(mf),offset);
      if (strlen(MF_COMPNAMES(mf))!=2*offset[NMATTYPES])
      {
        PrintErrorMessage('E',"newformat","number of matrix comp names != number of comps");
        return (1);
      }

      /* check next args for subm */
      while ((i+1<argc) && (strncmp(argv[i+1],"sub",3)==0))
      {
        i++;

        if (MF_NSUB(mf)>=MAX_SUB)
        {
          PrintErrorMessage('E',"newformat","max number of matrix subs exceeded");
          return (1);
        }
        subm = MF_SUB(mf,MF_NSUB(mf));
        MF_NSUB(mf)++;

        /* subm name */
        token = strtok(argv[i]+3,BLANKS);
        if (token==NULL)
        {
          PrintErrorMessage('E',"newformat","specify name of subm");
          return (1);
        }
        strcpy(SUBM_NAME(subm),token);

        /* subm comps */
        for (type=0; type<NMATTYPES; type++) nsc[type] = 0;
        while ((token=strtok(NULL,BLANKS))!=NULL)
        {
          /* scan size */
          if (sscanf(token,"%dx%d",&nr,&nc)!=2)
          {
            PrintErrorMessage('E',"newformat","specify size of subm");
            return (1);
          }

          while ((token=strtok(NULL,BLANKS))!=NULL)
          {
            if (strlen(token)!=2)
            {
              PrintErrorMessage('E',"newformat","specify two chars per subm comp");
              return (1);
            }
            for (p=MF_COMPNAMES(mf); *p!='\0'; p+=2)
              if ((p[0]==token[0]) && (p[1]==token[1]))
                break;
            if (*p=='\0')
            {
              PrintErrorMessage('E',"newformat","wrong subm comp");
              return (1);
            }
            n = (p - MF_COMPNAMES(mf))/2;
            for (type=0; type<NMATTYPES; type++)
              if (n<offset[type+1]) break;

            if (nsc[type]>=MAX_MAT_COMP)
            {
              PrintErrorMessage('E',"newformat","max number of subm comps exceeded");
              return (1);
            }
            if (nsc[type]==0)
              currtype = type;
            else if (type!=currtype)
            {
              PrintErrorMessage('E',"newformat","wrong comp type for subm");
              return (1);
            }
            nsc[type]++;
            SUBM_COMP(subm,type) = n-offset[type];

            if (nsc[type]==nr*nc) break;
          }
          SUBM_RCOMP(subm,type) = nr;
          SUBM_CCOMP(subm,type) = nc;
        }
      }
      break;

    case 'd' :
      if (sscanf(argv[i],"d %cx%c%d",&rt,&ct,&depth)!=3)
      {
        PrintErrorMessage('E',"newformat","could not read connection depth");
        return (1);
      }
      switch (rt)
      {
      case 'n' : rtype = NODEVECTOR; break;
      case 'k' : rtype = EDGEVECTOR; break;
      case 'e' : rtype = ELEMVECTOR; break;
      case 's' : rtype = SIDEVECTOR; break;
      default :
        PrintErrorMessage('E',"newformat","specify n,k,e,s for the row type");
        return (1);
      }
      switch (ct)
      {
      case 'n' : ctype = NODEVECTOR; break;
      case 'k' : ctype = EDGEVECTOR; break;
      case 'e' : ctype = ELEMVECTOR; break;
      case 's' : ctype = SIDEVECTOR; break;
      default :
        PrintErrorMessage('E',"newformat","specify n,k,e,s for the col type");
        return (1);
      }
      ConnDepth[MTP(rtype,ctype)] = depth;
      break;

    case 'I' :
      /* read types and sizes of Interpolation matrix */
      token = strtok(argv[i]+1,BLANKS);
      while (token!=NULL)
      {
        if (sscanf(token,"%c%d",&tp,&n)!=2)
        {
          PrintErrorMessage('E',"newformat","could not scan type and size");
          return (1);
        }
        switch (tp)
        {
        case 'n' : type = NODEVECTOR; break;
        case 'k' : type = EDGEVECTOR; break;
        case 'e' : type = ELEMVECTOR; break;
        case 's' : type = SIDEVECTOR; break;
        default :
          PrintErrorMessage('E',"newformat","specify n,k,e,s for the type (or change config to include type)");
          return (1);
        }
        ImatTypes[type] = n;
        token = strtok(NULL,BLANKS);
      }
      break;

    case 'e' :
      if (sscanf(argv[i],"e %d",&n) == 1)
        edata = n;
      break;

    case 'n' :
      if (sscanf(argv[i],"e %d",&n) == 1)
        ndata = n;
      break;

    case 'N' :
      if (argv[i][1] == 'E')
        nodeelementlist = TRUE;
      break;

    default :
      sprintf(buffer,"(invalid option '%s')",argv[i]);
      PrintErrorMessage('E',"newformat",buffer);
      return (1);
    }

  if ((ndata == TRUE) && (nodeelementlist == TRUE)) {
    PrintErrorMessage('E',"newformat","specify $n or $NE");
    return (5);
  }

  /* now we are ready to create the format */

  /* fill degrees of freedom needed */
  nvd = 0;
  for (type=0; type<NVECTYPES; type++)
    if (FirstVecComp[type]>0)
    {
      vd[nvd].pos   = type;
      vd[nvd].size  = FirstVecComp[type]*sizeof(DOUBLE);
      vd[nvd].print = PrintVectorDataPtr[type];
      nvd++;
    }

  if (nodeelementlist || ndata) {
    for (i=0; i<nvd; i++)
      if (vd[i].pos == NODEVECTOR)
        break;
    if (i == nvd) {
      PrintErrorMessage('E',"newformat","node data requires node vector");
      return (5);
    }
  }

  /* fill connections needed */
  nmd = 0;
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype; ctype<NVECTYPES; ctype++)
    {
      type = MTP(rtype,ctype);
      size = MAX(FirstMatComp[MTP(rtype,ctype)],FirstMatComp[MTP(ctype,rtype)]);

      if (size<= 0) continue;

      depth = MAX(ConnDepth[MTP(rtype,ctype)],ConnDepth[MTP(ctype,rtype)]);

      md[nmd].from  = rtype;
      md[nmd].to    = ctype;
      md[nmd].size  = size*sizeof(DOUBLE);
      md[nmd].depth = depth;
      md[nmd].print = PrintMatrixDataPtr[type];
      nmd++;
    }

  /* create format */
  newFormat = CreateFormat(formatname,0,0,
                           (ConversionProcPtr)NULL,(ConversionProcPtr)NULL,(ConversionProcPtr)NULL,
                           nvd,vd,nmd,md);
  if (newFormat==NULL)
  {
    PrintErrorMessage('E',"newformat","failed creating the format");
    return (3);
  }

#ifdef __INTERPOLATION_MATRIX__
  for (i=0; i<MAXVECTORS; i++)
    for (j=0; j<MAXVECTORS; j++)
      newFormat->IMatrixSizes[MatrixType[i][j]]
        = ImatTypes[i] * ImatTypes[j] * sizeof(DOUBLE);
#endif

  newFormat->nodeelementlist = nodeelementlist;
  newFormat->elementdata = edata;
  newFormat->nodedata = ndata;

  /* move tempaltes into the new directory */
  dir = ChangeEnvDir("/newformat");
  if (dir == NULL)
    return (1);
  for (item=ENVITEM_DOWN(dir); item!=NULL; item=NEXT_ENVITEM(item))
    if (MoveEnvItem (item,dir,(ENVDIR *)newFormat))     {
      PrintErrorMessage('E',"newformat","failed moving template");
      return (4);
    }
  if (RemoveEnvDir((ENVITEM *)dir))
    PrintErrorMessage('W',"InitFormats","could not remove newformat dir");

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  InitFormats	                                                */
/*                                                                          */
/* Purpose:   calls all inits of format definitions                         */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    INT 0: everything ok                                          */
/*            INT 1: fatal error (not enough env. space, file not found...  */
/*                                                                          */
/****************************************************************************/

INT InitFormats ()
{
  /* init printing routine ptrs */
        #ifdef __TWODIM__
  PrintVectorDataPtr[NODEVECTOR] = PrintNodeVectorData;
  PrintVectorDataPtr[ELEMVECTOR] = PrintElemVectorData;
  PrintVectorDataPtr[EDGEVECTOR] = PrintEdgeVectorData;

  PrintMatrixDataPtr[MTP(NODEVECTOR,NODEVECTOR)] = PrintNodeNodeMatrixData;
  PrintMatrixDataPtr[MTP(NODEVECTOR,ELEMVECTOR)] =
    PrintMatrixDataPtr[MTP(ELEMVECTOR,NODEVECTOR)] = PrintNodeElemMatrixData;
  PrintMatrixDataPtr[MTP(NODEVECTOR,EDGEVECTOR)] =
    PrintMatrixDataPtr[MTP(EDGEVECTOR,NODEVECTOR)] = PrintNodeEdgeMatrixData;
  PrintMatrixDataPtr[MTP(ELEMVECTOR,ELEMVECTOR)] = PrintElemElemMatrixData;
  PrintMatrixDataPtr[MTP(ELEMVECTOR,EDGEVECTOR)] =
    PrintMatrixDataPtr[MTP(EDGEVECTOR,ELEMVECTOR)] = PrintElemEdgeMatrixData;
  PrintMatrixDataPtr[MTP(EDGEVECTOR,EDGEVECTOR)] = PrintEdgeEdgeMatrixData;
        #endif
        #ifdef __THREEDIM__
  PrintVectorDataPtr[NODEVECTOR] = PrintNodeVectorData;
  PrintVectorDataPtr[ELEMVECTOR] = PrintElemVectorData;
  PrintVectorDataPtr[EDGEVECTOR] = PrintEdgeVectorData;
  PrintVectorDataPtr[SIDEVECTOR] = PrintSideVectorData;

  PrintMatrixDataPtr[MTP(NODEVECTOR,NODEVECTOR)] = PrintNodeNodeMatrixData;
  PrintMatrixDataPtr[MTP(NODEVECTOR,ELEMVECTOR)] =
    PrintMatrixDataPtr[MTP(ELEMVECTOR,NODEVECTOR)] = PrintNodeElemMatrixData;
  PrintMatrixDataPtr[MTP(NODEVECTOR,EDGEVECTOR)] =
    PrintMatrixDataPtr[MTP(EDGEVECTOR,NODEVECTOR)] = PrintNodeEdgeMatrixData;
  PrintMatrixDataPtr[MTP(NODEVECTOR,SIDEVECTOR)] =
    PrintMatrixDataPtr[MTP(SIDEVECTOR,NODEVECTOR)] = PrintNodeSideMatrixData;
  PrintMatrixDataPtr[MTP(ELEMVECTOR,ELEMVECTOR)] = PrintElemElemMatrixData;
  PrintMatrixDataPtr[MTP(ELEMVECTOR,EDGEVECTOR)] =
    PrintMatrixDataPtr[MTP(EDGEVECTOR,ELEMVECTOR)] = PrintElemEdgeMatrixData;
  PrintMatrixDataPtr[MTP(ELEMVECTOR,SIDEVECTOR)] =
    PrintMatrixDataPtr[MTP(SIDEVECTOR,ELEMVECTOR)] = PrintElemSideMatrixData;
  PrintMatrixDataPtr[MTP(EDGEVECTOR,EDGEVECTOR)] = PrintEdgeEdgeMatrixData;
  PrintMatrixDataPtr[MTP(EDGEVECTOR,SIDEVECTOR)] =
    PrintMatrixDataPtr[MTP(SIDEVECTOR,EDGEVECTOR)] = PrintEdgeSideMatrixData;
  PrintMatrixDataPtr[MTP(SIDEVECTOR,SIDEVECTOR)] = PrintSideSideMatrixData;
        #endif

  theNewFormatDirID = GetNewEnvDirID();
  theVecVarID = GetNewEnvVarID();
  theMatVarID = GetNewEnvVarID();

  return(0);
}
