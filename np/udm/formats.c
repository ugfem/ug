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

#define MAX_PRINT_SYM           5

/* format for PrintVectorData and PrintMatrixData */
#define VFORMAT                         " %c=%11.4lE"
#define MFORMAT                         " %c%c=%11.4lE"

/* seperators */
#define NAMESEP                         ':'
#define BLANKS                          " \t"
#define LIST_SEP                        " \t,"
#define IN_PARTS                        "in"

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static char default_type_names[MAXVECTORS];

/* printing routine ptrs */
static ConversionProcPtr PrintVectorDataPtr[NVECTYPES];
static ConversionProcPtr PrintMatrixDataPtr[NMATTYPES];

/* print symbol counters and lists */
static INT NPrintVectors=0;
static INT NPrintMatrixs=0;
static VECDATA_DESC *PrintVector[MAX_PRINT_SYM];
static MATDATA_DESC *PrintMatrix[MAX_PRINT_SYM];

/* environment dir and var ids */
static INT theNewFormatDirID;                   /* env type for NewFormat dir           */
static INT theVecVarID;                                 /* env type for VEC_TEMPLATE vars       */
static INT theMatVarID;                                 /* env type for MAT_TEMPLATE vars       */

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*																			*/
/* functions to set, display and change the printing format			                */
/*																			*/
/****************************************************************************/

INT DisplayPrintingFormat ()
{
  INT i;

  if (NPrintVectors==0)
    UserWrite("no vector symbols printed\n");
  else
  {
    UserWrite("printed vector symbols\n");
    for (i=0; i<NPrintVectors; i++)
      UserWriteF("   '%s'\n",ENVITEM_NAME(PrintVector[i]));
  }

  if (NPrintMatrixs==0)
    UserWrite("\nno matrix symbols printed\n");
  else
  {
    UserWrite("\nprinted matrix symbols\n");
    for (i=0; i<NPrintMatrixs; i++)
      UserWriteF("   '%s'\n",ENVITEM_NAME(PrintMatrix[i]));
  }

  return (NUM_OK);
}

/********************************************************/
/* for the following function							*/
/* please keep help comment in commands.c up to date	*/
/********************************************************/

INT SetPrintingFormatCmd (const MULTIGRID *mg, INT argc, char **argv)
{
  VECDATA_DESC *vd;
  MATDATA_DESC *md;
  INT i,j,add,vec;
  char *token;

  for (i=1; i<argc; i++)
    switch (argv[i][0])
    {
    case 'V' :
    case 'M' :
      if (strchr("0+-",argv[i][1])==NULL)
      {
        PrintErrorMessage('E',"setpf","specify 0,+ or - after V or M option");
        REP_ERR_RETURN (1);
      }
      vec = (argv[i][0]=='V');
      if (argv[i][1]=='0')
      {
        if (vec)
          NPrintVectors = 0;
        else
          NPrintMatrixs = 0;
        break;
      }
      add = (argv[i][1]=='+');
      token = strtok(argv[i]+1,BLANKS);                                 /* forget about first token (0,+ or -) */
      while ((token=strtok(NULL,BLANKS))!=NULL)
      {
        if (add)
        {
          if (vec)
          {
            if (NPrintVectors>=MAX_PRINT_SYM)
            {
              PrintErrorMessage('E',"setpf","max number of print vetor symbols exceeded");
              REP_ERR_RETURN (1);
            }
            for (j=0; j<NPrintVectors; j++)
              if (strcmp(token,ENVITEM_NAME(PrintVector[j]))==0)
                break;
            if (j<NPrintVectors) continue;                                                      /* already in list */
            if ((vd=GetVecDataDescByName(mg,token))==NULL)
            {
              PrintErrorMessage('E',"setpf","vector symbol not found");
              REP_ERR_RETURN (1);
            }
            PrintVector[NPrintVectors++] = vd;
          }
          else
          {
            if (NPrintMatrixs>=MAX_PRINT_SYM)
            {
              PrintErrorMessage('E',"setpf","max number of print vetor symbols exceeded");
              REP_ERR_RETURN (1);
            }
            for (j=0; j<NPrintMatrixs; j++)
              if (strcmp(token,ENVITEM_NAME(PrintMatrix[j]))==0)
                break;
            if (j<NPrintMatrixs) continue;                                                      /* already in list */
            if ((md=GetMatDataDescByName(mg,token))==NULL)
            {
              PrintErrorMessage('E',"setpf","matrix symbol not found");
              REP_ERR_RETURN (1);
            }
            PrintMatrix[NPrintMatrixs++] = md;
          }
        }
        else
        {
          if (vec)
          {
            for (j=0; j<NPrintVectors; j++)
              if (strcmp(token,ENVITEM_NAME(PrintVector[j]))==0)
                break;
            if (j==NPrintVectors)
            {
              PrintErrorMessage('W',"setpf","vector symbol not in list");
              continue;
            }
            for (j++; j<NPrintVectors; j++)
              PrintVector[j-1] = PrintVector[j];
            NPrintVectors--;
          }
          else
          {
            for (j=0; j<NPrintMatrixs; j++)
              if (strcmp(token,ENVITEM_NAME(PrintMatrix[j]))==0)
                break;
            if (j==NPrintMatrixs)
            {
              PrintErrorMessage('W',"setpf","matrix symbol not in list");
              continue;
            }
            for (j++; j<NPrintMatrixs; j++)
              PrintMatrix[j-1] = PrintMatrix[j];
            NPrintMatrixs--;
          }
        }
      }
      break;

    default :
      PrintErrorMessageF('E',"setpf","(invalid option '%s')",argv[i]);
      REP_ERR_RETURN (1);
    }

  DisplayPrintingFormat();

  return (NUM_OK);
}

static char *DisplayVecDD (const VECDATA_DESC *vd, INT type, const DOUBLE *data, const char *indent, char *s)
{
  INT i,n,off;

  n = VD_NCMPS_IN_TYPE(vd,type);
  if (n==0) return (s);

  off = VD_OFFSET(vd,type);

  s += sprintf(s,"%s%s:",indent,ENVITEM_NAME(vd));
  for (i=0; i<n; i++)
    s += sprintf(s,VFORMAT,VM_COMP_NAME(vd,off+i),data[VD_CMP_OF_TYPE(vd,type,i)]);

  *(s++) = '\n';

  return (s);
}

/****************************************************************************/
/*D
        PrintTypeVectorData - print selected vector user data for the 'nsr' format

        SYNOPSIS:
        static INT PrintTypeVectorData (INT type, void *data, const char *indent, char *s)

    PARAMETERS:
   .   type - consider only this type
   .   data - user data
   .   indent - is printed at the beginning of lines
   .   s - output string

        DESCRIPTION:
        Print selected vector user data for the 'nsr' format.

        RETURN VALUE:
        INT
   .n   0: ok

        SEE ALSO:
        setformat, showformat
   D*/
/****************************************************************************/

static INT PrintTypeVectorData (INT type, void *data, const char *indent, char *s)
{
  INT i;

  for (i=0; i<NPrintVectors; i++)
    s = DisplayVecDD(PrintVector[i],type,data,indent,s);

  /* remove last \n */
  *s = '\0';

  return(0);
}

static char *DisplayMatDD (const MATDATA_DESC *md, INT type, const DOUBLE *data, const char *indent, char *s)
{
  INT i,j,nr,nc,off;

  nr = MD_ROWS_IN_MTYPE(md,type);
  nc = MD_COLS_IN_MTYPE(md,type);
  if (nr==0) return (s);

  off = MD_MTYPE_OFFSET(md,type);

  for (i=0; i<nr; i++)
  {
    s += sprintf(s,"%s%s:",indent,ENVITEM_NAME(md));
    for (j=0; j<nc; j++)
      s += sprintf(s,MFORMAT,VM_COMP_NAME(md,2*(off+i*nc+j)),VM_COMP_NAME(md,2*(off+i*nc+j)+1),data[MD_IJ_CMP_OF_MTYPE(md,type,i,j)]);
    *(s++) = '\n';
  }

  return (s);
}

/****************************************************************************/
/*D
        PrintTypeMatrixData - print selected matrix user data for the 'nsr' format

        SYNOPSIS:
        static INT PrintTypeMatrixData (INT type, void *data, const char *indent, char *s)

    PARAMETERS:
   .   type - consider this mat type
   .   data - user data
   .   indent - is printed at the beginning of lines
   .   s - output string

        DESCRIPTION:
        Print selected matrix user data for the 'nsr' format.

        RETURN VALUE:
        INT
   .n   0: ok

        SEE ALSO:
        setformat, showformat
   D*/
/****************************************************************************/

static INT PrintTypeMatrixData (INT type, void *data, const char *indent, char *s)
{
  INT i;

  for (i=0; i<NPrintMatrixs; i++)
    s = DisplayMatDD(PrintMatrix[i],type,data,indent,s);

  /* remove last \n */
  *s = '\0';

  return(0);
}

VEC_TEMPLATE *GetVectorTemplate (const FORMAT *theFmt, const char *template)
{
  ENVITEM *item,*dir;
  VEC_TEMPLATE *first;

  if (ChangeEnvDir("/Formats") == NULL) REP_ERR_RETURN (NULL);
  dir = (ENVITEM *)ChangeEnvDir(ENVITEM_NAME(theFmt));
  if (dir == NULL) REP_ERR_RETURN (NULL);
  if (template != NULL)
    for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theVecVarID)
        if (strcmp(ENVITEM_NAME(item),template)==0)
          return ((VEC_TEMPLATE *)item);

  /* no (or nonexisting) template specified: just take first vector template */
  for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theVecVarID)
    {
      first = (VEC_TEMPLATE*) item;

      /* check if there are further vector templates */
      for (item = NEXT_ENVITEM(item); item != NULL; item = NEXT_ENVITEM(item))
        if (ENVITEM_TYPE(item) == theVecVarID)
        {
          PrintErrorMessage('W',"GetVectorTemplate","taking first of several vector templates");
          break;
        }
      return (first);
    }

  REP_ERR_RETURN (NULL);
}

/****************************************************************************/
/*D
        CreateVecDescOfTemplate - create a VECDATA_DESC according to a given template

        SYNOPSIS:
        VECDATA_DESC *CreateVecDescOfTemplate (MULTIGRID *theMG,
                                                                           const char *name, const char *template)

    PARAMETERS:
   .   theMG		- multigrid
   .   name		- name of the VECDATA_DESC
   .   template	- template name (if NULL a template called "name" is taken)

        DESCRIPTION:
        Create a VECDATA_DESC according to a given template created with the format
        of the multigrid.

        RETURN VALUE:
        VECDATA_DESC *
   .n   ptr to VECDATA_DESC if ok
   .n		NULL if an error occured
   D*/
/****************************************************************************/

VECDATA_DESC *CreateVecDescOfTemplate (MULTIGRID *theMG,
                                       const char *name, const char *template)
{
  VECDATA_DESC *vd,*svd;
  VEC_TEMPLATE *vt;
  SUBVEC *subv;
  SHORT *Comp,SubComp[MAX_VEC_COMP],*offset;
  char SubName[MAX_VEC_COMP];
  INT i,j,k,nc,cmp,type;
  char buffer[NAMESIZE];

  if (template != NULL)
    vt = GetVectorTemplate(MGFORMAT(theMG),template);
  else
    vt = GetVectorTemplate(MGFORMAT(theMG),name);
  if (vt == NULL) {
    PrintErrorMessage('E',"CreateVecDescOfTemplate",
                      "no vector template");
    REP_ERR_RETURN(NULL);
  }
  vd = CreateVecDesc(theMG,name,VT_COMPNAMES(vt),VT_COMPS(vt));
  if (vd == NULL) {
    PrintErrorMessage('E',"CreateVecDescOfTemplate",
                      "cannot create vector descriptor");
    REP_ERR_RETURN(NULL);
  }
  VM_LOCKED(vd) = 1;

  /* now create sub vec descs */
  offset = VD_OFFSETPTR(vd);
  Comp = VM_COMPPTR(vd);
  for (i=0; i<VT_NSUB(vt); i++) {
    subv = VT_SUB(vt,i);
    strcpy(buffer,SUBV_NAME(subv));
    strcat(buffer,name);

    /* compute sub components */
    k = 0;
    for (type=0; type<NVECTYPES; type++)
    {
      nc = SUBV_NCOMP(subv,type);
      for (j=0; j<nc; j++)
      {
        cmp = offset[type]+SUBV_COMP(subv,type,j);
        SubComp[k] = Comp[cmp];
        SubName[k] = VT_COMPNAME(vt,cmp);
        k++;
      }
    }
    svd = CreateSubVecDesc(theMG,buffer,SUBV_NCOMPS(subv),SubComp,SubName);
    if (svd == NULL) {
      PrintErrorMessage('E',"CreateVecDescOfTemplate",
                        "cannot create subvector descriptor");
      REP_ERR_RETURN(NULL);
    }
    VM_LOCKED(svd) = 1;
  }

  return (vd);
}

INT CreateVecDescCmd (MULTIGRID *theMG, INT argc, char **argv)
{
  char *token,*template,buffer[NAMESIZE];

  if (ReadArgvChar("t",buffer,argc,argv))
    template = NULL;
  else
    template = buffer;
  token = strtok(argv[0],BLANKS);
  token = strtok(NULL,BLANKS);
  while (token!=NULL) {
    if (CreateVecDescOfTemplate(theMG,token,template) == NULL) {
      PrintErrorMessage('E'," CreateVecDescCmd",
                        "cannot create vector descriptor");
      REP_ERR_RETURN(1);
    }
    token = strtok(NULL,BLANKS);
  }

  return (NUM_OK);
}

MAT_TEMPLATE *GetMatrixTemplate (const FORMAT *theFmt, const char *template)
{
  ENVITEM *item,*dir;
  MAT_TEMPLATE *first;

  if (ChangeEnvDir("/Formats") == NULL) REP_ERR_RETURN (NULL);
  dir = (ENVITEM *)ChangeEnvDir(ENVITEM_NAME(theFmt));
  if (dir == NULL) REP_ERR_RETURN (NULL);
  if (template != NULL)
    for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theMatVarID)
        if (strcmp(ENVITEM_NAME(item),template)==0)
          return ((MAT_TEMPLATE *)item);

  /* no (or nonexisting) template specified: just take first matrix template */
  for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theMatVarID)
    {
      first = (MAT_TEMPLATE *)item;

      for (item = NEXT_ENVITEM(item); item != NULL; item = NEXT_ENVITEM(item))
        if (ENVITEM_TYPE(item) == theMatVarID)
        {
          PrintErrorMessage('W',"GetMatrixTemplate","taking first of several matrix templates");
          break;
        }
      return (first);
    }

  REP_ERR_RETURN (NULL);
}

/****************************************************************************/
/*D
        CreateMatDescOfTemplate - create a MATDATA_DESC according to a given template

        SYNOPSIS:
        MATDATA_DESC *CreateMatDescOfTemplate (MULTIGRID *theMG,
                                                                           const char *name, const char *template)

    PARAMETERS:
   .   theMG		- multigrid
   .   name		- name of the VECDATA_DESC
   .   template	- template name (if NULL a template called "name" is taken)

        DESCRIPTION:
        Create a MATDATA_DESC according to a given template created with the format
        of the multigrid.

        RETURN VALUE:
        MATDATA_DESC *
   .n   ptr to MATDATA_DESC if ok
   .n		NULL if an error occured
   D*/
/****************************************************************************/

MATDATA_DESC *CreateMatDescOfTemplate (MULTIGRID *theMG,
                                       const char *name, const char *template)
{
  MATDATA_DESC *md,*smd;
  MAT_TEMPLATE *mt;
  SUBMAT *subm;
  SHORT *Comp,SubComp[MAX_MAT_COMP],*offset;
  INT i,j,k,nc,cmp,type;
  char SubName[2*MAX_MAT_COMP];
  char buffer[NAMESIZE];

  if (template != NULL)
    mt = GetMatrixTemplate(MGFORMAT(theMG),template);
  else
    mt = GetMatrixTemplate(MGFORMAT(theMG),name);
  if (mt == NULL) {
    PrintErrorMessage('E',"CreateMatDescOfTemplate",
                      "no matrix template");
    REP_ERR_RETURN(NULL);
  }
  md = CreateMatDesc(theMG,name,MT_COMPNAMES(mt),
                     MT_RCOMPS(mt),MT_CCOMPS(mt));
  if (md == NULL) {
    PrintErrorMessage('E',"CreateMatDescOfTemplate",
                      "cannot create matrix descriptor");
    REP_ERR_RETURN(NULL);
  }
  VM_LOCKED(md) = 1;

  /* now create sub mat descs */
  offset = MD_OFFSETPTR(md);
  Comp = VM_COMPPTR(md);
  for (i=0; i<MT_NSUB(mt); i++) {
    subm = MT_SUB(mt,i);
    strcpy(buffer,SUBM_NAME(subm));
    strcat(buffer,name);

    /* compute sub components */
    k = 0;
    for (type=0; type<NMATTYPES; type++)
    {
      nc = SUBM_RCOMP(subm,type)*SUBM_CCOMP(subm,type);
      for (j=0; j<nc; j++)
      {
        ASSERT(k<MAX_MAT_COMP);

        cmp = offset[type]+SUBM_COMP(subm,type,j);
        SubComp[k] = Comp[cmp];
        SubName[2*k]   = MT_COMPNAME(mt,2*cmp);
        SubName[2*k+1] = MT_COMPNAME(mt,2*cmp+1);
        k++;
      }
    }
    smd = CreateSubMatDesc(theMG,buffer,SUBM_RCOMPS(subm),SUBM_CCOMPS(subm),SubComp,SubName);
    if (smd == NULL) {
      PrintErrorMessage('E',"CreateMatDescOfTemplate",
                        "cannot create submatrix descriptor");
      REP_ERR_RETURN(NULL);
    }
    VM_LOCKED(smd) = 1;
  }

  return (md);
}

INT CreateMatDescCmd (MULTIGRID *theMG, INT argc, char **argv)
{
  char *token,*template,buffer[NAMESIZE];

  if (ReadArgvChar("t",buffer,argc,argv))
    template = NULL;
  else
    template = buffer;
  token = strtok(argv[0],BLANKS);
  token = strtok(NULL,BLANKS);
  while (token!=NULL) {
    if (CreateMatDescOfTemplate(theMG,token,template) == NULL) {
      PrintErrorMessage('E'," CreateMatDescCmd",
                        "cannot create matrix descriptor");
      REP_ERR_RETURN(1);
    }
    token = strtok(NULL,BLANKS);
  }

  return (NUM_OK);
}

/****************************************************************************/
/*D
        CreateVecTemplate - create a VEC_TEMPLATE

        SYNOPSIS:
        VEC_TEMPLATE *CreateVecTemplate (const char *name)

    PARAMETERS:
   .   theMG		- multigrid
   .   name		- name of the VEC_TEMPLATE

        DESCRIPTION:
        Create a VEC_TEMPLATE in the /newformat directoy of the environment.

        RETURN VALUE:
        VEC_TEMPLATE *
   .n   ptr to VEC_TEMPLATE if ok
   .n		NULL if an error occured
   D*/
/****************************************************************************/

static VEC_TEMPLATE *CreateVecTemplate (const char *name)
{
  VEC_TEMPLATE *vt;
  char *token;
  INT j;

  if (name == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir("/newformat")==NULL)
    REP_ERR_RETURN(NULL);

  vt = (VEC_TEMPLATE *) MakeEnvItem (name,theVecVarID,sizeof(VEC_TEMPLATE));
  if (vt==NULL) REP_ERR_RETURN (NULL);
  VT_NSUB(vt) = 0;
  token = DEFAULT_NAMES;
  for (j=0; j<MAX(MAX_VEC_COMP,strlen(DEFAULT_NAMES)); j++)
    VT_COMPNAME(vt,j) = token[j];

  return (vt);
}

/****************************************************************************/
/*D
        CreateMatTemplate - create a MAT_TEMPLATE

        SYNOPSIS:
        MAT_TEMPLATE *CreateMatTemplate (const char *name)

    PARAMETERS:
   .   theMG		- multigrid
   .   name		- name of the MAT_TEMPLATE

        DESCRIPTION:
        Create a MAT_TEMPLATE in the /newformat directoy of the environment.

        RETURN VALUE:
        MAT_TEMPLATE *
   .n   ptr to MAT_TEMPLATE if ok
   .n		NULL if an error occured
   D*/
/****************************************************************************/

static MAT_TEMPLATE *CreateMatTemplate (const char *name)
{
  MAT_TEMPLATE *mt;
  INT j;

  if (name == NULL) REP_ERR_RETURN (NULL);
  if (ChangeEnvDir("/newformat")==NULL)
    REP_ERR_RETURN(NULL);
  mt = (MAT_TEMPLATE *) MakeEnvItem (name,theMatVarID,sizeof(MAT_TEMPLATE));
  if (mt==NULL) REP_ERR_RETURN (NULL);
  MT_NSUB(mt) = 0;
  for (j=0; j<2*MAX_MAT_COMP; j++)
    MT_COMPNAME(mt,j) = ' ';

  return (mt);
}

/****************************************************************************/
/*D
        VDmatchesVT - check whether VECDATA_DESC and VEC_TEMPLATE match

        SYNOPSIS:
        INT VDmatchesVT (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt)

    PARAMETERS:
   .   vd			- vec data descriptor
   .   vt			- vector template

        DESCRIPTION:
        This function checks whether a VECDATA_DESC and a VEC_TEMPLATE match, i.e.
        the number of components per type coincide.

        RETURN VALUE:
        INT
   .n   YES: vd and vt match
   .n   NO:  vd and vt do not match
   D*/
/****************************************************************************/

INT VDmatchesVT (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt)
{
  INT tp;

  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_NCMPS_IN_TYPE(vd,tp)!=VT_COMP(vt,tp))
      return (NO);

  return (YES);
}

/****************************************************************************/
/*D
        VDsubDescFromVT - create a VECDATA_DESC as a sub descriptor from a vector template

        SYNOPSIS:
        INT VDsubDescFromVT (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt, INT sub, CONST_VECDATA_DESC_PTR *subvd)

    PARAMETERS:
   .   vd			- make a sub desc of this VECDATA_DESC
   .   vt			- template containing sub descriptor
   .   sub			- index of sub descriptor in template
   .   subvd		- handle to created sub descriptor

        DESCRIPTION:
        This function creates a sub descriptor to a given VECDATA_DESC according to the given
        sub descriptor of a template.

        RETURN VALUE:
        INT
   .n   0: ok
   .n      n: if an error occured
   D*/
/****************************************************************************/

INT VDsubDescFromVT (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt, INT sub, CONST_VECDATA_DESC_PTR *subvd)
{
  FORMAT *fmt;
  SUBVEC *subv;
  SHORT SubComp[MAX_VEC_COMP];
  const SHORT *offset,*Comp;
  INT i,k,l,type,nc,nn,cmp;
  char SubName[MAX_VEC_COMP];

  if (!VDmatchesVT(vd,vt))
    REP_ERR_RETURN(1);

  fmt = MGFORMAT(VD_MG(vd));

  subv = VT_SUB(vt,sub);
  offset = VD_OFFSETPTR(vd);
  Comp = VM_COMPPTR(vd);

  /* compute sub components */
  k = 0;
  for (type=0; type<NVECTYPES; type++)
  {
    nc = SUBV_NCOMP(subv,type);
    nn = VD_NCMPS_IN_TYPE(vd,type);
    for (i=0; i<nc; i++)
    {
      l = SUBV_COMP(subv,type,i);
      if (l>=nn)
        REP_ERR_RETURN (1);
      cmp = offset[type]+l;
      SubComp[k] = Comp[cmp];
      SubName[k] = VT_COMPNAME(vt,cmp);
      k++;
    }
  }
  *subvd = CreateSubVecDesc(VD_MG(vd),NULL,SUBV_NCOMPS(subv),SubComp,SubName);
  if (*subvd == NULL)
    REP_ERR_RETURN (1);

  return (0);
}

/****************************************************************************/
/*D
        MDmatchesMT - check whether MATDATA_DESC and MAT_TEMPLATE match

        SYNOPSIS:
        INT MDmatchesMT (const MATDATA_DESC *md, const MAT_TEMPLATE *mt)

    PARAMETERS:
   .   md			- matrix data descriptor
   .   mt			- matrix template

        DESCRIPTION:
        This function checks whether a MATDATA_DESC and a MAT_TEMPLATE match, i.e.
        the number of row/col components per type coincide.

        RETURN VALUE:
        INT
   .n   YES: vd and vt match
   .n   NO:  vd and vt do not match
   D*/
/****************************************************************************/

INT MDmatchesMT (const MATDATA_DESC *md, const MAT_TEMPLATE *mt)
{
  INT tp;

  for (tp=0; tp<NMATTYPES; tp++)
  {
    if (MD_ROWS_IN_MTYPE(md,tp)!=MT_RCOMP(mt,tp))
      return (NO);
    if (MD_COLS_IN_MTYPE(md,tp)!=MT_CCOMP(mt,tp))
      return (NO);
  }

  return (YES);
}

/****************************************************************************/
/*D
        MDmatchesVT - check whether MATDATA_DESC and VEC_TEMPLATE match (tensor product)

        SYNOPSIS:
        INT MDmatchesVT (const MATDATA_DESC *md, const VEC_TEMPLATE *vt)

    PARAMETERS:
   .   md			- matrix data descriptor
   .   vt			- vector template

        DESCRIPTION:
        This function checks whether a MATDATA_DESC and a VEC_TEMPLATE match, i.e.
        the number of components per type coincide in the sense that md is a
        tensor product of vt.

        RETURN VALUE:
        INT
   .n   YES: md and vt match
   .n   NO:  md and vt do not match
   D*/
/****************************************************************************/

INT MDmatchesVT (const MATDATA_DESC *md, const VEC_TEMPLATE *vt)
{
  INT rt,ct,mt;

  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
    {
      mt = MTP(rt,ct);
      if (MD_ROWS_IN_MTYPE(md,mt)!=VT_COMP(vt,rt))
        return (NO);
      if (MD_COLS_IN_MTYPE(md,mt)!=VT_COMP(vt,ct))
        return (NO);
    }

  return (YES);
}

/****************************************************************************/
/*D
        MDsubDescFromVT - create a MATDATA_DESC as a sub descriptor from a vector template

        SYNOPSIS:
        INT MDsubDescFromVT (const MATDATA_DESC *md, const VEC_TEMPLATE *vt, INT sub, CONST_MATDATA_DESC_PTR *submd)

    PARAMETERS:
   .   md			- make a sub desc of this MATDATA_DESC
   .   vt			- template containing sub descriptor
   .   sub			- index of sub descriptor in template
   .   submd		- handle to created sub descriptor

        DESCRIPTION:
        This function creates a sub descriptor to a given MATDATA_DESC according to the given
        subv of a vector template. The matrix is composed as tensor product of the subv.

        RETURN VALUE:
        INT
   .n   0: ok
   .n      n: if an error occured
   D*/
/****************************************************************************/

INT MDsubDescFromVT (const MATDATA_DESC *md, const VEC_TEMPLATE *vt, INT sub, CONST_MATDATA_DESC_PTR *submd)
{
  FORMAT *fmt;
  SUBVEC *subv;
  SHORT SubComp[MAX_MAT_COMP];
  const SHORT *offset,*Comp;
  SHORT RComp[NMATTYPES];
  SHORT CComp[NMATTYPES];
  INT i,j,k,l,rt,ct,mt,nr,nc,NC,nn,cmp;
  char SubName[MAX_MAT_COMP];

  fmt    = MGFORMAT(MD_MG(md));

  subv   = VT_SUB(vt,sub);
  offset = MD_OFFSETPTR(md);
  Comp   = VM_COMPPTR(md);

  /* compute sub components */
  k  = 0;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
    {
      nr = SUBV_NCOMP(subv,rt);
      nc = SUBV_NCOMP(subv,ct);
      mt = MTP(rt,ct);
      if (nr*nc==0)
      {
        RComp[mt] = CComp[mt] = 0;
        continue;
      }
      RComp[mt] = nr;
      CComp[mt] = nc;
      NC = MD_COLS_IN_RT_CT(md,rt,ct);
      nn = MD_ROWS_IN_RT_CT(md,rt,ct) * NC;
      for (i=0; i<nr; i++)
        for (j=0; j<nc; j++)
        {
          l = SUBV_COMP(subv,rt,i) * NC + SUBV_COMP(subv,ct,j);
          if (l>=nn)
            REP_ERR_RETURN (1);
          cmp = offset[mt]+l;
          SubComp[k] = Comp[cmp];
          SubName[2*k]   = VM_COMP_NAME(md,2*cmp);
          SubName[2*k+1] = VM_COMP_NAME(md,2*cmp+1);
          k++;
        }
    }

  *submd = CreateSubMatDesc(MD_MG(md),NULL,RComp,RComp,SubComp,SubName);
  if (*submd == NULL)
    REP_ERR_RETURN (1);

  return (0);
}

/**********************************************************************************/
/*D
        newformat - init a format and allocate templates for vec and mat descriptors

        DESCRIPTION:
        The 'newformat' command enrols a format for multigrid user data.
        It also creates templates for vector and matrix descriptors.

        SYNTAX:
   .vb
   newformat <format name>
        {$T <type specifier>}*
        {$V <dofs per type list>: <vector template name> <total needed>
                [$comp <comp names>
                        [$sub <sub name> <vector comp name list>]*]}+
        {{$M implicit(<vt>): <matrix template name> <total needed>
                [$sub <sub name> <rows>x<cols> <matrix comp name list>]}
 |
         {$M <matrix size list>: <matrix template name> <total needed>
                [$comp <matrix comp names>
                        [$sub <sub name> <rows>x<cols> <matrix comp name list>]]}}+
    [$d <type name1>x<type name2><connection depth>]
    [$I <type name><mat_size>]
    [$NE]
    [$e <size>]
    [$n <size>]
   .ve

        Use T-option(s) for definition of types (may be omitted, s.b.):~
   .     <type~specifier>			- <type name> in <domain part list>: <object list>
   .     <type~name>				- <character>
   .     <domain~part~list>		- <int> {, <int>}*
   .     <object~list>				- <obj> {, <obj>}*
   .     <obj>						- nd | ed | el | si

        NB:: If no T-option is found at all it is assumed that default types are defined:~
   .n    $T n in 0,...: nd $T k in 0,...: ed $T e in 0,...: el $T s in 0,...: si
        to ensure downward compatibilty.


        Use V-options for definition of vector templates:~
   .     <dofs~per~type~list>		- {<type name><int>}+
   .     <total~needed>			- int for total number of this template needed
                                                                        simultaneously

        Use comp-option following a V-option to specify component names:~
   .     <comp~names>				- string of single chracters, one per dof

        Use sub-option(s) following a comp-option to define sub templates for
        vector templates:~
   .     <vector~comp~name~list>	- any combination of characters from <comp~names>

        Use M-option(s) to define matrix templates as tensor product <vt>x<vt>:~
   .     <vt>						- vector template names defined above

        Use sub-option similar to those for vector templates (component names are defined in
        canonic way):~
   .     <matrix~comp~name~list>	- {<matrix comp name> }+
   .     <matrix~comp~name>		- two characters referring to the component
                                                                        names of vt1/2 (which have to be given
                                                                        above!) indicating row and col

        A second syntax for M-option(s) is still supported:~
   .     <matrix~size~list>		- {<matrix size> }+
   .     <matrix~size>				- <dofs per type list>x<dofs per type list>

        To define sub-matrix-templates as above you then have to specify component
        names with a comp-option following the M-option (second format only).

        Further options:~
   .	  d							- specify connection depth other than 0 (inside element only)
                                                                  for <type name1>x<type name2> connections
   .	  I							- (capital i) specify interpolation matrices
   .	  NE						- require node element lists to be generated
   .	  e							- user data in elements in bytes
   .	  n							- user data in nodes in bytes

        EXAMPLE:
   .vb
   newformat myfmt
        $T a in 0: nd,ed $T b in 1: el		# defines 2 abstract types:
 #      a in domain part 0 including nd and ed objects
 #      b in domain part 1 including el objects

        $V a3b2: vt 5 $comp uvwxy			# template vt contains 3 dofs in nd and ed of part 0
 #					   2 dofs in el of part 1
 # storage reservation for 5 VECDATA_DESCs of vt is made
 # component names in nodes and edges of part 0 (type a) are uvw
 # component names in elements of part 1 (type b) are xy

    $M implicit(vt): mt 2;				# implicit generation of a matrix template mt
 # with storage reservation for 2 of those matrices
   .ve

        KEYWORDS:
        storage, format
   D*/
/**********************************************************************************/

static INT ScanVecOption (      INT argc, char **argv,                  /* option list						*/
                                INT *curropt,                                                           /* next option to scan				*/
                                INT po2t[][MAXVOBJECTS],                                    /* part-obj to type table			*/
                                INT MaxType,                                                            /* bound for type id				*/
                                const char TypeNames[],                                         /* names of types					*/
                                INT TypeUsed[],                                                         /* indicate whether type is used	*/
                                INT *nvec,                                                                      /* just an index for templates		*/
                                SHORT VecStorageNeeded[])                                       /* to accumulate storage per type	*/
{
  VEC_TEMPLATE *vt,*vv;
  SUBVEC *subv;
  INT j,type,nsc[NMATTYPES];
  INT opt;
  SHORT offset[NMATOFFSETS];
  char tpltname[NAMESIZE],*names,*token,tp;
  int n;

  opt = *curropt;

  /* find name seperator */
  if ((names=strchr(argv[opt],NAMESEP))==NULL)
  {
    PrintErrorMessageF('E',"newformat","seperate names by a colon ':' from the description (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  *(names++) = '\0';

  /* create a vector template with default name */
  if (sscanf(names,"%d",&n) == 1)
  {
    PrintErrorMessageF('E',"newformat","specifying a number only is not\n"
                       "supported anymore: see man pages (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  if (sscanf(names,"%s",tpltname) != 1)
  {
    PrintErrorMessageF('E',"newformat","no default name specified (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  (*nvec)++;
  vt = CreateVecTemplate(tpltname);
  if (vt == NULL) {
    PrintErrorMessageF('E',"newformat",
                       "could not allocate environment storage (in '$%s')",argv[opt]);
    REP_ERR_RETURN (2);
  }

  /* read types and sizes */
  for (type=0; type<NVECTYPES; type++) VT_COMP(vt,type) = 0;
  token = strtok(argv[opt]+1,BLANKS);
  while (token!=NULL) {
    if (sscanf(token,"%c%d",&tp,&n)!=2) {
      PrintErrorMessageF('E',"newformat",
                         "could not scan type and size (in '$%s')",argv[opt]);
      REP_ERR_RETURN (1);
    }
    for (type=0; type<MaxType; type++)
      if (tp==TypeNames[type])
        break;
    if (type>=MaxType)
    {
      PrintErrorMessageF('E',"newformat","no valid type name '%c' (in '$%s')",tp,argv[opt]);
      REP_ERR_RETURN (1);
    }
    TypeUsed[type] = TRUE;
    if (VT_COMP(vt,type) !=0 ) {
      PrintErrorMessageF('E',"newformat",
                         "double vector type specification (in '$%s')",argv[opt]);
      REP_ERR_RETURN (1);
    }
    VT_COMP(vt,type) = n;
    token = strtok(NULL,BLANKS);
  }

  /* check next arg for compnames */
  if (opt+1 < argc)
    if (strncmp(argv[opt+1],"comp",4)==0) {
      opt++;
      if (sscanf(argv[opt],"comp %s",VT_COMPNAMES(vt))!=1) {
        PrintErrorMessageF('E',"newformat",
                           "no vector comp names specified with comp option (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      ConstructVecOffsets(VT_COMPS(vt),offset);
      if (strlen(VT_COMPNAMES(vt))!=offset[NVECTYPES]) {
        PrintErrorMessageF('E',"newformat",
                           "number of vector comp names != number of comps (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      /* check next args for subv */
      while ((opt+1<argc) && (strncmp(argv[opt+1],"sub",3)==0)) {
        opt++;
        if (VT_NSUB(vt)>=MAX_SUB) {
          PrintErrorMessageF('E',"newformat",
                             "max number of vector subs exceeded (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }
        subv = AllocEnvMemory(sizeof(SUBVEC));
        if (subv == NULL) {
          PrintErrorMessageF('E',"newformat",
                             "could not allocate environment storage (in '$%s')",argv[opt]);
          REP_ERR_RETURN (2);
        }
        memset(subv,0,sizeof(SUBVEC));
        VT_SUB(vt,VT_NSUB(vt)) = subv;
        VT_NSUB(vt)++;

        /* subv name */
        token = strtok(argv[opt]+3,BLANKS);
        if (token==NULL) {
          PrintErrorMessageF('E',"newformat",
                             "specify name of subv (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }
        strcpy(SUBV_NAME(subv),token);

        /* subv comps */
        for (type=0; type<NVECTYPES; type++) nsc[type] = 0;
        while ((token=strtok(NULL,BLANKS))!=NULL) {
          if (strlen(token)!=1) {
            PrintErrorMessageF('E',"newformat",
                               "specify one char per subv comp (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }
          if (strchr(VT_COMPNAMES(vt),*token)==NULL) {
            PrintErrorMessageF('E',"newformat",
                               "wrong subv comp");
            REP_ERR_RETURN (1);
          }
          /* component relative to template */
          n = strchr(VT_COMPNAMES(vt),*token)
              - VT_COMPNAMES(vt);

          /* corresponding type */
          for (type=0; type<NVECTYPES; type++)
            if (n<offset[type+1]) break;
          if (nsc[type]>=MAX_VEC_COMP) {
            PrintErrorMessageF('E',"newformat",
                               "max number of subv comps exceeded (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }
          SUBV_COMP(subv,type,nsc[type]++) = n-offset[type];
        }
        for (type=0; type<NVECTYPES; type++)
          SUBV_NCOMP(subv,type) = nsc[type];
      }
    }

  /* read names of templates */
  if (sscanf(names,"%s %d",tpltname,&n) != 2)
  {
    /* old style: template list (should be avoided) */
    n = 1;
    token = strtok(names,BLANKS);
    token = strtok(NULL,BLANKS);                /* skip first (we already have that as default) */
    while (token!=NULL) {
      n++;
      (*nvec)++;
      vv = CreateVecTemplate(token);
      if (vv == NULL) {
        PrintErrorMessageF('E',"newformat",
                           "could not allocate environment storage (in '$%s')",argv[opt]);
        REP_ERR_RETURN (2);
      }
      /* everything but the name is as in the default template vt */
      for (type=0; type<NVECTYPES; type++)
        VT_COMP(vv,type) = VT_COMP(vt,type);
      for (j=0; j<MAX_VEC_COMP; j++)
        VT_COMPNAME(vv,j) = VT_COMPNAME(vt,j);
      VT_NSUB(vv) = VT_NSUB(vt);
      for (j=0; j<VT_NSUB(vt); j++)
        VT_SUB(vv,j) = VT_SUB(vt,j);
      token = strtok(NULL,BLANKS);
    }
  }
  /* compute storage needed */
  for (type=0; type<NVECTYPES; type++)
    VecStorageNeeded[type] += n * VT_COMP(vt,type);

  *curropt = opt;

  return (0);
}

static INT ScanMatOption (      INT argc, char **argv,                  /* option list						*/
                                INT *curropt,                                                           /* next option to scan				*/
                                INT po2t[][MAXVOBJECTS],                                    /* part-obj to type table			*/
                                INT MaxType,                                                            /* bound for type id				*/
                                const char TypeNames[],                                         /* names of types					*/
                                const INT TypeUsed[],                                           /* indicate whether type is used	*/
                                INT *nmat,                                                                      /* just an index for templates		*/
                                SHORT MatStorageNeeded[])                                       /* to accumulate storage per type	*/
{
  ENVDIR *dir;
  ENVITEM *item;
  VEC_TEMPLATE *vt;
  SUBVEC *subv;
  MAT_TEMPLATE *mt,*mm;
  SUBMAT *subm;
  INT opt,i,j,k,checksub,type,currtype,rtype,ctype,nsc[NMATTYPES];
  SHORT offset[NMATOFFSETS];
  char tpltname[NAMESIZE],subname[NAMESIZE],*names,*token,rt,ct,*p;
  int n,nr,nc,NC;

  opt = *curropt;

  /* find name seperator */
  if ((names=strchr(argv[opt],NAMESEP))==NULL) {
    PrintErrorMessageF('E',"newformat",
                       "seperate names by a colon ':' from the description (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  *(names++) = '\0';

  /* create a matrix template with default name */
  if (sscanf(names,"%d",&n) == 1)
  {
    PrintErrorMessageF('E',"newformat","specifying a number only is not\n"
                       "supported anymore: see man pages (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  if (sscanf(names,"%s",tpltname) != 1)
  {
    PrintErrorMessageF('E',"newformat","no default name specified (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  (*nmat)++;
  mt = CreateMatTemplate(tpltname);
  if (mt == NULL) {
    PrintErrorMessageF('E',"newformat",
                       "could not allocate environment storage (in '$%s')",argv[opt]);
    REP_ERR_RETURN (2);
  }

  /* read types and sizes */
  checksub = FALSE;
  for (type=0; type<NMATTYPES; type++)
    MT_RCOMP(mt,type) = MT_CCOMP(mt,type) = 0;
  token = strtok(argv[opt]+1,BLANKS);
  if (sscanf(token,"implicit(%s)",tpltname)==1)
  {
    /* define matrix template implicitly by (vt x vt) */

    /* strip trailing ')' and get template */
    if ((p=strchr(tpltname,')'))!=NULL)
      *p = '\0';
    if ((dir=ChangeEnvDir("/newformat"))==NULL)
      REP_ERR_RETURN(2);
    for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theVecVarID)
        if (strcmp(ENVITEM_NAME(item),tpltname)==0)
          break;
    if ((vt=(VEC_TEMPLATE *)item)==NULL)
    {
      PrintErrorMessageF('E',"newformat",
                         "vec template in 'implicit(%s)' not found (in '$%s')",tpltname,argv[opt]);
      REP_ERR_RETURN (2);
    }
    for (rtype=0; rtype<NVECTYPES; rtype++)
      for (ctype=0; ctype<NVECTYPES; ctype++)
      {
        nr = VT_COMP(vt,rtype);
        nc = VT_COMP(vt,ctype);
        if (nr*nc<=0) continue;
        type = MTP(rtype,ctype);
        MT_RCOMP(mt,type) = nr;
        MT_CCOMP(mt,type) = nc;
      }

    if (VT_COMPNAMES(vt)[0]!=' ')
    {
      /* define also compnames */
      ConstructVecOffsets(VT_COMPS(vt),offset);

      p = MT_COMPNAMES(mt);
      for (rtype=0; rtype<NVECTYPES; rtype++)
        for (ctype=0; ctype<NVECTYPES; ctype++)
        {
          nr = VT_COMP(vt,rtype);
          nc = VT_COMP(vt,ctype);
          for (j=0; j<nr; j++)
            for (k=0; k<nc; k++)
            {
              *(p++) = VT_COMPNAMES(vt)[offset[rtype]+j];
              *(p++) = VT_COMPNAMES(vt)[offset[ctype]+k];
            }
        }
      checksub = TRUE;
    }

    ConstructMatOffsets(MT_RCOMPS(mt),
                        MT_CCOMPS(mt),offset);
  }
  else
  {
    while (token!=NULL) {
      if (sscanf(token,"%c%dx%c%d",&rt,&nr,&ct,&nc)!=4) {
        PrintErrorMessageF('E',"newformat",
                           "could not scan type and size (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      for (rtype=0; rtype<MaxType; rtype++)
        if (rt==TypeNames[rtype])
          break;
      if (rtype>=MaxType)
      {
        PrintErrorMessageF('E',"newformat","no valid rtype name '%c' (in '$%s')",rt,argv[opt]);
        REP_ERR_RETURN (1);
      }
      if (!TypeUsed[rtype])
      {
        PrintErrorMessageF('W',"newformat","matrix defined in type '%c' without vector? (in '$%s'),",rt,argv[opt]);
      }
      for (ctype=0; ctype<MaxType; ctype++)
        if (ct==TypeNames[ctype])
          break;
      if (ctype>=MaxType)
      {
        PrintErrorMessageF('E',"newformat","no valid ctype name '%c' (in '$%s')",ct,argv[opt]);
        REP_ERR_RETURN (1);
      }
      if (!TypeUsed[ctype])
      {
        PrintErrorMessageF('W',"newformat","matrix defined in type '%c' without vector? (in '$%s'),",ct,argv[opt]);
      }
      type = MTP(rtype,ctype);
      if (MT_RCOMP(mt,type) !=0 ) {
        PrintErrorMessageF('E',"newformat",
                           "double matrix type specification (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      MT_RCOMP(mt,type) = nr;
      MT_CCOMP(mt,type) = nc;
      token = strtok(NULL,BLANKS);
    }

    /* check next arg for compnames */
    if (opt+1 < argc)
      if (strncmp(argv[opt+1],"comp",4) == 0) {
        opt++;
        if (sscanf(argv[opt],"comp %s",MT_COMPNAMES(mt))!=1) {
          PrintErrorMessageF('E',"newformat",
                             "no matrix comp names specified with comp option (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }
        ConstructMatOffsets(MT_RCOMPS(mt),
                            MT_CCOMPS(mt),offset);
        if (strlen(MT_COMPNAMES(mt))!=2*offset[NMATTYPES]) {
          PrintErrorMessageF('E',"newformat",
                             "number of matrix comp names != number of comps (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }
        checksub = TRUE;
      }
  }

  if (checksub)
    /* check next args for subm */
    while ((opt+1<argc) && (strncmp(argv[opt+1],"sub",3)==0)) {
      opt++;
      if (MT_NSUB(mt)>=MAX_SUB) {
        PrintErrorMessageF('E',"newformat",
                           "max number of matrix subs exceeded (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      subm = AllocEnvMemory(sizeof(SUBMAT));
      if (subm == NULL) {
        PrintErrorMessageF('E',"newformat",
                           "could not allocate environment storage (in '$%s')",argv[opt]);
        REP_ERR_RETURN (2);
      }
      memset(subm,0,sizeof(SUBMAT));
      MT_SUB(mt,MT_NSUB(mt)) = subm;
      MT_NSUB(mt)++;

      /* subm name */
      token = strtok(argv[opt]+3,BLANKS);
      if (token==NULL) {
        PrintErrorMessageF('E',"newformat",
                           "specify name of subm (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      strcpy(SUBM_NAME(subm),token);

      /* check next token first for implicit declaration */
      if ((token=strtok(NULL,BLANKS))==NULL)
      {
        PrintErrorMessageF('E',"newformat",
                           "implicit declaration or size expected (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      if (sscanf(token,expandfmt("implicit(%[a-zA_Z0-9_]/%[a-zA_Z0-9_])"),subname,tpltname)==2)
      {
        /* define sub-matrix implicitly by (s/vt x s/vt) */

        /* get vector template */
        if ((dir=ChangeEnvDir("/newformat"))==NULL)
          REP_ERR_RETURN(2);
        for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
          if (ENVITEM_TYPE(item) == theVecVarID)
            if (strcmp(ENVITEM_NAME(item),tpltname)==0)
              break;
        if ((vt=(VEC_TEMPLATE *)item)==NULL)
        {
          PrintErrorMessageF('E',"newformat",
                             "vec template in 'implicit(%s/%s)' not found (in '$%s')",subname,tpltname,argv[opt]);
          REP_ERR_RETURN (2);
        }

        /* get sub vector */
        for (i=0; i<VT_NSUB(vt); i++)
          if (strcmp(SUBV_NAME(VT_SUB(vt,i)),subname)==0)
            break;
        if (i>=VT_NSUB(vt))
        {
          PrintErrorMessageF('E',"newformat",
                             "sub vector '%s' of template '%s' not found (in '$%s')",subname,tpltname,argv[opt]);
          REP_ERR_RETURN (2);
        }
        subv = VT_SUB(vt,i);

        for (rtype=0; rtype<NVECTYPES; rtype++)
          for (ctype=0; ctype<NVECTYPES; ctype++)
          {
            nr = SUBV_NCOMP(subv,rtype);
            nc = SUBV_NCOMP(subv,ctype);
            if (nr*nc<=0) continue;
            type = MTP(rtype,ctype);
            SUBM_RCOMP(subm,type) = nr;
            SUBM_CCOMP(subm,type) = nc;

            NC = MT_CCOMP(mt,type);
            k  = 0;
            for (i=0; i<nr; i++)
              for (j=0; j<nc; j++)
                SUBM_COMP(subm,type,k++) = SUBV_COMP(subv,rtype,i)*NC + SUBV_COMP(subv,ctype,j);
          }

        continue;                               /* while ((opt+1<argc) && (strncmp(argv[opt+1],"sub",3)==0)) */
      }

      /* subm comps */
      for (type=0; type<NMATTYPES; type++) nsc[type] = 0;
      do
      {
        /* scan size */
        if (sscanf(token,"%dx%d",&nr,&nc)!=2) {
          PrintErrorMessageF('E',"newformat",
                             "specify size of subm (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }
        currtype = NOVTYPE;
        while ((token=strtok(NULL,BLANKS))!=NULL) {
          if (strlen(token)!=2) {
            PrintErrorMessageF('E',"newformat",
                               "specify two chars per subm comp (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }
          for (p=MT_COMPNAMES(mt); *p!='\0'; p+=2)
            if ((p[0]==token[0])&&(p[1]==token[1]))
              break;
          if (*p=='\0') {
            PrintErrorMessageF('E',"newformat",
                               "wrong subm comp (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }
          /* comp relative to template */
          n = (p - MT_COMPNAMES(mt))/2;

          /* corresponding type */
          for (type=0; type<NMATTYPES; type++)
            if (n<offset[type+1]) break;

          if (nsc[type]>=MAX_MAT_COMP) {
            PrintErrorMessageF('E',"newformat",
                               "max number of subm comps exceeded (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }
          if (currtype==NOVTYPE)
            currtype = type;
          else if (type!=currtype)
          {
            PrintErrorMessageF('E',"newformat",
                               "wrong comp type for subm (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }
          SUBM_COMP(subm,type,nsc[type]++) = n-offset[type];
          if (nsc[type]==nr*nc) break;
        }
        SUBM_RCOMP(subm,type) = nr;
        SUBM_CCOMP(subm,type) = nc;
      }
      while ((token=strtok(NULL,BLANKS))!=NULL);
    }


  /* read names of templates */
  if (sscanf(names,"%s %d",tpltname,&n) != 2)
  {
    /* old style: template list (should be avoided) */
    n = 1;
    token = strtok(names,BLANKS);
    token = strtok(NULL,BLANKS);                /* skip first (we already have that as default) */
    while (token!=NULL) {
      n++;
      (*nmat)++;
      mm = CreateMatTemplate(token);
      if (mm == NULL) {
        PrintErrorMessageF('E',"newformat",
                           "could not allocate environment storage (in '$%s')",argv[opt]);
        REP_ERR_RETURN (2);
      }
      /* everything but the name is as in the default template vt */
      for (type=0; type<NMATTYPES; type++) {
        MT_RCOMP(mm,type) = MT_RCOMP(mt,type);
        MT_CCOMP(mm,type) = MT_CCOMP(mt,type);
      }
      for (j=0; j<2*MAX_MAT_COMP; j++)
        MT_COMPNAME(mm,j) = MT_COMPNAME(mt,j);
      MT_NSUB(mm) = MT_NSUB(mt);
      for (j=0; j<MT_NSUB(mt); j++)
        MT_SUB(mm,j) = MT_SUB(mt,j);
      token = strtok(NULL,BLANKS);
    }
  }
  /* compute storage needed */
  for (type=0; type<NMATTYPES; type++)
    MatStorageNeeded[type] += n*MT_RCOMP(mt,type)*MT_CCOMP(mt,type);

  *curropt = opt;

  return (0);
}

static INT ScanDepthOption (INT argc, char **argv,                      /* option list						*/
                            INT *curropt,                                                               /* next option to scan				*/
                            INT MaxType,                                                                /* bound for type id				*/
                            const char TypeNames[],                                             /* names of types					*/
                            const INT TypeUsed[],                                               /* indicate whether type is used	*/
                            SHORT ConnDepth[])                                                          /* connection depth of matrices		*/
{
  INT opt,rtype,ctype;
  char rt,ct;
  int depth;

  opt = *curropt;

  if (sscanf(argv[opt],"d %cx%c%d",&rt,&ct,&depth)!=3)
  {
    PrintErrorMessageF('E',"newformat","could not read connection depth (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  for (rtype=0; rtype<MaxType; rtype++)
    if (rt==TypeNames[rtype])
      break;
  if (rtype>=MaxType)
  {
    PrintErrorMessageF('E',"newformat","no valid rtype name '%c' (in '$%s')",rt,argv[opt]);
    REP_ERR_RETURN (1);
  }
  if (!TypeUsed[rtype])
  {
    PrintErrorMessageF('W',"newformat","depth defined in type '%c' without vector? (in '$%s'),",rt,argv[opt]);
  }
  for (ctype=0; ctype<MaxType; ctype++)
    if (ct==TypeNames[ctype])
      break;
  if (ctype>=MaxType)
  {
    PrintErrorMessageF('E',"newformat","no valid ctype name '%c' (in '$%s')",ct,argv[opt]);
    REP_ERR_RETURN (1);
  }
  if (!TypeUsed[ctype])
  {
    PrintErrorMessageF('W',"newformat","depth defined in type '%c' without vector? (in '$%s'),",ct,argv[opt]);
  }

  ConnDepth[MTP(rtype,ctype)] = depth;
  *curropt = opt;

  return (0);
}

static INT ScanIMatOption (INT argc, char **argv,                       /* option list						*/
                           INT *curropt,                                                                /* next option to scan				*/
                           INT MaxType,                                                                 /* bound for type id				*/
                           const char TypeNames[],                                              /* names of types					*/
                           const INT TypeUsed[],                                                /* indicate whether type is used	*/
                           SHORT ImatTypes[])                                                           /* connection depth of matrices		*/
{
  INT opt,type;
  char tp,*token;
  int n;

  opt = *curropt;

  /* read types and sizes of Interpolation matrix */
  token = strtok(argv[opt]+1,BLANKS);
  while (token!=NULL)
  {
    if (sscanf(token,"%c%d",&tp,&n)!=2)
    {
      PrintErrorMessageF('E',"newformat","could not scan type and size (in '$%s')",argv[opt]);
      REP_ERR_RETURN (1);
    }
    for (type=0; type<MaxType; type++)
      if (tp==TypeNames[type])
        break;
    if (type>=MaxType)
    {
      PrintErrorMessageF('E',"newformat","no valid type name '%c' (in '$%s')",tp,argv[opt]);
      REP_ERR_RETURN (1);
    }
    ImatTypes[type] = n;
    token = strtok(NULL,BLANKS);
  }
  *curropt = opt;
  return (0);
}

static INT ScanTypeOptions (INT argc, char **argv, INT po2t[][MAXVOBJECTS], INT *MaxTypes, char TypeNames[])
{
  INT i,j,opt,found,max;
  INT nparts,partlist[MAXDOMPARTS],part;
  INT nobjs,objlist[MAXDOMPARTS];
  char *objstr,*partstr,c,*token;

  /* init po2t */
  for (i=0; i<MAXDOMPARTS; i++)
    for (j=0; j<MAXVOBJECTS; j++)
      po2t[i][j] = NOVTYPE;

  /* scan type specifications from option list */
  found = max = 0;
  for (opt=1; opt<argc; opt++)
    if (argv[opt][0]=='T')
    {
      found++;

      /* scan type name */
      if (sscanf(argv[opt],"T %c",&c)!=1)
      {
        PrintErrorMessageF('E',"newformat",
                           "type name not found (in '$%s')",argv[opt]);
        REP_ERR_RETURN(1);
      }
      for (i=0; i<max; i++)
        if (c==TypeNames[i])
        {
          PrintErrorMessageF('E',"newformat",
                             "duplicate type names '%c' (in '$%s')",c,argv[opt]);
          REP_ERR_RETURN(1);
        }
      if ((c<FROM_VTNAME) || (TO_VTNAME<c))
      {
        PrintErrorMessageF('E',"newformat",
                           "type name '%c' out of range [%c-%c] (in '$%s')",c,FROM_VTNAME,TO_VTNAME,argv[opt]);
        REP_ERR_RETURN(1);
      }
      TypeNames[max] = c;

      /* seperate object list */
      if ((objstr=strchr(argv[opt],NAMESEP))==NULL)
      {
        PrintErrorMessageF('E',"newformat",
                           "no type sperator ':' found in T-option (in '$%s')",argv[opt]);
        REP_ERR_RETURN(1);
      }
      *(objstr++) = '\0';

      /* scan part list */
      if ((partstr=strstr(argv[opt],IN_PARTS))==NULL)
      {
        PrintErrorMessageF('E',"newformat",
                           "no '%s' token found in T-option (in '$%s')",IN_PARTS,argv[opt]);
        REP_ERR_RETURN(1);
      }
      token = strtok(partstr+strlen(IN_PARTS),LIST_SEP);
      nparts = 0;
      while (token!=NULL)
      {
        if (sscanf(token,"%d",&part)!=1)
        {
          PrintErrorMessageF('E',"newformat",
                             "could not scan parts in part-list (in '$%s')",argv[opt]);
          REP_ERR_RETURN(1);
        }
        if ((part<0) || (MAXDOMPARTS<=part))
        {
          PrintErrorMessageF('E',"newformat",
                             "part out of range [%d-%d] (in '$%s')",0,MAXDOMPARTS-1,argv[opt]);
          REP_ERR_RETURN(1);
        }
        partlist[nparts++] = part;
        token = strtok(NULL,LIST_SEP);
      }

      /* scan object list */
      token = strtok(objstr,LIST_SEP);
      nobjs = 0;
      while (token!=NULL)
      {
        for (i=0; i<MAXVOBJECTS; i++)
          if (strcmp(token,ObjTypeName[i])==0)
            break;
        if (i>=MAXVOBJECTS)
        {
          PrintErrorMessageF('E',"newformat",
                             "could not scan object '%s' in object-list (in '$%s')",token,argv[opt]);
          REP_ERR_RETURN(1);
        }
        objlist[nobjs++] = i;
        token = strtok(NULL,LIST_SEP);
      }

      /* update po2t table */
      for (i=0; i<nparts; i++)
        for (j=0; j<nobjs; j++)
          po2t[partlist[i]][objlist[j]] = max;
      max++;
    }

  if (!found)
  {
    /* no T-option: set default types in part 0 */
    for (max=0; max<MAXVOBJECTS; max++)
    {
      TypeNames[max] = default_type_names[max];
      po2t[0][max] = max;
    }
  }

  *MaxTypes = max;

  return (0);
}

INT CreateFormatCmd (INT argc, char **argv)
{
  FORMAT *newFormat;
  ENVDIR *dir;
  VectorDescriptor vd[MAXVECTORS];
  MatrixDescriptor md[MAXMATRICES];
  INT opt,i,j,size,type,rtype,ctype,nvec,nmat,nvd,nmd;
  INT edata,ndata,nodeelementlist;
  INT po2t[MAXDOMPARTS][MAXVOBJECTS],MaxTypes,TypeUsed[MAXVECTORS];
  SHORT ConnDepth[NMATTYPES],ImatTypes[NVECTYPES];
  SHORT VecStorageNeeded[NVECTYPES],MatStorageNeeded[NMATTYPES];
  char formatname[NAMESIZE],TypeNames[NVECTYPES];
  int n,depth;

  /* scan name of format */
  if ((sscanf(argv[0],expandfmt(CONCAT3(" newformat %",NAMELENSTR,"[ -~]")),formatname)!=1) || (strlen(formatname)==0)) {
    PrintErrorMessage('E',"newformat","no format name specified");
    REP_ERR_RETURN (1);
  }
  if (GetFormat(formatname) != NULL) {
    PrintErrorMessage('W',"newformat","format already exists");
    return (NUM_OK);
  }

  /* install the /newformat directory */
  if (ChangeEnvDir("/")==NULL) {
    PrintErrorMessage('F',"InitFormats","could not changedir to root");
    REP_ERR_RETURN(__LINE__);
  }
  if (MakeEnvItem("newformat",theNewFormatDirID,sizeof(ENVDIR))==NULL)
  {
    PrintErrorMessage('F',"InitFormats",
                      "could not install '/newformat' dir");
    REP_ERR_RETURN(__LINE__);
  }

  /* init */
  for (type=0; type<NVECTYPES; type++)
    ImatTypes[type] = VecStorageNeeded[type] = TypeUsed[type] = 0;
  for (type=0; type<NMATTYPES; type++)
    MatStorageNeeded[type] = ConnDepth[type] = 0;
  for (type=0; type<NMATTYPES; type++) ConnDepth[type] = 0;
  nvec = nmat = 0;
  edata = ndata = nodeelementlist = 0;

  /* scan type option or set default po2t */
  if (ScanTypeOptions(argc,argv,po2t,&MaxTypes,TypeNames))
    REP_ERR_RETURN(1);

  /* scan other options */
  for (opt=1; opt<argc; opt++)
    switch (argv[opt][0])
    {
    case 'T' :
      /* this case hase been handled before */
      break;

    case 'V' :
      if (ScanVecOption(argc,argv,&opt,po2t,MaxTypes,TypeNames,TypeUsed,&nvec,VecStorageNeeded))
        REP_ERR_RETURN(1);
      break;

    case 'M' :
      if (ScanMatOption(argc,argv,&opt,po2t,MaxTypes,TypeNames,TypeUsed,&nmat,MatStorageNeeded))
        REP_ERR_RETURN(1);
      break;

    case 'd' :
      if (ScanDepthOption(argc,argv,&opt,MaxTypes,TypeNames,TypeUsed,ConnDepth))
        REP_ERR_RETURN(1);
      break;

    case 'I' :
      if (ScanIMatOption(argc,argv,&opt,MaxTypes,TypeNames,TypeUsed,ImatTypes))
        REP_ERR_RETURN(1);
      break;

    case 'e' :
      if (sscanf(argv[opt],"e %d",&n) == 1)
        edata = n;
      break;

    case 'n' :
      if (sscanf(argv[opt],"n %d",&n) == 1)
        ndata = n;
      break;

    case 'N' :
      if (argv[opt][1] == 'E')
        nodeelementlist = TRUE;
      break;

    default :
      PrintErrorMessageF('E',"newformat","(invalid option '%s')",argv[opt]);
      REP_ERR_RETURN (1);
    }

  if ((ndata == TRUE) && (nodeelementlist == TRUE)) {
    PrintErrorMessage('E',"newformat","specify either $n or $NE");
    REP_ERR_RETURN (5);
  }

  /* remove types not needed from po2t */
  for (i=0; i<MAXDOMPARTS; i++)
    for (j=0; j<MAXVOBJECTS; j++)
      if (!TypeUsed[po2t[i][j]])
        po2t[i][j] = NOVTYPE;

  /* now we are ready to create the format */

  /* fill degrees of freedom needed */
  nvd = 0;
  for (type=0; type<NVECTYPES; type++)
    if (VecStorageNeeded[type]>0)
    {
      vd[nvd].tp    = type;
      vd[nvd].size  = VecStorageNeeded[type]*sizeof(DOUBLE);
      vd[nvd].name  = TypeNames[type];
      nvd++;
    }

  if (nodeelementlist || ndata) {
    for (opt=0; opt<nvd; opt++)
      if (vd[opt].tp == NODEVEC)
        break;
    if (opt == nvd) {
      PrintErrorMessage('E',"newformat","node data requires node vector");
      REP_ERR_RETURN (5);
    }
  }

  /* fill connections needed */
  nmd = 0;
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=rtype; ctype<NVECTYPES; ctype++)
    {
      type = MTP(rtype,ctype);
      size = MAX(MatStorageNeeded[MTP(rtype,ctype)],MatStorageNeeded[MTP(ctype,rtype)]);

      if (size<= 0) continue;

      depth = MAX(ConnDepth[MTP(rtype,ctype)],ConnDepth[MTP(ctype,rtype)]);

      md[nmd].from  = rtype;
      md[nmd].to    = ctype;
      md[nmd].size  = size*sizeof(DOUBLE);
      md[nmd].depth = depth;
      nmd++;
    }

  /* create format */
  newFormat = CreateFormat(formatname,0,0,
                           (ConversionProcPtr)NULL,(ConversionProcPtr)NULL,(ConversionProcPtr)NULL,
                           PrintTypeVectorData,PrintTypeMatrixData,
                           nvd,vd,nmd,md,ImatTypes,po2t,nodeelementlist,edata,ndata);
  if (newFormat==NULL)
  {
    PrintErrorMessage('E',"newformat","failed creating the format");
    REP_ERR_RETURN (3);
  }

  /* move templates into the new directory */
  dir = ChangeEnvDir("/newformat");
  if (dir == NULL) {
    PrintErrorMessage('E',"newformat","failed moving template");
    REP_ERR_RETURN (4);
  }
  if (ENVITEM_DOWN((ENVDIR *)newFormat) != NULL) {
    PrintErrorMessage('E',"newformat","failed moving template");
    REP_ERR_RETURN (4);
  }
  ENVITEM_DOWN((ENVDIR *)newFormat) = ENVITEM_DOWN(dir);
  ENVITEM_DOWN(dir) = NULL;
  ENVITEM_LOCKED(dir) = 0;
  ChangeEnvDir("/");
  if (RemoveEnvDir((ENVITEM *)dir))
    PrintErrorMessage('W',"InitFormats","could not remove newformat dir");

  return (NUM_OK);
}

/****************************************************************************/
/*D
        RemoveFormatWithSubs - remove format including sub descriptors (iff)

        SYNOPSIS:
        INT RemoveFormatWithSubs (const char *name)

    PARAMETERS:
   .   name - format name

        DESCRIPTION:
        Remove format including sub descriptors (iff). It is not sufficient to
        call DeleteFormat for formats allocated using the CreateFormatCmd of
        formats.c since sub descriptors are allocated directly from the environment
        heap. Calling this function cleans everything and all memory is released.

        RETURN VALUE:
        INT
   .n   0: ok

        SEE ALSO:
        CreateFormatCmd
   D*/
/****************************************************************************/

static INT RemoveTemplateSubs (FORMAT *fmt)
{
  ENVITEM *item;
  VEC_TEMPLATE *vt;
  MAT_TEMPLATE *mt;
  INT i;

  for (item=ENVITEM_DOWN(fmt); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theVecVarID)
    {
      vt = (VEC_TEMPLATE*) item;

      for (i=0; i<VT_NSUB(vt); i++)
        FreeEnvMemory(VT_SUB(vt,i));
      VT_NSUB(vt) = 0;
    }
    else if (ENVITEM_TYPE(item) == theMatVarID)
    {
      mt = (MAT_TEMPLATE*) item;

      for (i=0; i<MT_NSUB(mt); i++)
        FreeEnvMemory(MT_SUB(mt,i));
      MT_NSUB(mt) = 0;
    }
  return (0);
}

INT RemoveFormatWithSubs (const char *name)
{
  FORMAT *fmt;

  fmt = GetFormat(name);
  if (fmt==NULL)
  {
    PrintErrorMessageF('W',"RemoveFormatWithSubs","format '%s' doesn't exist",name);
    return (GM_OK);
  }
  if (RemoveTemplateSubs(fmt))
    REP_ERR_RETURN (1);
  if (DeleteFormat(name))
    REP_ERR_RETURN (1);

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
  INT tp;

  theNewFormatDirID = GetNewEnvDirID();
  theVecVarID = GetNewEnvVarID();
  theMatVarID = GetNewEnvVarID();

  /* init default type names */
  for (tp=0; tp<MAXVECTORS; tp++)
    switch (tp) {
    case NODEVEC : default_type_names[tp] = 'n'; break;
    case EDGEVEC : default_type_names[tp] = 'k'; break;
    case ELEMVEC : default_type_names[tp] = 'e'; break;
    case SIDEVEC : default_type_names[tp] = 's'; break;
    default :
      PrintErrorMessage('E',"newformat","Huh");
      return (__LINE__);
    }

  return (0);
}
