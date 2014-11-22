// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      formats.c                                                     */
/*                                                                          */
/* Purpose:   definition of user data and symbols                           */
/*                                                                          */
/* Author:    Henrik Rentz-Reichert                                         */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   27.03.95 begin, ug version 3.0                                */
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

#include <config.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "ugdevices.h"
#include "enrol.h"
#include "ugtypes.h"
#include "misc.h"
#include "gm.h"
#include "ugenv.h"
#include "ugstruct.h"
#include "ugm.h"
#include "algebra.h"
#include "helpmsg.hh"
#include "general.h"

#include "np.h"

#include "formats.h"

USING_UG_NAMESPACES

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAX_PRINT_SYM                                   5

/** @name Format for PrintVectorData and PrintMatrixData */
/*@{*/
#define VFORMAT                                                 " %c=%11.4E"
#define MFORMAT                                                 " %c%c=%11.4E"
/*@}*/

/** @name Separators */
/*@{*/
#define NAMESEP                                                 ':'
#define BLANKS                                                  " \t"
#define LIST_SEP                                                " \t,"
#define IN_PARTS                                                "in"
/*@}*/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static char default_type_names[MAXVECTORS];

/** @name Print symbol counters and lists */
/*@{*/
static INT NPrintVectors=0;
static INT NPrintMatrixs=0;
static VECDATA_DESC *PrintVector[MAX_PRINT_SYM];
static MATDATA_DESC *PrintMatrix[MAX_PRINT_SYM];
/*@}*/

/** @name Environment dir and var ids */
/*@{*/
static INT theNewFormatDirID;                   /* env type for NewFormat dir           */
static INT theVecVarID;                                 /* env type for VEC_TEMPLATE vars       */
static INT theMatVarID;                                 /* env type for MAT_TEMPLATE vars       */
/*@}*/

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*                                                                          */
/* functions to set, display and change the printing format                 */
/*                                                                          */
/****************************************************************************/

INT NS_DIM_PREFIX DisplayPrintingFormat ()
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

/****************************************************************************/
/** \brief No data will be printed

        After call of this function no data will be printed.
        Do this when closing a multigrid since all descriptors will go out of scope then.

        \return
        0: ok

        \sa
        setformat, showformat
 */
/****************************************************************************/

INT NS_DIM_PREFIX ResetPrintingFormat (void)
{
  NPrintVectors = NPrintMatrixs = 0;
  return (0);
}

/********************************************************/
/* for the following function                           */
/* please keep help comment in commands.c up to date    */
/********************************************************/

INT NS_DIM_PREFIX SetPrintingFormatCmd (const MULTIGRID *mg, INT argc, char **argv)
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
/** \brief Print selected vector user data for the 'nsr' format

   \param type - consider only this type
   \param data - user data
   \param indent - is printed at the beginning of lines
   \param s - output string

        Print selected vector user data for the 'nsr' format.

        \return
        0: ok

        \sa
        setformat, showformat
 */
/****************************************************************************/

static INT PrintTypeVectorData (INT type, void *data, const char *indent, char *s)
{
  INT i;

  for (i=0; i<NPrintVectors; i++)
    s = DisplayVecDD(PrintVector[i],type,(double*)data,indent,s);

  /* remove last \n */
  *s = '\0';

  return(0);
}

static char *DisplayMatDD (const MATDATA_DESC *md, INT type, const DOUBLE *data, const char *indent, char *s)
{
  INT comp,i,j,nr,nc,off;

  nr = MD_ROWS_IN_MTYPE(md,type);
  nc = MD_COLS_IN_MTYPE(md,type);
  if (nr==0) return (s);

  /* diagonals get the same name */
  off = MD_MTYPE_OFFSET(md,MTP(MTYPE_RT(type),MTYPE_CT(type)));

  for (i=0; i<nr; i++)
  {
    s += sprintf(s,"%s%s:",indent,ENVITEM_NAME(md));
    for (j=0; j<nc; j++)
    {
      comp = MD_IJ_CMP_OF_MTYPE(md,type,i,j);
      if (comp<0)
        s += sprintf(s,MFORMAT,VM_COMP_NAME(md,2*(off+i*nc+j)),
                     VM_COMP_NAME(md,2*(off+i*nc+j)+1),
                     0.0);
      else
        s += sprintf(s,MFORMAT,VM_COMP_NAME(md,2*(off+i*nc+j)),
                     VM_COMP_NAME(md,2*(off+i*nc+j)+1),
                     data[MD_IJ_CMP_OF_MTYPE(md,type,i,j)]);
    }
    *(s++) = '\n';
  }

  return (s);
}

/****************************************************************************/
/** \brief Print selected matrix user data for the 'nsr' format

   \param type - consider this mat type
   \param data - user data
   \param indent - is printed at the beginning of lines
   \param s - output string

        Print selected matrix user data for the 'nsr' format.

        \return
        0: ok

        \sa
        setformat, showformat
 */
/****************************************************************************/

static INT PrintTypeMatrixData (INT type, void *data, const char *indent, char *s)
{
  INT i;

  for (i=0; i<NPrintMatrixs; i++)
    s = DisplayMatDD(PrintMatrix[i],type,(double*)data,indent,s);

  /* remove last \n */
  *s = '\0';

  return(0);
}

VEC_TEMPLATE *NS_DIM_PREFIX GetVectorTemplate (const FORMAT *theFmt, const char *theTmplt)
{
  ENVITEM *item,*dir;
  VEC_TEMPLATE *first;

  if (ChangeEnvDir("/Formats") == NULL) REP_ERR_RETURN_PTR (NULL);
  dir = (ENVITEM *)ChangeEnvDir(ENVITEM_NAME(theFmt));
  if (dir == NULL) REP_ERR_RETURN_PTR (NULL);
  if (theTmplt != NULL)
    for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theVecVarID)
        if (strcmp(ENVITEM_NAME(item),theTmplt)==0)
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
          UserWriteF("Trying to handle the vector '%s', but...\n", theTmplt);
          PrintErrorMessage('E',"GetVectorTemplate","there are several vector templates - specify!");
          REP_ERR_RETURN_PTR (NULL);
        }
      return (first);
    }

  REP_ERR_RETURN_PTR (NULL);
}

/****************************************************************************/
/** \brief Create a VECDATA_DESC according to a given template

   \param theMG		- multigrid
   \param name		- name of the VECDATA_DESC
   \param theTmplt	- template name (if NULL a template called "name" is taken)

        Create a VECDATA_DESC according to a given template created with the format
        of the multigrid.

        \return <ul>
        <li> Pointer to VECDATA_DESC if ok </li>
       <li> NULL if an error occured </li>
       </ul>
 */
/****************************************************************************/

VECDATA_DESC * NS_DIM_PREFIX CreateVecDescOfTemplate (MULTIGRID *theMG,
                                                      const char *name, const char *theTmplt)
{
  VECDATA_DESC *vd,*svd;
  VEC_TEMPLATE *vt;
  SUBVEC *subv;
  SHORT *Comp,SubComp[MAX_VEC_COMP],*offset;
  char SubName[MAX_VEC_COMP];
  INT i,j,k,nc,cmp,type;
  char buffer[NAMESIZE];

  if (theTmplt != NULL)
    vt = GetVectorTemplate(MGFORMAT(theMG),theTmplt);
  else
    vt = GetVectorTemplate(MGFORMAT(theMG),name);
  if (vt == NULL) {
    PrintErrorMessage('E',"CreateVecDescOfTemplate",
                      "no vector template");
    REP_ERR_RETURN_PTR(NULL);
  }
  vd = CreateVecDesc(theMG,name,VT_COMPNAMES(vt),VT_COMPS(vt),VT_NID(vt),VT_IDENT_PTR(vt));
  if (vd == NULL) {
    PrintErrorMessage('E',"CreateVecDescOfTemplate",
                      "cannot create vector descriptor");
    REP_ERR_RETURN_PTR(NULL);
  }
  if (LockVD(theMG,vd)) REP_ERR_RETURN_PTR(NULL);

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
      REP_ERR_RETURN_PTR(NULL);
    }
    if (LockVD(theMG,svd)) REP_ERR_RETURN_PTR(NULL);
  }

  return (vd);
}

INT NS_DIM_PREFIX CreateVecDescCmd (MULTIGRID *theMG, INT argc, char **argv)
{
  char *token,*theTmplt,buffer[NAMESIZE];

  if (ReadArgvChar("t",buffer,argc,argv))
    theTmplt = NULL;
  else
    theTmplt = buffer;
  token = strtok(argv[0],BLANKS);
  token = strtok(NULL,BLANKS);
  while (token!=NULL) {
    if (CreateVecDescOfTemplate(theMG,token,theTmplt) == NULL) {
      PrintErrorMessage('E'," CreateVecDescCmd",
                        "cannot create vector descriptor");
      REP_ERR_RETURN(1);
    }
    token = strtok(NULL,BLANKS);
  }

  return (NUM_OK);
}

MAT_TEMPLATE *NS_DIM_PREFIX GetMatrixTemplate (const FORMAT *theFmt, const char *theTmplt)
{
  ENVITEM *item,*dir;
  MAT_TEMPLATE *first;

  if (ChangeEnvDir("/Formats") == NULL) REP_ERR_RETURN_PTR (NULL);
  dir = (ENVITEM *)ChangeEnvDir(ENVITEM_NAME(theFmt));
  if (dir == NULL) REP_ERR_RETURN_PTR (NULL);
  if (theTmplt != NULL)
    for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theMatVarID)
        if (strcmp(ENVITEM_NAME(item),theTmplt)==0)
          return ((MAT_TEMPLATE *)item);

  /* no (or nonexisting) template specified: just take first matrix template */
  for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theMatVarID)
    {
      first = (MAT_TEMPLATE *)item;
      /* check if there are other mat templates */
      for (item = NEXT_ENVITEM(item); item != NULL; item = NEXT_ENVITEM(item))
        if (ENVITEM_TYPE(item) == theMatVarID)
        {
          PrintErrorMessage('W',"GetMatrixTemplate","taking first of several matrix templates");
          REP_ERR_RETURN_PTR (NULL);
        }
      return (first);
    }

  REP_ERR_RETURN_PTR (NULL);
}

/****************************************************************************/
/** \brief Create a MATDATA_DESC according to a given template

   \param theMG		- multigrid
   \param name		- name of the VECDATA_DESC
   \param theTmplt	- template name (if NULL a template called "name" is taken)

        Create a MATDATA_DESC according to a given template created with the format
        of the multigrid.

        \return <ul>
        <li> Pointer to MATDATA_DESC if ok </li>
        <li> NULL if an error occured </li>
        </ul>
 */
/****************************************************************************/

MATDATA_DESC * NS_DIM_PREFIX CreateMatDescOfTemplate (MULTIGRID *theMG,
                                                      const char *name, const char *theTmplt)
{
  MATDATA_DESC *md,*smd;
  MAT_TEMPLATE *mt;
  SUBMAT *subm;
  SHORT *CmpsInType[NMATTYPES],SubComp[MAX_MAT_COMP];
  INT i,j,k,nc,cmp,type;
  char SubName[2*MAX_MAT_COMP];
  char buffer[NAMESIZE];

  if (theTmplt != NULL)
    mt = GetMatrixTemplate(MGFORMAT(theMG),theTmplt);
  else
    mt = GetMatrixTemplate(MGFORMAT(theMG),name);
  if (mt == NULL) {
    PrintErrorMessage('E',"CreateMatDescOfTemplate",
                      "no matrix template");
    REP_ERR_RETURN_PTR(NULL);
  }
  md = CreateMatDesc(theMG,name,MT_COMPNAMES(mt),
                     MT_RCOMPS(mt),MT_CCOMPS(mt),MT_MCMPPTR(mt));
  if (md == NULL) {
    PrintErrorMessage('E',"CreateMatDescOfTemplate",
                      "cannot create matrix descriptor");
    REP_ERR_RETURN_PTR(NULL);
  }
  if (LockMD(md)) REP_ERR_RETURN_PTR(NULL);

  /* now create sub mat descs */
  for (i=0; i<MT_NSUB(mt); i++) {
    subm = MT_SUB(mt,i);
    strcpy(buffer,SUBM_NAME(subm));
    strcat(buffer,name);

    /* compute actual sub components in Subcomp, index it with CmpsInType */
    k = 0;
    for (type=0; type<NMATTYPES; type++)
    {
      nc = SUBM_RCOMP(subm,type)*SUBM_CCOMP(subm,type);
      CmpsInType[type] = &(SubComp[k]);
      for (j=0; j<nc; j++)
      {
        ASSERT(k<MAX_MAT_COMP);

        cmp = SUBM_MCMP_OF_MTYPE(subm,type,j);
        SubComp[k] = MD_MCMP_OF_MTYPE(md,type,cmp);
        SubName[2*k]   = MT_COMPNAME(mt,2*cmp);
        SubName[2*k+1] = MT_COMPNAME(mt,2*cmp+1);
        k++;
      }
    }
    smd = CreateSubMatDesc(theMG,buffer,SubName,SUBM_RCOMPS(subm),SUBM_CCOMPS(subm),CmpsInType);
    if (smd == NULL) {
      PrintErrorMessage('E',"CreateMatDescOfTemplate",
                        "cannot create submatrix descriptor");
      REP_ERR_RETURN_PTR(NULL);
    }
    if (LockMD(smd)) REP_ERR_RETURN_PTR(NULL);
  }

  return (md);
}

INT NS_DIM_PREFIX CreateMatDescCmd (MULTIGRID *theMG, INT argc, char **argv)
{
  char *token,*theTmplt,buffer[NAMESIZE];

  if (ReadArgvChar("t",buffer,argc,argv))
    theTmplt = NULL;
  else
    theTmplt = buffer;
  token = strtok(argv[0],BLANKS);
  token = strtok(NULL,BLANKS);
  while (token!=NULL) {
    if (CreateMatDescOfTemplate(theMG,token,theTmplt) == NULL) {
      PrintErrorMessage('E'," CreateMatDescCmd",
                        "cannot create matrix descriptor");
      REP_ERR_RETURN(1);
    }
    token = strtok(NULL,BLANKS);
  }

  return (NUM_OK);
}

INT NS_DIM_PREFIX FreeMatDescCmd (MULTIGRID *theMG, INT argc, char **argv)
{
  char *token;
  MATDATA_DESC *theMD;

  token = strtok(argv[0],BLANKS);
  token = strtok(NULL,BLANKS);
  while (token!=NULL) {
    theMD = GetMatDataDescByName (theMG, token);
    if (theMD == NULL) {
      PrintErrorMessage('E',"FreeMatDescCmd",
                        "could not find MD");
      REP_ERR_RETURN(-1);
    }
    UnlockMD(theMD);
    if (FreeMD (theMG, 0, TOPLEVEL(theMG),theMD)) {
      PrintErrorMessage('E',"FreeMatDescCmd",
                        "could not free MD");
      REP_ERR_RETURN(-1);
    }
    token = strtok(NULL,BLANKS);
  }

  return (NUM_OK);
}

/****************************************************************************/
/** \brief Create a VEC_TEMPLATE

   \param theMG		- multigrid
   \param name		- name of the VEC_TEMPLATE

        Create a VEC_TEMPLATE in the /newformat directory of the environment.

        \return <ul>
        <li> Pointer to VEC_TEMPLATE if ok </li>
        <li> NULL if an error occured </li>
        </ul>
 */
/****************************************************************************/

static VEC_TEMPLATE *CreateVecTemplate (const char *name)
{
  VEC_TEMPLATE *vt;
  const char *token;
  INT j;

  if (name == NULL) REP_ERR_RETURN_PTR (NULL);
  if (ChangeEnvDir("/newformat")==NULL)
    REP_ERR_RETURN_PTR(NULL);

  vt = (VEC_TEMPLATE *) MakeEnvItem (name,theVecVarID,sizeof(VEC_TEMPLATE));
  if (vt==NULL) REP_ERR_RETURN_PTR (NULL);
  VT_NSUB(vt) = 0;
  VT_NID(vt) = NO_IDENT;
  token = DEFAULT_NAMES;
  for (j=0; j<MAX(MAX_VEC_COMP,strlen(DEFAULT_NAMES)); j++)
    VT_COMPNAME(vt,j) = token[j];

  return (vt);
}

/****************************************************************************/
/** \brief Create a MAT_TEMPLATE

   \param theMG		- multigrid
   \param name		- name of the MAT_TEMPLATE

   Create a MAT_TEMPLATE in the /newformat directory of the environment.

   \return <ul>
   <li> Pointer to MAT_TEMPLATE if ok </li>
   <li> NULL if an error occured </li>
   </ul>
 */
/****************************************************************************/

static MAT_TEMPLATE *CreateMatTemplate (const char *name)
{
  MAT_TEMPLATE *mt;
  INT j;

  if (name == NULL) REP_ERR_RETURN_PTR (NULL);
  if (ChangeEnvDir("/newformat")==NULL)
    REP_ERR_RETURN_PTR(NULL);
  mt = (MAT_TEMPLATE *) MakeEnvItem (name,theMatVarID,sizeof(MAT_TEMPLATE));
  if (mt==NULL) REP_ERR_RETURN_PTR (NULL);
  MT_NSUB(mt) = 0;
  for (j=0; j<2*MAX_MAT_COMP; j++)
    MT_COMPNAME(mt,j) = ' ';

  return (mt);
}

/****************************************************************************/
/** \brief Check whether VECDATA_DESC and VEC_TEMPLATE match

   \param vd			- vec data descriptor
   \param vt			- vector template

   This function checks whether a VECDATA_DESC and a VEC_TEMPLATE match, i.e.
   the number of components per type coincide.

   \return <ul>
   <li> YES: vd and vt match </li>
   <li> NO:  vd and vt do not match </li>
   </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX VDmatchesVT (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt)
{
  INT tp;

  for (tp=0; tp<NVECTYPES; tp++)
    if (VD_NCMPS_IN_TYPE(vd,tp)!=VT_COMP(vt,tp))
      return (NO);

  return (YES);
}

/****************************************************************************/
/** \brief Create a VECDATA_DESC as a sub descriptor from a vector template

   \param vd			- make a sub desc of this VECDATA_DESC
   \param vt			- template containing sub descriptor
   \param sub			- index of sub descriptor in template
   \param subvd		- handle to created sub descriptor

   This function creates a sub descriptor to a given VECDATA_DESC according to the given
   sub descriptor of a template.

   \return <ul>
   <li>         0: ok
   <li>     n: if an error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX VDsubDescFromVT (const VECDATA_DESC *vd, const VEC_TEMPLATE *vt, INT sub, VECDATA_DESC **subvd)
{
  FORMAT *fmt;
  SUBVEC *subv;
  SHORT SubComp[MAX_VEC_COMP];
  const SHORT *offset,*Comp;
  INT i,k,l,type,nc,nn,cmp;
  char SubName[MAX_VEC_COMP],buffer[NAMESIZE];

  if (!VDmatchesVT(vd,vt))
    REP_ERR_RETURN(1);

  ASSERT(sub<VT_NSUB(vt));

  subv = VT_SUB(vt,sub);

  /* generate name and see if desc already exists */
  strcpy(buffer,SUBV_NAME(subv));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,ENVITEM_NAME(vd));
  *subvd = GetVecDataDescByName(MD_MG(vd),buffer);
  if (*subvd != NULL) {
    if (TransmitLockStatusVD(vd,*subvd))
      REP_ERR_RETURN(1);
    return(0);
  }

  fmt = MGFORMAT(VD_MG(vd));

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
  *subvd = CreateSubVecDesc(VD_MG(vd),buffer,SUBV_NCOMPS(subv),SubComp,SubName);
  if (*subvd == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusVD(vd,*subvd))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/** \brief
        VDsubDescFromVS - create a VECDATA_DESC as a vector sub descriptor

        SYNOPSIS:
        INT VDsubDescFromVS (const VECDATA_DESC *vd, const SUBVEC *subv, VECDATA_DESC **subvd)

    PARAMETERS:
   \param vd			- make a sub desc of this VECDATA_DESC
   \param subv		- sub vector descriptor
   \param subvd		- handle to created sub descriptor

        DESCRIPTION:
        This function creates a sub descriptor to a given VECDATA_DESC according to the given
        vector sub descriptor. If a template is available, 'VDsubDescFromVT' should be
        preferred.

        \return <ul>
        INT
   <li>         0: ok
   <li>     n: if an error occured
 */
/****************************************************************************/

INT NS_DIM_PREFIX VDsubDescFromVS (const VECDATA_DESC *vd, const SUBVEC *subv, VECDATA_DESC **subvd)
{
  FORMAT *fmt;
  SHORT SubComp[MAX_VEC_COMP];
  const SHORT *offset,*Comp;
  INT i,k,l,type,nc,nn,cmp;
  char SubName[MAX_VEC_COMP],buffer[NAMESIZE];

  /* generate name and see if desc already exists */
  strcpy(buffer,SUBV_NAME(subv));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,ENVITEM_NAME(vd));
  *subvd = GetVecDataDescByName(MD_MG(vd),buffer);
  if (*subvd != NULL) {
    if (TransmitLockStatusVD(vd,*subvd))
      REP_ERR_RETURN(1);
    return(0);
  }

  fmt = MGFORMAT(VD_MG(vd));

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
      SubName[k] = VM_COMP_NAME(vd,cmp);
      k++;
    }
  }
  *subvd = CreateSubVecDesc(VD_MG(vd),buffer,SUBV_NCOMPS(subv),SubComp,SubName);
  if (*subvd == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusVD(vd,*subvd))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/** \brief
        CompMDwithMT - check whether MATDATA_DESC and MAT_TEMPLATE match

        SYNOPSIS:
        INT CompMDwithMT (const MATDATA_DESC *md, const MAT_TEMPLATE *mt)

    PARAMETERS:
   \param md			- matrix data descriptor
   \param mt			- matrix template

        DESCRIPTION:
        This function checks whether a MATDATA_DESC and a MAT_TEMPLATE match, i.e.
        the number of row/col components per type coincide.

        \return <ul>
        INT
   <li>         0: vd and vt match
   <li>         1: sizes do not match
   <li>     2: sparse structure does not match
 */
/****************************************************************************/

INT NS_DIM_PREFIX CompMDwithMT (const MATDATA_DESC *md, const MAT_TEMPLATE *mt)
{
  return (CompMatDesc(md,MT_RCOMPS(mt),MT_CCOMPS(mt),mt->CmpsInType));
}

/****************************************************************************/
/** \brief
        MDmatchesVT - check whether MATDATA_DESC and VEC_TEMPLATE match (tensor product)

        SYNOPSIS:
        INT MDmatchesVT (const MATDATA_DESC *md, const VEC_TEMPLATE *vt)

    PARAMETERS:
   \param md			- matrix data descriptor
   \param vt			- vector template

        DESCRIPTION:
        This function checks whether a MATDATA_DESC and a VEC_TEMPLATE match, i.e.
        the number of components per type coincide in the sense that md is a
        tensor product of vt.

        \return <ul>
        INT
   <li>         YES: md and vt match
   <li>         NO:  md and vt do not match
 */
/****************************************************************************/

INT NS_DIM_PREFIX MDmatchesVT (const MATDATA_DESC *md, const VEC_TEMPLATE *vt)
{
  INT rt,ct,mt,nr,nc;

  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
    {
      nr = VT_COMP(vt,rt);
      nc = VT_COMP(vt,ct);
      if (nr*nc==0)
        nr = nc = 0;

      mt = MTP(rt,ct);
      if (MD_ROWS_IN_MTYPE(md,mt)!=nr)
        return (NO);
      if (MD_COLS_IN_MTYPE(md,mt)!=nc)
        return (NO);
    }

  return (YES);
}

/****************************************************************************/
/** \brief Check whether MATDATA_DESC and VEC_TEMPLATE x VEC_TEMPLATE match (tensor product)

   \param md			- matrix data descriptor
   \param rvt			- row vector template
   \param cvt			- col vector template

        This function checks whether a MATDATA_DESC and a VEC_TEMPLATE x VEC_TEMPLATE
        match, i.e. the number of components per type coincide in the sense that md is a
        tensor product of rvt and cvt.

        \return <ul>
        <li> YES: md and vt match </li>
        <li> NO:  md and vt do not match </li>
        </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX MDmatchesVTxVT (const MATDATA_DESC *md, const VEC_TEMPLATE *rvt, const VEC_TEMPLATE *cvt)
{
  INT rt,ct,mt,nr,nc;

  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
    {
      nr = VT_COMP(rvt,rt);
      nc = VT_COMP(cvt,ct);
      if (nr*nc==0)
        nr = nc = 0;

      mt = MTP(rt,ct);
      if (MD_ROWS_IN_MTYPE(md,mt)!=nr)
        return (NO);
      if (MD_COLS_IN_MTYPE(md,mt)!=nc)
        return (NO);
    }

  return (YES);
}

static INT MTmatchesVTxVT (const MAT_TEMPLATE *mt, const VEC_TEMPLATE *rvt, const VEC_TEMPLATE *cvt)
{
  INT rt,ct,mtp,nr,nc;

  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
    {
      nr = VT_COMP(rvt,rt);
      nc = VT_COMP(cvt,ct);
      if (nr*nc==0)
        nr = nc = 0;

      mtp = MTP(rt,ct);
      if (MT_RCOMP(mt,mtp)!=nr)
        return (NO);
      if (MT_CCOMP(mt,mtp)!=nc)
        return (NO);
    }

  return (YES);
}

/****************************************************************************/
/** \brief Create a MATDATA_DESC as a sub descriptor from a matrix template

   \param vd			- make a sub desc of this MATDATA_DESC
   \param vt			- template containing sub descriptor
   \param sub			- index of sub descriptor in template
   \param subvd		- handle to created sub descriptor

        This function creates a sub descriptor to a given MATDATA_DESC according to the given
        sub descriptor of a matrix template.

        \return <ul>
        <li>    0: ok </li>
        <li>     n: if an error occured </li>
        </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX MDsubDescFromMT (const MATDATA_DESC *md, const MAT_TEMPLATE *mt, INT sub, MATDATA_DESC **submd)
{
  FORMAT *fmt;
  SUBMAT *subm;
  SHORT *CmpsInType[NMATTYPES],SubComp[MAX_MAT_COMP];
  INT i,k,l,type,nc,nn,cmp;
  char SubName[2*MAX_MAT_COMP],buffer[NAMESIZE];
  const SHORT *offset = MD_OFFSETPTR(md);
  const SHORT *Comp = VM_COMPPTR(md);

  if (CompMDwithMT(md,mt))
    REP_ERR_RETURN(1);

  subm = MT_SUB(mt,sub);

  /* generate name and see if desc already exists */
  strcpy(buffer,SUBM_NAME(subm));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,ENVITEM_NAME(md));
  *submd = GetMatDataDescByName(MD_MG(md),buffer);
  if (*submd != NULL) {
    if (TransmitLockStatusMD(md,*submd))
      REP_ERR_RETURN(1);
    return(0);
  }
  fmt = MGFORMAT(MD_MG(md));

  /* compute sub components */
  k = 0;
  for (type=0; type<NMATTYPES; type++)
  {
    nc = SUBM_RCOMP(subm,type)*SUBM_CCOMP(subm,type);
    nn = MD_NCMPS_IN_MTYPE(md,type);
    CmpsInType[type] = &(SubComp[k]);
    for (i=0; i<nc; i++)
    {
      l = SUBM_MCMP_OF_MTYPE(subm,type,i);
      if (l>=nn)
        REP_ERR_RETURN (1);
      cmp = offset[type]+l;
      SubComp[k] = Comp[cmp];
      SubName[2*k]   = MT_COMPNAME(mt,2*cmp);
      SubName[2*k+1] = MT_COMPNAME(mt,2*cmp+1);
      k++;
    }
  }
  *submd = CreateSubMatDesc(MD_MG(md),buffer,SubName,SUBM_RCOMPS(subm),
                            SUBM_CCOMPS(subm),CmpsInType);
  if (*submd == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusMD(md,*submd))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/** \brief Create a MATDATA_DESC as a sub descriptor from a vector template

   \param md			- make a sub desc of this MATDATA_DESC
   \param vt			- template containing sub descriptor
   \param sub			- index of sub descriptor in template
   \param submd		- handle to created sub descriptor

   This function creates a sub descriptor to a given MATDATA_DESC according to the given
   subv of a vector template. The matrix is composed as tensor product of the subv.

   \return <ul>
   <li>         0: ok </li>
   <li>  n: if an error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX MDsubDescFromVT (const MATDATA_DESC *md, const VEC_TEMPLATE *vt, INT sub, MATDATA_DESC **submd)
{
  FORMAT *fmt;
  SUBVEC *subv;
  SHORT SubComp[MAX_MAT_COMP];
  SHORT RComp[NMATTYPES];
  SHORT CComp[NMATTYPES];
  INT i,j,k,l,rt,ct,mt,nr,nc,NC,nn,cmp;
  char SubName[MAX_MAT_COMP],buffer[NAMESIZE];

  PrintErrorMessageF('E',"MDsubDescFromVT","not yet implemented");
  REP_ERR_RETURN (1);

  fmt    = MGFORMAT(MD_MG(md));

  ASSERT(sub<VT_NSUB(vt));

  subv   = VT_SUB(vt,sub);

  /* generate name and see if desc already exists */
  strcpy(buffer,SUBV_NAME(subv));
  strcat(buffer,SUBV_NAME(subv));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,ENVITEM_NAME(md));
  *submd = GetMatDataDescByName(MD_MG(md),buffer);
  if (*submd != NULL) {
    if (TransmitLockStatusMD(md,*submd))
      REP_ERR_RETURN(1);
    return(0);
  }

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
          /* cmp = #@@@ */
          PrintErrorMessageF('E',"MDsubDescFromVT","not yet implemented");
          REP_ERR_RETURN (1);
          SubComp[k] = MD_MCMP_OF_MTYPE(md,mt,cmp);
          SubName[2*k]   = VM_COMP_NAME(md,2*cmp);
          SubName[2*k+1] = VM_COMP_NAME(md,2*cmp+1);
          k++;
        }
    }

  /* #@@@@
   *submd = CreateSubMatDesc(MD_MG(md),buffer,SubName,RComp,CComp,SubComp); */
  if (*submd == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusMD(md,*submd))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/** \brief Create a MATDATA_DESC as a sub descriptor from a vector template

   \param md			- make a sub desc of this MATDATA_DESC
   \param rvt			- template containing row sub descriptor
   \param rsub		- index of row sub descriptor in template
   \param cvt			- template containing col sub descriptor
   \param csub		- index of col sub descriptor in template
   \param submd		- handle to created sub descriptor

   This function creates a sub descriptor to a given MATDATA_DESC according to the given
   subv of a vector template. The matrix is composed as tensor product of the subv descriptors.

   \return <ul>
   <li>         0: ok </li>
   <li>  n: if an error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX MDsubDescFromVTxVT (const MATDATA_DESC *md, const VEC_TEMPLATE *rvt, INT rsub,
                                      const VEC_TEMPLATE *cvt, INT csub,
                                      MATDATA_DESC **submd)
{
  FORMAT *fmt;
  const VEC_TEMPLATE *vt;
  SUBVEC *rsubv,*csubv,*subv;
  SHORT SubComp[MAX_MAT_COMP];
  const SHORT *offset,*Comp;
  SHORT RComp[NMATTYPES];
  SHORT CComp[NMATTYPES];
  INT i,j,k,l,rt,ct,mt,nr,nc,NC,nn,cmp,type,n;
  char SubName[MAX_MAT_COMP],buffer[NAMESIZE];

  PrintErrorMessageF('E',"MDsubDescFromVTxVT","not yet implemented");
  return (1);

  fmt    = MGFORMAT(MD_MG(md));

  if (!MDmatchesVTxVT(md,rvt,cvt))
    REP_ERR_RETURN(1);

  if ((rsub==FULL_TPLT) && (csub==FULL_TPLT))
    REP_ERR_RETURN(1);

  if ((rsub==FULL_TPLT) || (csub==FULL_TPLT))
  {
    /* create subv identical to template */
    subv = (SUBVEC*)AllocEnvMemory(sizeof(SUBVEC));
    if (subv==NULL)
      REP_ERR_RETURN(1);

    memset(subv,0,sizeof(SUBVEC));

    if (rsub==FULL_TPLT)
    {
      vt = rvt;
      rsubv = subv;
      ASSERT(csub<VT_NSUB(cvt));
      csubv   = VT_SUB(cvt,csub);
    }
    else
    {
      vt = cvt;
      csubv = subv;
      ASSERT(rsub<VT_NSUB(rvt));
      rsubv   = VT_SUB(rvt,rsub);
    }
    for (type=0; type<NVECTYPES; type++)
    {
      n = SUBV_NCOMP(subv,type) = VT_COMP(vt,type);
      for (i=0; i<n; i++)
        SUBV_COMP(subv,type,i) = i;
    }
  }
  else
  {
    ASSERT(rsub<VT_NSUB(rvt));
    rsubv   = VT_SUB(rvt,rsub);
    ASSERT(csub<VT_NSUB(cvt));
    csubv   = VT_SUB(cvt,csub);
  }

  /* generate name and see if desc already exists */
  if (rsub==FULL_TPLT)
    strcpy(buffer,ENVITEM_NAME(rvt));
  else
    strcpy(buffer,SUBV_NAME(rsubv));
  if (csub==FULL_TPLT)
    strcpy(buffer,ENVITEM_NAME(cvt));
  else
    strcat(buffer,SUBV_NAME(csubv));
  strcat(buffer,GENERATED_NAMES_SEPERATOR);
  strcat(buffer,ENVITEM_NAME(md));
  *submd = GetMatDataDescByName(MD_MG(md),buffer);
  if (*submd != NULL) {
    if (TransmitLockStatusMD(md,*submd))
      REP_ERR_RETURN(1);
    return(0);
  }

  offset = MD_OFFSETPTR(md);
  Comp   = VM_COMPPTR(md);

  /* compute sub components */
  k  = 0;
  for (rt=0; rt<NVECTYPES; rt++)
    for (ct=0; ct<NVECTYPES; ct++)
    {
      nr = SUBV_NCOMP(rsubv,rt);
      nc = SUBV_NCOMP(csubv,ct);
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
          l = SUBV_COMP(rsubv,rt,i) * NC + SUBV_COMP(csubv,ct,j);
          if (l>=nn)
            REP_ERR_RETURN (1);
          cmp = offset[mt]+l;
          SubComp[k] = Comp[cmp];
          SubName[2*k]   = VM_COMP_NAME(md,2*cmp);
          SubName[2*k+1] = VM_COMP_NAME(md,2*cmp+1);
          k++;
        }
    }

  if ((rsub==FULL_TPLT) || (csub==FULL_TPLT))
    FreeEnvMemory(subv);
  /* #@@@
   *submd = CreateSubMatDesc(MD_MG(md),buffer,SubName,RComp,CComp,SubComp);*/
  if (*submd == NULL)
    REP_ERR_RETURN (1);
  if (TransmitLockStatusMD(md,*submd))
    REP_ERR_RETURN(1);

  return (0);
}


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
        if (VT_SUB(vt,i)!=NULL)
          FreeEnvMemory(VT_SUB(vt,i));
      VT_NSUB(vt) = 0;
    }
    else if (ENVITEM_TYPE(item) == theMatVarID)
    {
      mt = (MAT_TEMPLATE*) item;

      for (i=0; i<MT_NSUB(mt); i++)
        if (MT_SUB(mt,i)!=NULL)
          FreeEnvMemory(MT_SUB(mt,i));
      MT_NSUB(mt) = 0;
    }
  return (0);
}


/****************************************************************************/
/** \brief Remove format including sub descriptors (iff)

   \param name - format name

   Remove format including sub descriptors (iff). It is not sufficient to
   call DeleteFormat for formats allocated using the CreateFormatCmd of
   formats.c since sub descriptors are allocated directly from the environment
   heap. Calling this function cleans everything and all memory is released.

   \return
   0: ok

   \sa
   CreateFormatCmd
 */
/****************************************************************************/

INT NS_DIM_PREFIX RemoveFormatWithSubs (const char *name)
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

/**********************************************************************************/
/** \page newformat newformat - Init a format and allocate templates for vector and matrix descriptors

        The 'newformat' command enrols a format for multigrid user data.
        It also creates templates for vector and matrix descriptors.

   \section Syntax
   \verbatim
   newformat <format name>
        {$T <type specifier>}*
        {$V <dofs per type list>: <vector template name> <total needed>
                [$comp <comp names>
                        [$ident <comp identification>]
                        [$sub <sub name> <vector comp name list>]*]
        }+
        {{$M implicit(<rvt>[,<cvt>][|<sparse_matrix_name>]): <matrix template name> <total needed>
 |
         $M <matrix size list>: <matrix template name> <total needed>
                [$comp <matrix comp names>]
        }
                [$sub <sub name> {<rows>x<cols> <matrix comp name list>} | {implicit(<rsv>/<rvt>[,<csv>/<cvt>])} [$alloc <n>]]
        }+
    [$d <type name1>x<type name2><connection depth>]
    [$I <type name><mat_size>]
    [$NE]
    [$e <size>]
    [$n <size>]
   \endverbatim

        Use T-option(s) for definition of types (may be omitted, see below):
        <ul>
        <li> \<type~specifier\>   - \<type name\> in \<domain part list\>: \<object list\> </li>
        <li> \<type~name\>        - \<character\> </li>
        <li> \<domain~part~list\> - \<int\> {, \<int\>}* </li>
        <li> \<object~list\>      - \<obj\> {, \<obj\>}* </li>
        <li> \<obj\>              - nd | ed | el | si </li>
        </ul>

        NB: If no T-option is found at all it is assumed that default types are defined:
        <ul>
        <li>
        $T n in 0,...: nd $T k in 0,...: ed $T e in 0,...: el $T s in 0,...: si
        </li>
        </ul>
        to ensure downward compatibility.


        Use V-options for definition of vector templates:~
   .     <dofs~per~type~list>		- {<type name><int>}+
   .     <total~needed>			- int for total number of this template needed
                                                                        simultaneously

        Use comp-option following a V-option to specify component names:~
   .     <comp~names>				- string of single chracters, one per dof

        Use ident-option following a comp-option to specify component identification:~
   .     <comp~identification>		- string of single chracters, one per dof, multiple occurence
                                                                        of a character will invoke identification of the resp.
                                                                        components for the convergence test as well as output of
                                                                        defect and defect-reduction in solvers

        Use sub-option(s) following a comp-option to define sub templates for
        vector templates:~
   .     <vector~comp~name~list>	- any combination of characters from <comp~names>

        Use M implicit-option(s) to define matrix templates as tensor product <rvt>x<cvt>:~
   .     <rvt>						- row vector template name (defined above in a V-option)
   .     <cvt>						- col vector template name (defined above in a V-option)
        If ',<cvt' is omitted, the tensor product <rvt>x<rvt> is take.

        Use sub-option similar to those for vector templates (component names are defined in
        canonic way):~
   .     <matrix~comp~name~list>	- {<matrix comp name> }+
   .     <matrix~comp~name>		- two characters referring to the component
                                                                        names of vt1/2 (which have to be given
                                                                        above!) indicating row and col

        Alternatively specify sub-matrices by an implicit declaration deriving the sub-matrix
        as tensor product of one (or two) sub-vectors for certain vector templates:~
   .n    implicit(<rsv>/<rvt>[,<csv>/<cvt>])
   .		rsv - row sub vector name of
   .		rvt - row vector template name
   .		csv - col sub vector name of
   .		cvt - col vector template name
        If ',<csv>/<cvt>' is omitted, the tensor product of <rsv>x<rsv> is take.

        If additional storage allocation for submatrices is desired use option
   .n    alloc~<n> - allocate storage for <n> extra submatrices
        following immediately after a sub-matrix declaration.

        A second syntax for M-option(s) is still supported:~
   .     <matrix~size~list>		- {<matrix size> }+
   .     <matrix~size>				- <dofs per type list>x<dofs per type list>

        To define sub-matrix-templates as above you then have to specify component
        names with a comp-option following the M-option (second format only).

        \subsection fo Further Options
   .	  d							- specify connection depth other than 0 (inside element only)
                                                                  for <type name1>x<type name2> connections
   .	  I							- (capital i) specify interpolation matrices
   .	  NE						- require node element lists to be generated
   .	  e							- user data in elements in bytes
   .	  n							- user data in nodes in bytes

   \subsection moreInfo More Information
        The sparse_matrix_name must be a structure in
        ':SparseFormats' that contains for every combination of types a
        string 'T<type~name><type~name>' (three characters).  Since
        usually a modified sparse structure will be needed for "diagonal"
        matrices (here we mean matrices pointing to the vector itself), a
        string 'D<type~name><type~name>' has to be supplied as well.
        This string consists of '0's and '*'s meaning zero resp. non-zero
        entries.  Additionally, non-zero entries may be identified by
        using the characters 'a-z' at the appropriate places.

   \section Example
   \verbatim
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
   \endverbatim

   \section Keywords
        storage, format
 */
/**********************************************************************************/

static INT ScanVecOption (INT argc, char **argv,                        /* option list						*/
                          INT *curropt,                                 /* next option to scan				*/
                          INT po2t[][MAXVOBJECTS],          /* part-obj to type table			*/
                          INT MaxType,                                  /* bound for type id				*/
                          const char TypeNames[],                       /* names of types					*/
                          INT TypeUsed[],                                       /* indicate whether type is used	*/
                          INT *nvec,                                            /* just an index for templates		*/
                          SHORT VecStorageNeeded[])             /* to accumulate storage per type	*/
{
  VEC_TEMPLATE *vt,*vv;
  SUBVEC *subv;
  INT i,j,type,nsc[NMATTYPES];
  INT opt;
  SHORT offset[NMATOFFSETS];
  char tpltname[NAMESIZE],*names,*token,*p,tp;
  char ident[V_COMP_NAMES];
  int n;

  opt = *curropt;

  /* find name separator */
  if ((names=strchr(argv[opt],NAMESEP))==NULL)
  {
    PrintErrorMessageF('E',"newformat","separate names by a colon ':' from the description (in '$%s')",argv[opt]);
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
  if (strstr(tpltname,GENERATED_NAMES_SEPERATOR)!=NULL)
  {
    PrintErrorMessageF('E',"newformat",
                       "vector template name '%s' is not allowed to contain '%s' (in '$%s')",
                       tpltname,GENERATED_NAMES_SEPERATOR,argv[opt]);
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
      /* check uniqueness */
      for (p=VT_COMPNAMES(vt); *p!='\0'; p++)
        if (strchr(p+1,*p)!=NULL)
        {
          PrintErrorMessageF('E',"newformat",
                             "vec component names are not unique (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }

      /* check next arg for ident */
      if (opt+1 < argc)
        if (strncmp(argv[opt+1],"ident",5)==0)
        {
          opt++;
          if (sscanf(argv[opt],"ident %s",ident)!=1) {
            PrintErrorMessageF('E',"newformat",
                               "no vector comp names specified with ident option (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }
          if (strlen(ident)!=offset[NVECTYPES]) {
            PrintErrorMessageF('E',"newformat",
                               "number of ident comp names != number of comps (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }

          /* compute identification table */
          VT_NID(vt) = 0;
          for (i=0; i<offset[NVECTYPES]; i++)
            for (j=0; j<=i; j++)
              if (ident[i]==ident[j])
              {
                VT_IDENT(vt,i) = j;
                if (i==j)
                  VT_NID(vt)++;
                break;
              }
        }

      /* check next args for subv */
      while ((opt+1<argc) && (strncmp(argv[opt+1],"sub",3)==0)) {
        opt++;
        if (VT_NSUB(vt)>=MAX_SUB) {
          PrintErrorMessageF('E',"newformat",
                             "max number of vector subs exceeded (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }
        subv = (SUBVEC*)AllocEnvMemory(sizeof(SUBVEC));
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
        if (strstr(token,GENERATED_NAMES_SEPERATOR)!=NULL)
        {
          PrintErrorMessageF('E',"newformat",
                             "sub vector name '%s' is not allowed to contain '%s' (in '$%s')",
                             token,GENERATED_NAMES_SEPERATOR,argv[opt]);
          REP_ERR_RETURN (1);
        }
        strcpy(SUBV_NAME(subv),token);

        /* check uniqueness of name */
        for (i=0; i<VT_NSUB(vt)-1; i++)
          if (strcmp(SUBV_NAME(VT_SUB(vt,i)),SUBV_NAME(subv))==0)
          {
            PrintErrorMessageF('E',"newformat",
                               "subv name not unique (in '$%s')",argv[opt]);
            REP_ERR_RETURN (1);
          }

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

static INT ParseImplicitMTDeclaration (const char *str, INT MaxType,
                                       const char TypeNames[], MAT_TEMPLATE *mt)
{
  ENVDIR *dir;
  ENVITEM *item;
  VEC_TEMPLATE *rvt,*cvt;
  INT j,k,n,off,type,rtype,ctype;
  INT nr,nc;
  SHORT roffset[NMATOFFSETS],coffset[NMATOFFSETS];
  const char *p;
  char *t,tpltname[NAMESIZE];
  char txx[4];        /* field for building the sparse format name */
  STRVAR * strvar;

  /* parse row template in implicit(rvt[,cvt]) */
  if ((p=strchr(str,'('))==NULL)
    REP_ERR_RETURN(1);
  for (t=tpltname, p++; *p!='\0' && *p!=',' && *p!='|' && *p!=')'; p++)
    *(t++) = *p;
  *t = '\0';

  if ((dir=ChangeEnvDir("/newformat"))==NULL)
    REP_ERR_RETURN(2);
  for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theVecVarID)
      if (strcmp(ENVITEM_NAME(item),tpltname)==0)
        break;
  if ((rvt=(VEC_TEMPLATE *)item)==NULL)
  {
    PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                       "row vec template in '%s' not found (in '%s')",tpltname,str);
    REP_ERR_RETURN (2);
  }

  if (*p==',')
  {
    /* parse col template in implicit(rvt[,cvt]) */
    for (t=tpltname, p++; *p!='\0' && *p!='|' && *p!=')'; p++)
      *(t++) = *p;
    *t = '\0';

    if ((dir=ChangeEnvDir("/newformat"))==NULL)
      REP_ERR_RETURN(2);
    for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theVecVarID)
        if (strcmp(ENVITEM_NAME(item),tpltname)==0)
          break;
    if ((cvt=(VEC_TEMPLATE *)item)==NULL)
    {
      PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                         "col vec template in '%s' not found (in '%s')",tpltname,str);
      REP_ERR_RETURN (2);
    }
  }
  else
    cvt = rvt;

  /* define matrix template implicitly by (rvt x cvt) */
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=0; ctype<NVECTYPES; ctype++)
    {
      nr = VT_COMP(rvt,rtype);
      nc = VT_COMP(cvt,ctype);
      if (nr*nc<=0) nr=nc=0;
      type = MTP(rtype,ctype);
      MT_RCOMP(mt,type) = nr;
      MT_CCOMP(mt,type) = nc;
      if (rtype==ctype)
      {
        /* define also the size of the diagonal types */
        type = DMTP(rtype);
        MT_RCOMP(mt,type) = nr;
        MT_CCOMP(mt,type) = nc;
      }
    }

  MT_COMPNAMES(mt)[0] = '\0';
  if ((VT_COMPNAMES(rvt)[0]!=' ') && VT_COMPNAMES(cvt)[0]!=' ')
  {
    /* define also compnames */
    ConstructVecOffsets(VT_COMPS(rvt),roffset);
    ConstructVecOffsets(VT_COMPS(cvt),coffset);

    t = MT_COMPNAMES(mt);
    for (rtype=0; rtype<NVECTYPES; rtype++)
      for (ctype=0; ctype<NVECTYPES; ctype++)
      {
        nr = VT_COMP(rvt,rtype);
        nc = VT_COMP(cvt,ctype);
        for (j=0; j<nr; j++)
          for (k=0; k<nc; k++)
          {
            *(t++) = VT_COMPNAMES(rvt)[roffset[rtype]+j];
            *(t++) = VT_COMPNAMES(cvt)[coffset[ctype]+k];
          }
      }

    *t = '\0';
  }

  if (*p=='|')
  {
    /* read name of sparse matrix format structure */
    for (t=tpltname, p++; *p!='\0' && *p!=')'; p++)
      *(t++) = *p;
    *t = '\0';

    if ((dir=FindStructure(NULL,"SparseFormats"))==NULL)
    {
      PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                         ":SparseFormats does not exist");
      REP_ERR_RETURN (2);
    }
    if ((dir=FindStructure(dir,tpltname))==NULL)
    {
      PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                         ":SparseFormats:%s does not exist",tpltname);
      REP_ERR_RETURN (2);
    }

    /* for the different types read the sparse structure */
    off=0;
    for (type=0; type<NMATTYPES; type++)
    {
      rtype = MTYPE_RT(type); if (rtype>=MaxType) continue;
      ctype = MTYPE_CT(type); if (ctype>=MaxType) continue;

      /* we don't need sparsity formats, where there are no data
         in the vector template */
      if (VT_COMP(rvt,rtype)==0) continue;
      if (VT_COMP(cvt,ctype)==0) continue;

      if (type<NMATTYPES_NORMAL) txx[0] = 'T';else txx[0] = 'D';
      txx[1] = TypeNames[rtype];
      txx[2] = TypeNames[ctype];
      txx[3] = '\0';

      if ((strvar=FindStringVar(dir,txx))==NULL)
      {
        PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                           "sparse format '%s' not found",txx);
        REP_ERR_RETURN (2);
      }

      MT_MCMPPTR_OF_MTYPE(mt,type) = &(MT_COMP(mt,off));
      if ((n=MT_RCOMP(mt,type)*MT_CCOMP(mt,type)) != 0)
      {
        if (off+n>MAX_MAT_COMP_TOTAL)
        {
          PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                             "too many matrix entries per MAT_DATA_DESC");
          REP_ERR_RETURN (2);
        }

        if (String2SMArray(n,strvar->s,MT_MCMPPTR_OF_MTYPE(mt,type)) != 0)
        {
          PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                             "could not read '%s' as sparse matrix array",txx);
          REP_ERR_RETURN (2);
        }
        off+=n;
      }
    }
  }
  else
  {
    /* set standard format */
    off=0;
    for (type=0; type<NMATTYPES; type++)
    {
      MT_MCMPPTR_OF_MTYPE(mt,type) = &(MT_COMP(mt,off));
      if ((n=MT_RCOMP(mt,type)*MT_CCOMP(mt,type)) != 0)
      {
        if (off+n>MAX_MAT_COMP)
        {
          PrintErrorMessageF('E',"ParseImplicitMTDeclaration",
                             "too many matrix entries per MAT_DATA_DESC");
          REP_ERR_RETURN (2);
        }

        for (k=0; k<n; k++)
          MT_MCMP_OF_MTYPE(mt,type,k) = k;
        off+=n;
      }
    }
  }

  return (0);
}


static INT ParseImplicitSMDeclaration (const char *str, const MAT_TEMPLATE *mt, SUBMAT *subm)
{
  ENVDIR *dir;
  ENVITEM *item;
  VEC_TEMPLATE *rvt,*cvt,*vt;
  SUBVEC *rsubv,*csubv,*subv;
  INT i,j,k,type,rtype,ctype;
  INT n,nr,nc,NC,r_sub,c_sub;
  const char *p;
  char *t,tpltname[NAMESIZE],subname[NAMESIZE];

  /* parse row sub in implicit(<rsv>/<rvt>[,<csv>/<cvt>]) */
  if ((p=strchr(str,'('))==NULL)
  {
    PrintErrorMessageF('E',"ParseImplicitSMDeclaration",
                       "left bracket missing (in '%s')",str);
    REP_ERR_RETURN (2);
  }
  for (t=subname, p++; *p!='\0' && *p!=',' && *p!=')' && *p!='/'; p++)
    *(t++) = *p;
  *t = '\0';

  if (*p=='/')
  {
    /* parse row template in implicit(<rsv>/<rvt>[,<csv>/<cvt>]) */
    for (t=tpltname, p++; *p!='\0' && *p!=',' && *p!=')'; p++)
      *(t++) = *p;
    *t = '\0';
    r_sub = TRUE;
  }
  else
  {
    /* subname actually is meant as a tpltname */
    strcpy(tpltname,subname);
    r_sub = FALSE;
  }

  /* get vector template */
  if ((dir=ChangeEnvDir("/newformat"))==NULL)
    REP_ERR_RETURN(2);
  for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == theVecVarID)
      if (strcmp(ENVITEM_NAME(item),tpltname)==0)
        break;
  if ((rvt=(VEC_TEMPLATE *)item)==NULL)
  {
    PrintErrorMessageF('E',"newformat",
                       "vec template in '%s' not found (in '%s')",tpltname,str);
    REP_ERR_RETURN (2);
  }

  if (r_sub)
  {
    /* get sub vector */
    for (i=0; i<VT_NSUB(rvt); i++)
      if (strcmp(SUBV_NAME(VT_SUB(rvt,i)),subname)==0)
        break;
    if (i>=VT_NSUB(rvt))
    {
      PrintErrorMessageF('E',"ParseImplicitSMDeclaration",
                         "sub vector '%s' of template '%s' not found (in '%s')",subname,tpltname,str);
      REP_ERR_RETURN (2);
    }
    rsubv = VT_SUB(rvt,i);
  }

  if (*p==',')
  {
    /* parse col sub in implicit(<rsv>/<rvt>[,<csv>/<cvt>]) */
    for (t=subname, p++; *p!='\0' && *p!=',' && *p!=')' && *p!='/'; p++)
      *(t++) = *p;
    *t = '\0';

    if (*p=='/')
    {
      /* parse col template in implicit(<rsv>/<rvt>[,<csv>/<cvt>]) */
      for (t=tpltname, p++; *p!='\0' && *p!=',' && *p!=')'; p++)
        *(t++) = *p;
      *t = '\0';
      c_sub = TRUE;
    }
    else
    {
      /* subname actually is meant as a tpltname */
      strcpy(tpltname,subname);
      c_sub = FALSE;
    }

    /* get vector template */
    if ((dir=ChangeEnvDir("/newformat"))==NULL)
      REP_ERR_RETURN(2);
    for (item=ENVITEM_DOWN(dir); item != NULL; item = NEXT_ENVITEM(item))
      if (ENVITEM_TYPE(item) == theVecVarID)
        if (strcmp(ENVITEM_NAME(item),tpltname)==0)
          break;
    if ((cvt=(VEC_TEMPLATE *)item)==NULL)
    {
      PrintErrorMessageF('E',"newformat",
                         "col vec template in '%s' not found (in '%s')",tpltname,str);
      REP_ERR_RETURN (2);
    }

    if (c_sub)
    {
      /* get sub vector */
      for (i=0; i<VT_NSUB(cvt); i++)
        if (strcmp(SUBV_NAME(VT_SUB(cvt,i)),subname)==0)
          break;
      if (i>=VT_NSUB(cvt))
      {
        PrintErrorMessageF('E',"ParseImplicitSMDeclaration",
                           "col sub vector '%s' of col template '%s' not found (in '%s')",subname,tpltname,str);
        REP_ERR_RETURN (2);
      }
      csubv = VT_SUB(cvt,i);
    }
  }
  else
  {
    cvt = rvt;
    csubv = rsubv;
    c_sub = r_sub;
  }

  if (!r_sub && !c_sub)
  {
    PrintErrorMessageF('E',"ParseImplicitSMDeclaration",
                       "neither row nor col sub specified: matrix sub would be identical to matrix template (in '%s')",str);
    REP_ERR_RETURN (2);
  }

  /* check compatibility of vec templates with mat templates */
  if (!MTmatchesVTxVT(mt,rvt,cvt))
  {
    PrintErrorMessageF('E',"ParseImplicitSMDeclaration",
                       "row template and col template do not match matrix template (in '%s')",str);
    REP_ERR_RETURN(1);
  }

  if (!r_sub || !c_sub)
  {
    /* create subv identical to template */
    subv = (SUBVEC*)AllocEnvMemory(sizeof(SUBVEC));
    if (subv==NULL)
      REP_ERR_RETURN(1);
    memset(subv,0,sizeof(SUBVEC));

    if (!r_sub)
    {
      vt = rvt;
      rsubv = subv;
    }
    else
    {
      vt = cvt;
      csubv = subv;
    }
    for (type=0; type<NVECTYPES; type++)
    {
      n = SUBV_NCOMP(subv,type) = VT_COMP(vt,type);
      for (i=0; i<n; i++)
        SUBV_COMP(subv,type,i) = i;
    }
  }

  /* fill sub matrix template (RCOMP, CCOMP, CmpsInType, Comps) */
  k  = 0;
  for (rtype=0; rtype<NVECTYPES; rtype++)
    for (ctype=0; ctype<NVECTYPES; ctype++)
    {
      type = MTP(rtype,ctype);
      SUBM_MCMPPTR_OF_MTYPE(subm,type) = &SUBM_COMP(subm,k);
      nr = SUBV_NCOMP(rsubv,rtype);
      nc = SUBV_NCOMP(csubv,ctype);
      if (nr*nc<=0)
        nr = nc = 0;
      SUBM_RCOMP(subm,type) = nr;
      SUBM_CCOMP(subm,type) = nc;

      NC = MT_CCOMP(mt,type);
      for (i=0; i<nr; i++)
        for (j=0; j<nc; j++)
          SUBM_COMP(subm,k++) = SUBV_COMP(rsubv,rtype,i)*NC + SUBV_COMP(csubv,ctype,j);
    }

  /* and for the diagonal types (note that blocks are not necessarily quadratic) */
  for (rtype=0; rtype<NVECTYPES; rtype++)
  {
    type = DMTP(rtype);
    SUBM_MCMPPTR_OF_MTYPE(subm,type) = &(SUBM_COMP(subm,k));
    nr = SUBV_NCOMP(rsubv,rtype);
    nc = SUBV_NCOMP(csubv,rtype);
    if (nr*nc<=0)
      nr = nc = 0;
    SUBM_RCOMP(subm,type) = nr;
    SUBM_CCOMP(subm,type) = nc;

    NC = MT_CCOMP(mt,type);
    for (i=0; i<nr; i++)
      for (j=0; j<nc; j++)
        SUBM_COMP(subm,k++) = SUBV_COMP(rsubv,rtype,i)*NC + SUBV_COMP(csubv,rtype,j);
  }

  if (!r_sub || !c_sub)
    FreeEnvMemory(subv);

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
  MAT_TEMPLATE *mt,*mm;
  SUBMAT *subm;
  INT opt,i,j,checksub,type,currtype,rtype,ctype,nsc[NMATTYPES];
  SHORT offset[NMATOFFSETS];
  char tpltname[NAMESIZE],*names,*token,rt,ct,*p;
  int n,nr,nc;
  SHORT N,Nred;

  opt = *curropt;

  /* find name seperator */
  if ((names=strchr(argv[opt],NAMESEP))==NULL) {
    PrintErrorMessageF('E',"newformat",
                       "separate names by a colon ':' from the description (in '$%s')",argv[opt]);
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
  if (strstr(tpltname,GENERATED_NAMES_SEPERATOR)!=NULL)
  {
    PrintErrorMessageF('E',"newformat",
                       "matrix template name '%s' is not allowed to contain '%s' (in '$%s')",
                       tpltname,GENERATED_NAMES_SEPERATOR,argv[opt]);
    REP_ERR_RETURN (1);
  }
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
  if (token==NULL)
  {
    PrintErrorMessageF('E',"newformat","empty definition in matrix template declaration (in '$%s')",argv[opt]);
    REP_ERR_RETURN (1);
  }
  if (strncmp(token,"implicit",8)==0)
  {
    if (ParseImplicitMTDeclaration(token,MaxType,TypeNames,mt))
      REP_ERR_RETURN(1);

    ConstructMatOffsets(MT_RCOMPS(mt),
                        MT_CCOMPS(mt),offset);

    if (strlen(MT_COMPNAMES(mt))==2*offset[NMATTYPES_NORMAL])
      checksub = TRUE;
  }
  else
  {
    PrintErrorMessageF('E',"ScanMatOption",
                       "old style not yet implemented",argv[opt]);
    REP_ERR_RETURN (1);

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
        if (strlen(MT_COMPNAMES(mt))!=2*offset[NMATTYPES_NORMAL]) {
          PrintErrorMessageF('E',"newformat",
                             "number of matrix comp names != number of comps (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }
        /* check uniqueness */
        for (i=0; i<strlen(MT_COMPNAMES(mt)); i+=2)
          for (j=i+2; j<strlen(MT_COMPNAMES(mt)); j+=2)
            if (MT_COMPNAMES(mt)[i]==MT_COMPNAMES(mt)[j])
              if (MT_COMPNAMES(mt)[i+1]==MT_COMPNAMES(mt)[j+1])
              {
                PrintErrorMessageF('E',"newformat",
                                   "mat component names are not unique (in '$%s')",argv[opt]);
                REP_ERR_RETURN (1);
              }
        checksub = TRUE;
      }
  }

  if (checksub)
    /* check next args for subm */
    while ((opt+1<argc) && (strncmp(argv[opt+1],"sub",3)==0))
    {
      int ns;

      opt++;
      if (MT_NSUB(mt)>=MAX_SUB) {
        PrintErrorMessageF('E',"newformat",
                           "max number of matrix subs exceeded (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      subm = (SUBMAT*)AllocEnvMemory(sizeof(SUBMAT));
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
      if (strstr(token,GENERATED_NAMES_SEPERATOR)!=NULL)
      {
        PrintErrorMessageF('E',"newformat",
                           "sub matrix name '%s' is not allowed to contain '%s' (in '$%s')",
                           token,GENERATED_NAMES_SEPERATOR,argv[opt]);
        REP_ERR_RETURN (1);
      }
      strcpy(SUBM_NAME(subm),token);

      /* check uniqueness of name */
      for (i=0; i<MT_NSUB(mt)-1; i++)
        if (strcmp(SUBM_NAME(MT_SUB(mt,i)),SUBM_NAME(subm))==0)
        {
          PrintErrorMessageF('E',"newformat",
                             "subm name not unique (in '$%s')",argv[opt]);
          REP_ERR_RETURN (1);
        }

      /* check next token first for implicit declaration */
      if ((token=strtok(NULL,BLANKS))==NULL)
      {
        PrintErrorMessageF('E',"newformat",
                           "implicit declaration or size expected (in '$%s')",argv[opt]);
        REP_ERR_RETURN (1);
      }
      if (strncmp(token,"implicit",8)==0)
      {
        if (ParseImplicitSMDeclaration(token,mt,subm))
          REP_ERR_RETURN(1);

        /* check next arg for storage allocation */
        if (opt+1<argc)
          if (sscanf(argv[opt+1],"alloc %d",&ns)==1)
          {
            opt++;
            for (type=0; type<NMATTYPES; type++)
              MatStorageNeeded[type] += ns*SUBM_RCOMP(subm,type)*SUBM_CCOMP(subm,type);
          }
        continue;                               /* while ((opt+1<argc) && (strncmp(argv[opt+1],"sub",3)==0)) */
      }

      /* subm comps */
      PrintErrorMessageF('E',"ScanMatOption",
                         "old SUBM style not yet bugfree",argv[opt]);
      REP_ERR_RETURN (1);

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
          SUBM_MCMP_OF_MTYPE(subm,type,nsc[type]++) = n-offset[type];                               /* !!@ */
          if (nsc[type]==nr*nc) break;
        }
        SUBM_RCOMP(subm,type) = nr;
        SUBM_CCOMP(subm,type) = nc;
      }
      while ((token=strtok(NULL,BLANKS))!=NULL);

      /* check next arg for storage allocation */
      if (opt+1<argc)
        if (sscanf(argv[opt+1],"alloc %d",&ns)==1)
        {
          opt++;
          for (type=0; type<NMATTYPES; type++)
            MatStorageNeeded[type] += ns*SUBM_RCOMP(subm,type)*SUBM_CCOMP(subm,type);
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

  /* add needed storage */
  for (type=0; type<NMATTYPES; type++)
  {
    if ( ComputeSMSizeOfArray (MT_RCOMP(mt,type), MT_CCOMP(mt,type),
                               MT_MCMPPTR_OF_MTYPE(mt,type), &N, &Nred)
         !=0)
      REP_ERR_RETURN(-1);                    /* Error: ... in ComputeSMSizeOfArray */

    MatStorageNeeded[type] += n * Nred;
  }

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
  INT nparts,partlist[MAXDOMPARTS];
  INT nobjs,objlist[MAXDOMPARTS];
  char *objstr,*partstr,c,*token;
  int part;

  /* init po2t */
  for (i=0; i<MAXDOMPARTS; i++)
    for (j=0; j<MAXVOBJECTS; j++)
      po2t[i][j] = NOVTYPE;

  /* scan type specifications from option list */
  found = max = 0;
  for (opt=1; opt<argc; opt++)
    if (argv[opt][0]=='T')
    {
      if (max>=(1<<VTYPE_LEN))
      {
        printf("I would love to define another type for you, but control flags are rare... (in '$%s')",argv[opt]);
        ASSERT(FALSE);
        REP_ERR_RETURN(1);
      }

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

      /* separate object list */
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
          if (po2t[partlist[i]][objlist[j]]!=NOVTYPE)
          {
            PrintErrorMessageF('E',"newformat",
                               "the combination of obj %s in part %d is already defined (in '$%s')",ObjTypeName[objlist[j]],partlist[i],argv[opt]);
            REP_ERR_RETURN(1);
          }
          else
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

static INT CleanupTempDir (void)
{
  ENVDIR *dir;

  dir = ChangeEnvDir("/newformat");
  if (dir == NULL) {
    PrintErrorMessage('E',"CleanupTempDir","/newformat does not exist");
    REP_ERR_RETURN (1);
  }

  if (RemoveTemplateSubs((FORMAT *) dir))
    REP_ERR_RETURN (1);

  ChangeEnvDir("/");
  ENVITEM_LOCKED(dir) = 0;
  if (RemoveEnvDir((ENVITEM *) dir))
    REP_ERR_RETURN (1);

  return (0);
}

INT NS_DIM_PREFIX CreateFormatCmd (INT argc, char **argv)
{
  FORMAT *newFormat;
  ENVDIR *dir;
  VectorDescriptor vd[MAXVECTORS];
  MatrixDescriptor md[MAXMATRICES*MAXVECTORS];
  INT opt,i,j,size,type,type2,rtype,ctype,nvec,nmat,nvd,nmd;
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
  if (ScanTypeOptions(argc,argv,po2t,&MaxTypes,TypeNames)) {
    CleanupTempDir();
    REP_ERR_RETURN(1);
  }

  /* scan other options */
  for (opt=1; opt<argc; opt++)
    switch (argv[opt][0])
    {
    case 'T' :
      /* this case hase been handled before */
      break;

    case 'V' :
      if (ScanVecOption(argc,argv,&opt,po2t,MaxTypes,TypeNames,TypeUsed,&nvec,VecStorageNeeded)) {
        CleanupTempDir();
        REP_ERR_RETURN(1);
      }
      break;

    case 'M' :
      if (ScanMatOption(argc,argv,&opt,po2t,MaxTypes,TypeNames,TypeUsed,&nmat,MatStorageNeeded)) {
        CleanupTempDir();
        REP_ERR_RETURN(1);
      }
      break;

    case 'd' :
      if (ScanDepthOption(argc,argv,&opt,MaxTypes,TypeNames,TypeUsed,ConnDepth)) {
        CleanupTempDir();
        REP_ERR_RETURN(1);
      }
      break;

    case 'I' :
      if (ScanIMatOption(argc,argv,&opt,MaxTypes,TypeNames,TypeUsed,ImatTypes)) {
        CleanupTempDir();
        REP_ERR_RETURN(1);
      }
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
      CleanupTempDir();
      REP_ERR_RETURN (1);
    }

  if ((ndata == TRUE) && (nodeelementlist == TRUE)) {
    PrintErrorMessage('E',"newformat","specify either $n or $NE");
    CleanupTempDir();
    REP_ERR_RETURN (1);
  }

  /* remove types not needed from po2t */
  for (i=0; i<MAXDOMPARTS; i++)
    for (j=0; j<MAXVOBJECTS; j++)
      if (po2t[i][j] != NOVTYPE && !TypeUsed[po2t[i][j]])
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

      if (nvd > MAXVECTORS) {
        PrintErrorMessage('E',"newformat","increase MAXVECTORS");
        CleanupTempDir();
        REP_ERR_RETURN (1);
      }
    }

  if (nodeelementlist || ndata) {
    for (opt=0; opt<nvd; opt++)
      if (vd[opt].tp == NODEVEC)
        break;
    if (opt == nvd) {
      PrintErrorMessage('E',"newformat","node data requires node vector");
      CleanupTempDir();
      REP_ERR_RETURN (1);
    }
  }

  /* fill connections needed */
  nmd = 0;
  for (type=0; type<NMATTYPES; type++)
  {
    rtype = MTYPE_RT(type); ctype = MTYPE_CT(type);

    size = MatStorageNeeded[type];
    depth = ConnDepth[type];

    /***** ensure symmetry (could be circumvented)
            if (type<NMATTYPES_NORMAL)
            {
                    type2=MTP(ctype,rtype);
                    if (size<MatStorageNeeded[type2])
                            size=MatStorageNeeded[type2];
                    if (depth<ConnDepth[type2])
                            depth=ConnDepth[type2];
            } *****/

    if (ctype==rtype)
    {
      /* ensure diag/matrix coexistence (might not be necessary) */
      type2=(type<NMATTYPES_NORMAL) ? DMTP(rtype) : MTP(rtype,rtype);
      if ((size<=0) && (MatStorageNeeded[type2]<=0)) continue;
    }
    else
    {
      /* ensure symmetry of the matrix graph */
      type2=MTP(ctype,rtype);
      if ((size<=0) && (MatStorageNeeded[type2]<=0)) continue;
    }

    md[nmd].from  = rtype;
    md[nmd].to    = ctype;
    md[nmd].diag  = (type>=NMATTYPES_NORMAL);
    md[nmd].size  = size*sizeof(DOUBLE);
    md[nmd].depth = depth;
    nmd++;

    if (nmd > MAXMATRICES*MAXVECTORS) {
      PrintErrorMessage('E',"newformat","increase MAXMATRICES");
      CleanupTempDir();
      REP_ERR_RETURN (1);
    }
  }

  /* create format */
  newFormat = CreateFormat(formatname,0,0,
                           (ConversionProcPtr)NULL,(ConversionProcPtr)NULL,(ConversionProcPtr)NULL,
                           PrintTypeVectorData,PrintTypeMatrixData,
                           nvd,vd,nmd,md,ImatTypes,po2t,nodeelementlist,edata,ndata);
  if (newFormat==NULL)
  {
    PrintErrorMessage('E',"newformat","failed creating the format");
    REP_ERR_RETURN (CleanupTempDir());
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

INT NS_DIM_PREFIX InitFormats ()
{
  INT tp;

  theNewFormatDirID = GetNewEnvDirID();
  theVecVarID = GetNewEnvVarID();
  theMatVarID = GetNewEnvVarID();

  if (MakeStruct(":SparseFormats")!=0) return(__LINE__);

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
