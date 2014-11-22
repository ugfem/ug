// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  data_io.c	                                                                                                */
/*																			*/
/* Purpose:   input/output of data			                                                                */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   16.12.96 begin,												*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <config.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ugtypes.h"
#include "heaps.h"
#include "defaults.h"
#include "ugstruct.h"
#include "general.h"
#include "debug.h"
#include "ugdevices.h"
#include "gm.h"
#include "elements.h"
#include "algebra.h"
#include "misc.h"
#include "bio.h"
#include "dio.h"
#include "num.h"
#include "ugm.h"
#include "fileopen.h"

#include "data_io.h"
#include "ppif_namespace.h"

USING_UG_NAMESPACES
  USING_PPIF_NAMESPACE

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define DTIO_PARFILE    (nparfiles > 1)
#define DTIO_BLP(p,s,i) ((DTIO_BLOCK*)(((char*)(p))+((s)*(i))))

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {

  INT nElem;
  VECTOR *theV;
  DOUBLE data[1];
} DTIO_BLOCK;

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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static int ElementCompare (ELEMENT **e0, ELEMENT **e1)
{
  INT n0 = CORNERS_OF_ELEM(*e0);
  INT n1 = CORNERS_OF_ELEM(*e0);
  INT i;

  for (i=0; i<MIN(n0,n1); i++)
    if      (ID(CORNER(*e0,i))>ID(CORNER(*e1,i))) return(1);
    else if (ID(CORNER(*e0,i))<ID(CORNER(*e1,i))) return(-1);

  if (n0 == n1) return(0);

  /* error case */
  assert(0);
}

static INT LoadElementData (MULTIGRID *theMG)
{
  HEAP *Heap = MGHEAP(theMG);
  INT m = EDATA_DEF_IN_MG(theMG) / sizeof(DOUBLE);
  INT level;

  return(0);

  if (m == 0) return(0);

  for (level=0; level<=TOPLEVEL(theMG); level++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,level);
    INT i,n;
    ELEMENT **EList,*e;
    INT MarkKey;

    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e)) n++;
    MarkTmpMem(Heap,&MarkKey);
    EList = (ELEMENT**) GetTmpMem(Heap,sizeof(ELEMENT *)*n,MarkKey);
    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e)) EList[n++] = e;
    qsort(EList,n,sizeof(ELEMENT*),(int (*)(const void *, const void *))ElementCompare);
    for (i=0; i<n; i++)
      if (Bio_Read_mdouble(m,(double*)EDATA(EList[i])))
        return (1);
    ReleaseTmpMem(Heap,MarkKey);
  }
    #ifdef ModelP
  a_elementdata_consistent(theMG,0,TOPLEVEL(theMG));
    #endif

  return(0);
}

static INT SaveElementData (MULTIGRID *theMG)
{
  HEAP *Heap = MGHEAP(theMG);
  INT m = EDATA_DEF_IN_MG(theMG) / sizeof(DOUBLE);
  INT level;

  if (m == 0) return(0);

  for (level=0; level<=TOPLEVEL(theMG); level++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,level);
    INT i,n;
    ELEMENT **EList,*e;
    INT MarkKey;

    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e)) n++;
    MarkTmpMem(Heap,&MarkKey);
    EList = (ELEMENT**) GetTmpMem(Heap,sizeof(ELEMENT *)*n,MarkKey);
    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e)) EList[n++] = e;
    qsort(EList,n,sizeof(ELEMENT*),(int (*)(const void *, const void *))ElementCompare);
    for (i=0; i<n; i++)
      if (Bio_Write_mdouble(m,(double*)EDATA(EList[i])))
        return (1);
    ReleaseTmpMem(Heap,MarkKey);
  }
    #ifdef ModelP
  a_elementdata_consistent(theMG,0,TOPLEVEL(theMG));
    #endif

  return(0);
}

static int EdgeCompare (EDGE **ed0, EDGE **ed1)
{
  INT n0 = ID(NBNODE(LINK0(*ed0)));
  INT n1 = ID(NBNODE(LINK0(*ed1)));

  if (n0 < n1)
    return(1);
  else if (n0 > n1)
    return(-1);

  if (ID(NBNODE(LINK1(*ed0))) < ID(NBNODE(LINK1(*ed1))))
    return(1);

  return(-1);
}

static INT LoadEdgeData (MULTIGRID *theMG, VECDATA_DESC *v)
{
  HEAP *Heap = MGHEAP(theMG);
  INT m = VD_NCMPS_IN_TYPE(v,EDGEVEC);
  INT level;

  if (m == 0) return(0);

  for (level=0; level<=TOPLEVEL(theMG); level++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,level);
    INT i,n;
    ELEMENT *e;
    EDGE **ed;
    INT MarkKey;
    INT c = VD_CMP_OF_TYPE(v,EDGEVEC,0);

    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      for (i=0; i<EDGES_OF_ELEM(e); i++)
        /*  if (ID(CORNER(e,CORNER_OF_EDGE(e,i,0)))
                      < ID(CORNER(e,CORNER_OF_EDGE(e,i,1)))) */
        n++;
    MarkTmpMem(Heap,&MarkKey);
    ed = (EDGE**) GetTmpMem(Heap,sizeof(EDGE *)*n,MarkKey);
    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      for (i=0; i<EDGES_OF_ELEM(e); i++)
        /* if (ID(CORNER(e,CORNER_OF_EDGE(e,i,0)))
                    < ID(CORNER(e,CORNER_OF_EDGE(e,i,1)))) */
        ed[n++] =
          GetEdge(CORNER(e,CORNER_OF_EDGE(e,i,0)),
                  CORNER(e,CORNER_OF_EDGE(e,i,1)));
    qsort(ed,n,sizeof(EDGE*),
          (int (*)(const void *, const void *))EdgeCompare);
    for (i=0; i<n; i++)
      if (Bio_Read_mdouble(m,(double*)VVALUEPTR(EDVECTOR(ed[i]),c)))
        return (1);
    ReleaseTmpMem(Heap,MarkKey);
        #ifdef ModelP
    l_ghostvector_consistent(theGrid,v);
        #endif
  }

  return(0);
}

static INT SaveEdgeData (MULTIGRID *theMG, VECDATA_DESC *v)
{
  HEAP *Heap = MGHEAP(theMG);
  INT m = VD_NCMPS_IN_TYPE(v,EDGEVEC);
  INT level;

  if (m == 0) return(0);

  for (level=0; level<=TOPLEVEL(theMG); level++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,level);
    INT i,n;
    ELEMENT *e;
    EDGE **ed;
    INT MarkKey;
    INT c = VD_CMP_OF_TYPE(v,EDGEVEC,0);

    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      for (i=0; i<EDGES_OF_ELEM(e); i++)
        /*   if (ID(CORNER(e,CORNER_OF_EDGE(e,i,0)))
                       < ID(CORNER(e,CORNER_OF_EDGE(e,i,1)))) */
        n++;
    MarkTmpMem(Heap,&MarkKey);
    ed = (EDGE**) GetTmpMem(Heap,sizeof(EDGE *)*n,MarkKey);
    for (n=0, e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      for (i=0; i<EDGES_OF_ELEM(e); i++)
        /* if (ID(CORNER(e,CORNER_OF_EDGE(e,i,0)))
                      < ID(CORNER(e,CORNER_OF_EDGE(e,i,1)))) */
        ed[n++] =
          GetEdge(CORNER(e,CORNER_OF_EDGE(e,i,0)),
                  CORNER(e,CORNER_OF_EDGE(e,i,1)));
    qsort(ed,n,sizeof(EDGE*),
          (int (*)(const void *, const void *))EdgeCompare);
    for (i=0; i<n; i++)
      if (Bio_Write_mdouble(m,(double*)VVALUEPTR(EDVECTOR(ed[i]),c)))
        return (1);
    ReleaseTmpMem(Heap,MarkKey);
  }

  return(0);
}

/****************************************************************************/
/*D
   LoadData - reads data

   SYNOPSIS:
   INT LoadData (MULTIGRID *theMG, char *FileName, INT n, VECDATA_DESC **theVDList);

   PARAMETERS:
   .  theMG - multigrid
   .  FileName - name of file
   .  n - nb of vecdata desc
   .  theVDList - list of vecdata

   DESCRIPTION:
   loads data from file

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 when error occured.

   SEE ALSO:
   D*/
/****************************************************************************/

MULTIGRID * NS_DIM_PREFIX OpenMGFromDataFile (MULTIGRID *theMG, INT number, char *type, char *DataFileName, MEM heapSize)
{
  MULTIGRID *mg;
  DIO_GENERAL dio_general;
  char FileName[NAMESIZE],*mgname,*mgtype,NumberString[8],*p;
  INT close,load,nparfiles;
  char buf[64];

  if (me == master)
  {
    /* open file */
    strcpy(FileName,DataFileName);
    if (number!=-1)
    {
      sprintf(NumberString,".%06d",(int)number);
      strcat(FileName,NumberString);
    }
    strcat(FileName,".ug.data.");
    strcat(FileName,type);
    nparfiles = 1;
    if (DTIO_filetype(FileName) == FT_DIR)
    {
      sprintf(buf,"/data.%04d",(int)me);
      strcat(FileName,buf);
      if (Read_OpenDTFile (FileName))                                                                         {nparfiles = -1;}
      else
      if (Read_DT_General(&dio_general))                                                              {nparfiles = -1;}
      nparfiles = dio_general.nparfiles;
      if (nparfiles>procs)                                                                                            {UserWrite("ERROR: too many processors needed\n"); nparfiles = -1;}
      assert(dio_general.me == me);
    }
    else if(DTIO_filetype(FileName) == FT_FILE)
    {
      if (Read_OpenDTFile (FileName))                                                                         {nparfiles = -1;}
      else
      if (Read_DT_General(&dio_general))                                                              {nparfiles = -1;}
    }
    else
      nparfiles = -1;

    CloseDTFile();
  }

  Broadcast(&nparfiles,sizeof(INT));
  if (nparfiles == -1) return(NULL);

  Broadcast(&dio_general,sizeof(DIO_GENERAL));

  if (theMG==NULL)
  {
    close = 0;
    load = 1;
  }
  else if (!MG_SAVED(theMG) || dio_general.magic_cookie!=MG_MAGIC_COOKIE(theMG))
  {
    close = 1;
    load = 1;
  }
  else
  {
    close = 0;
    load = 0;
    mg = theMG;
  }

  /* dispose multigrid */
  if (close)
    if (DisposeMultiGrid(theMG))
      return (NULL);

  /* reload multigrid */
  if (load)
  {
    p = strstr(dio_general.mgfile,".ug.mg.");
    if (p==NULL) return (NULL);
    else
    {
      p[0] = '\0';
      p+=7;
      p[3]='\0';
      mgtype = p;
      mgname = dio_general.mgfile;
    }
    mg = LoadMultiGrid (NULL,mgname,mgtype,NULL,NULL,heapSize,0,0,0);
  }

  return (mg);
}

INT NS_DIM_PREFIX LoadData (MULTIGRID *theMG, char *name, char *type, INT number, INT n, VECDATA_DESC **theVDList)
{
  INT i,*entry,copied_until,copy_until,still_to_read,read,nvec,id,nparfiles;
  unsigned long ncomp, j, m, s;
  DIO_GENERAL dio_general;
#ifdef ModelP
  DIO_GENERAL temp_dio_general;
#endif
  HEAP *theHeap;
  double *data;
  VECTOR *theV, **VectorList;
  GRID *theGrid;
  NODE *theNode;
  SHORT *cp[DIO_VDMAX];
  INT ncmp[DIO_VDMAX];
  INT MarkKey;
  char FileName[NAMESIZE], NumberString[8];
  char buf[64];

  if (theMG==NULL) return (1);
  theHeap = MGHEAP(theMG);
  if (n<1) return (1);

  /* check if only NODEVECTORs and thereby get components */
  for (i=0; i<n; i++)
  {
    if (theVDList[i]==NULL) continue;
    cp[i] = VD_ncmp_cmpptr_of_otype(theVDList[i],NODEVEC,ncmp+i);
    if (ncmp[i]<=0)
    {
      PrintErrorMessageF('E',"LoadData","vd mismatch for io (no %d)",i);
      REP_ERR_RETURN(1);
    }
  }

  /* open file */
  strcpy(FileName,name);
  if (number!=-1)
  {
    sprintf(NumberString,".%06d",(int)number);
    strcat(FileName,NumberString);
  }
  strcat(FileName,".ug.data.");
  strcat(FileName,type);
#ifdef ModelP
  if (me == master)
  {
#endif
  nparfiles = 1;
  if (DTIO_filetype(FileName) == FT_DIR)
  {
    sprintf(buf,"/data.%04d",(int)me);
    strcat(FileName,buf);
    if (Read_OpenDTFile (FileName))                                                                         {nparfiles = -1;}
    else
    if (Read_DT_General(&dio_general))                                                              {CloseDTFile (); nparfiles = -1;}
    nparfiles = dio_general.nparfiles;
    if (nparfiles>procs)                                                                                            {UserWrite("ERROR: too many processors needed\n");CloseDTFile (); nparfiles = -1;}
    assert(dio_general.me == me);
  }
  else if(DTIO_filetype(FileName) == FT_FILE)
  {
    if (Read_OpenDTFile (FileName))                                                                         {nparfiles = -1;}
    else
    if (Read_DT_General(&dio_general))                                                              {CloseDTFile (); nparfiles = -1;}
  }
  else
    nparfiles = -1;
#ifdef ModelP
  Broadcast(&nparfiles,sizeof(int));
}
else
{
  Broadcast(&nparfiles,sizeof(int));
  if (me < nparfiles)
  {
    sprintf(buf,"/data.%04d",(int)me);
    strcat(FileName,buf);
    if (Read_OpenDTFile (FileName))                                                                         {nparfiles = -1;}
    else
    if (Read_DT_General(&dio_general))                                                              {CloseDTFile (); nparfiles = -1;}
  }
}
nparfiles = UG_GlobalMinINT(nparfiles);
#endif
  if (nparfiles == -1) return(1);

#ifdef ModelP
  if (procs>nparfiles)
  {
    /* for me >nparfiles do: dio_general{me}:=dio_general{master} */

    if (me==master)
      temp_dio_general = dio_general;
    Broadcast(&temp_dio_general,sizeof(DIO_GENERAL));                   /* for ease of communication use Broadcast and avoid overwriting of dio_general for me<nparfiles */
    if (me >= nparfiles)
      dio_general = temp_dio_general;
  }
#endif

  if (SetStringValue(":IO:TIME",dio_general.time))                {CloseDTFile(); return (1);}
  if (SetStringValue(":IO:DT",dio_general.dt))                    {CloseDTFile(); return (1);}
  if (SetStringValue(":IO:NDT",dio_general.ndt))                  {CloseDTFile(); return (1);}
  if (SetStringVar(":IDENTIFICATION",dio_general.ident))  {CloseDTFile(); return (1);}
  if (me >= nparfiles) return (0);

  /* read general information */
  if (strcmp(dio_general.version,DIO_VERSION)!=0)                 {CloseDTFile(); UserWrite("ERROR: wrong version\n"); return (1);}
  if (0 && dio_general.magic_cookie != MG_MAGIC_COOKIE(theMG))    {CloseDTFile(); UserWrite("m-c-error"); return (1);}
  ncomp = 0;
  for (i=0; i<dio_general.nVD; i++)
  {
    if (i<n && theVDList[i]!=NULL)
      if (dio_general.VDncomp[i] != ncmp[i]) {CloseDTFile(); UserWrite("vd-comp do not match\n"); return (1);}
    ncomp += dio_general.VDncomp[i];
  }

  /* create list of vectors */
  MarkTmpMem(theHeap,&MarkKey);
  nvec=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
    nvec += NN(GRID_ON_LEVEL(theMG,i));
  if( nvec == 0 )
  {
    CloseDTFile();
    return (0);
  }
  VectorList = (VECTOR **)GetTmpMem(theHeap,nvec*sizeof(VECTOR*),MarkKey);
  if (VectorList==NULL)                                                                   {CloseDTFile(); UserWrite("ERROR: cannot allocate memoriy for VectorList\n"); return (1);}
  for (i=0; i<nvec; i++) VectorList[i] = NULL;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      theV = NVECTOR(theNode);
      id = ID(theNode);
      if (id >= nvec)                                                                 {CloseDTFile(); UserWrite("internal ERROR: id is out of range\n"); return (1);}
      if (VectorList[id]!=NULL)
      {CloseDTFile(); UserWrite("internal ERROR: VectorList is set twice\n"); return (1);}
      VectorList[id] = theV;
    }
  }

  /* load data */
  entry = (INT *)GetTmpMem(theHeap,ncomp*sizeof(INT),MarkKey);
  if (entry==NULL)                                                                                {CloseDTFile(); return (1);}
  s=0;
  for (i=0; i<dio_general.nVD; i++)
  {
    if (i<n && theVDList[i]!=NULL)
    {
      for (j=0; j<dio_general.VDncomp[i]; j++)
        entry[s++] = cp[i][j];
    }
    else
      for (j=0; j<dio_general.VDncomp[i]; j++)
        entry[s++] = -1;
  }
  if (s!=ncomp)                                                                                   {CloseDTFile(); return (1);}
  m = theHeap->size - theHeap->used;
  m = m - m%sizeof(double); m /= sizeof(double);
  while(m>=ncomp)
  {
    data = (double *)GetTmpMem(theHeap,m*sizeof(double),MarkKey);
    if (data!=NULL) break;
    m *= 0.5;
  }
  if (data==NULL)                                                                                 {CloseDTFile(); UserWrite("ERROR: cannot allocate tmp mem\n"); return (1);}
  if (m<ncomp)                                                                                    {CloseDTFile(); UserWrite("ERROR: tmp mem too small\n"); return (1);}

  copied_until = copy_until=0; still_to_read=dio_general.ndata;
  for (id=0; id<nvec; id++)
  {
    if (copied_until>=copy_until)
    {
      if (still_to_read<=m) read=still_to_read;
      else read=m;
      read -= read%ncomp;
      if (Bio_Read_mdouble(read,data))                        {CloseDTFile(); UserWrite("ERROR: tmp mem too small\n"); return (1);}
      still_to_read -= read;
      copy_until += read;
      s=0;
    }
    theV = VectorList[id];
    for (j=0; j<ncomp; j++)
      if (entry[j]>=0)
        VVALUE(theV,entry[j]) = data[s++];
      else
        s++;
    copied_until += ncomp;
  }
  ReleaseTmpMem(theHeap,MarkKey);

  /* load the data field of the elements */
  if (EDATA_DEF_IN_MG(theMG))
    if (LoadElementData(theMG))
    {CloseDTFile(); UserWrite("ERROR: load element data failed\n"); return (1);}

  /* load edge vector data (if exists) */
  for (i=0; i<n; i++)
  {
    if (theVDList[i]==NULL) continue;
    if (VD_NCMPS_IN_TYPE(theVDList[i],EDGEVEC) > 0)
      if (LoadEdgeData(theMG,theVDList[i])) {
        CloseDTFile();
        UserWrite("ERROR: save element data failed\n");
        return (1);
      }
  }

  /* close file */
  if (CloseDTFile()) return (1);

  return (0);
}

/****************************************************************************/
/*D
   SaveData - writes data

   SYNOPSIS:
   INT SaveData (MULTIGRID *theMG, char *FileName, INT n, VECDATA_DESC **theVDList);

   PARAMETERS:
   .  theMG - multigrid
   .  FileName - name of file
   .  n - nb of vecdata desc
   .  theVDList - list of vecdata

   DESCRIPTION:
   saves data on file

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 when error occured.

   SEE ALSO:
   D*/
/****************************************************************************/

INT NS_DIM_PREFIX SaveData (MULTIGRID *theMG, char *name, INT rename, INT save_without_mg, char *type, INT number, DOUBLE time, DOUBLE dt, DOUBLE ndt, INT n, VECDATA_DESC **theVDList, EVALUES **theEVal, EVECTOR **theEVec, char **NameList)
{
  INT i,j,k,l,ncomp,s,t,*entry,nNode,store_from_eval,id,tag,coe,q,mode,nparfiles;
  DIO_GENERAL dio_general;
  HEAP *theHeap;
  GRID *theGrid;
  NODE *theNode;
  ELEMENT *theElement;
  const DOUBLE *x[MAX_CORNERS_OF_ELEM];
  DOUBLE value,fnblock;
  DOUBLE_VECTOR vector;
  char *p,FileName[NAMESIZE],NumberString[8];
  SHORT *cp[DIO_VDMAX];
  INT ncmp[DIO_VDMAX],blocksize,free,fb,lb,nblock,saved;
  INT MarkKey;
  DTIO_BLOCK *block,*bptr;
#ifdef ModelP
  INT error;
  char buf[64];
#endif

  /* init */
  if (theMG==NULL) return (1);
  if (save_without_mg) {
    if (RenumberMultiGrid(theMG,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0)!=GM_OK) {
      UserWrite("ERROR: cannot renumber multigrid\n");
      return (1);
    }
  } else {
    saved = MG_SAVED(theMG);
#ifdef ModelP
    saved = UG_GlobalMinINT(saved);
#endif
    if (!saved)
    {
      if (SaveMultiGrid (theMG,NULL,NULL,NULL,1,rename))
      {
        UserWrite("ERROR: cannot autosave multigrid\n");
        return (1);
      }
      if (!MG_SAVED(theMG))
      {
        UserWrite("ERROR: autosave of multigrid failed\n");
        return (1);
      }
    }
  }
  theHeap = MGHEAP(theMG);
  if (n<1) return (1);

  /* check if only NODEVECTORs and thereby get components */
  for (i=0; i<n; i++)
  {
    if (theVDList[i]==NULL) continue;
    cp[i] = VD_ncmp_cmpptr_of_otype(theVDList[i],NODEVEC,ncmp+i);
    if (ncmp[i]<=0)
    {
      PrintErrorMessageF('E',"SaveData","vd mismatch for io (no %d)",i);
      REP_ERR_RETURN(1);
    }
  }

  /* open file */
  nparfiles = procs;
  if (strcmp(type,"xdr")==0) mode = BIO_XDR;
  else if (strcmp(type,"asc")==0) mode = BIO_ASCII;
  else if (strcmp(type,"bin")==0) mode = BIO_BIN;
  else return (1);
  strcpy(FileName,name);
  if (number!=-1)
  {
    sprintf(NumberString,".%06d",(int)number);
    strcat(FileName,NumberString);
  }
  strcat(FileName,".ug.data.");
  strcat(FileName,type);
#ifdef ModelP
  error = 0;
  if (me == master)
    if (DTIO_PARFILE)
      if (DTIO_dircreate(FileName,(int)rename))
        error = -1;

  Broadcast(&error,sizeof(int));
  if (error == -1)
  {
    UserWriteF("SaveData(): error during file/directory creation\n");
    return(1);
  }
  if (DTIO_PARFILE)
  {
    sprintf(buf,"/data.%04d",(int)me);
    strcat(FileName,buf);
  }
#endif
  if (Write_OpenDTFile (FileName,(int)rename)) return (1);

  /* write general information */
  dio_general.mode = mode;
  strcpy(dio_general.version,DIO_VERSION);
  if (save_without_mg) strcpy(dio_general.mgfile,"saved_without_mg");
  else strcpy(dio_general.mgfile,MG_FILENAME(theMG));
  p = GetStringVar (":IDENTIFICATION");
  if (p!=NULL) strcpy(dio_general.ident,p);
  else strcpy(dio_general.ident,"---");
  if (number!=-1)
  {
    dio_general.time = time;
    dio_general.dt = dt;
    dio_general.ndt = ndt;
  }
  else
  {
    dio_general.time = -1.0;
    dio_general.dt = -1.0;
    dio_general.ndt = -1.0;
  }
  dio_general.nparfiles           = procs;
  dio_general.me                  = me;
  dio_general.magic_cookie = MG_MAGIC_COOKIE(theMG);
  dio_general.nVD = n;
  ncomp = 0; store_from_eval = 0;
  for (i=0; i<n; i++)
  {
    if (theVDList[i]!=NULL)
    {
      if (NameList==NULL) strcpy(dio_general.VDname[i],ENVITEM_NAME(theVDList[i]));
      else strcpy(dio_general.VDname[i],NameList[i]);
      ncomp += dio_general.VDncomp[i] = ncmp[i];
      if (dio_general.VDncomp[i]==1)
        dio_general.VDtype[i] = DIO_SCALAR;
      else
        dio_general.VDtype[i] = DIO_MULTIPLE_SCALAR;
      for (j=0; j<dio_general.VDncomp[i]; j++)
        dio_general.VDcompNames[i][j] = theVDList[i]->compNames[j];
      dio_general.VDcompNames[i][j] = '\0';
    }
    else if (theEVal[i]!=NULL)
    {
      if (NameList==NULL) strcpy(dio_general.VDname[i],ENVITEM_NAME(theEVal[i]));
      else strcpy(dio_general.VDname[i],NameList[i]);
      dio_general.VDncomp[i] = 1;
      ncomp += 1;
      store_from_eval = 1;
      dio_general.VDtype[i] = DIO_SCALAR;
      strcpy(dio_general.VDcompNames[i],"u");
    }
    else if (theEVec[i]!=NULL)
    {
      if (NameList==NULL) strcpy(dio_general.VDname[i],ENVITEM_NAME(theEVec[i]));
      else strcpy(dio_general.VDname[i],NameList[i]);
      dio_general.VDncomp[i] = DIM;
      ncomp += DIM;
      store_from_eval = 1;
      dio_general.VDtype[i] = DIO_VECTOR;
      strcpy(dio_general.VDcompNames[i],"xyz");
    }
    else
    {CloseDTFile(); return (1);}
  }
  nNode=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    nNode += NN(theGrid);
  }
  dio_general.ndata = nNode*ncomp;
  if (Write_DT_General (&dio_general))                                    {CloseDTFile(); return (1);}

  /* save data: entries */
  MarkTmpMem(theHeap,&MarkKey);
  entry = (INT *)GetTmpMem(theHeap,ncomp*sizeof(INT),MarkKey);
  if (entry==NULL)                                                                                {CloseDTFile(); return (1);}
  s=0;
  for (i=0; i<n; i++)
  {
    if (theVDList[i]!=NULL)
      for (j=0; j<dio_general.VDncomp[i]; j++)
        entry[s++] = cp[i][j];
    else if (theEVal[i]!=NULL)
      entry[s++] = -1;
    else if (theEVec[i]!=NULL)
      for (j=0; j<DIM; j++)
        entry[s++] = -1;
    else
      return (1);
  }
  if (s!=ncomp)                                                                                   {CloseDTFile(); return (1);}

  /* prepare eval/evec and call each, to give a change to allocate tmpmem */
  theGrid = GRID_ON_LEVEL(theMG,0);
  for (k=0; k<n; k++)
  {
    if (theVDList[k]!=NULL) continue;
    if (theEVal[k]!=NULL)
    {
      if (theEVal[k]->PreprocessProc!=NULL)
        if ((*(theEVal[k]->PreprocessProc))(ENVITEM_NAME(theEVal[k]),theMG))
          return (1);

      for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      {
        coe = CORNERS_OF_ELEM(theElement);
        for (l=0; l<coe; l++)
          x[l] = CVECT(MYVERTEX(CORNER(theElement,l)));
        (*(theEVal[k]->EvalProc))(theElement,x,LOCAL_COORD_OF_ELEM(theElement,0));
      }
    }
    else if (theEVec[k]!=NULL)
    {
      if (theEVec[k]->PreprocessProc!=NULL)
        if ((*(theEVec[k]->PreprocessProc))(ENVITEM_NAME(theEVec[k]),theMG))
          return (1);

      for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      {
        coe = CORNERS_OF_ELEM(theElement);
        for (l=0; l<coe; l++)
          x[l] = CVECT(MYVERTEX(CORNER(theElement,l)));
        (*(theEVec[k]->EvalProc))(theElement,x,LOCAL_COORD_OF_ELEM(theElement,0),vector);
      }
    }
    else
    {CloseDTFile(); return (1);}
  }

  /* save data: temporary heap */
  blocksize = sizeof(DTIO_BLOCK) + (ncomp-1)*sizeof(DOUBLE);
  free = theHeap->size - theHeap->used - 1024;
  if (free<=0)                                                                                    {CloseDTFile(); return (1);}
  fnblock = (DOUBLE)free/(DOUBLE)blocksize;
  nblock = (INT)floor(fnblock);
  if (nblock<1)                                                                                   {CloseDTFile(); return (1);}
  block = (DTIO_BLOCK *)GetTmpMem(theHeap,nblock*blocksize,MarkKey);
  if (block==NULL)                                                                                {CloseDTFile(); return (1);}
  if (100*nblock<nNode)
    UserWrite("WARNING: save will take long due to lack of temporary memory\n");

  /* save data */
  for (fb=0,lb=MIN(nblock,nNode); fb<nNode; fb+=nblock,lb=MIN(lb+nblock,nNode))
  {
    /* init blocks */
    for (i=0; i<=TOPLEVEL(theMG); i++)
    {
      theGrid = GRID_ON_LEVEL(theMG,i);
      for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
        if (ID(theNode)>=fb && ID(theNode)<lb)
        {
          bptr = DTIO_BLP(block,blocksize,ID(theNode)-fb);
          bptr->nElem = 0;
          bptr->theV = NVECTOR(theNode);
        }
    }
    for (i=0; i<=TOPLEVEL(theMG); i++)
    {
      theGrid = GRID_ON_LEVEL(theMG,i);
      for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      {
        tag = TAG(theElement);
        coe = CORNERS_OF_ELEM(theElement);
        for (l=0; l<coe; l++)
        {
          id = ID(CORNER(theElement,l));
          if (id<fb || id>=lb) continue;
          bptr = DTIO_BLP(block,blocksize,id-fb);
          bptr->nElem++;
        }
      }
    }

    /* save data from symbols */
    for (i=0; i<lb-fb; i++)
    {
      bptr = DTIO_BLP(block,blocksize,i);
      for (j=0; j<ncomp; j++)
        if (entry[j]>=0.0)
          bptr->data[j] = VVALUE(bptr->theV,entry[j]);
        else
          bptr->data[j] = 0.0;
    }

    /* save data from eval/evec-procs */
    t=0;
    for (k=0; k<n; k++)
    {
      if (theVDList[k]!=NULL) continue;
      while(entry[t]>=0.0) t++;
      if (theEVal[k]!=NULL)
      {
        /* save data from eval procs */
        for (i=0; i<=TOPLEVEL(theMG); i++)
        {
          theGrid = GRID_ON_LEVEL(theMG,i);
          for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
          {
            tag = TAG(theElement);
            coe = CORNERS_OF_ELEM(theElement);
            for (l=0; l<coe; l++)
              x[l] = CVECT(MYVERTEX(CORNER(theElement,l)));
            for (l=0; l<coe; l++)
            {
              id = ID(CORNER(theElement,l));
              if (id<fb || id>=lb) continue;
              bptr = DTIO_BLP(block,blocksize,id-fb);
              value = (*(theEVal[k]->EvalProc))(theElement,x,LOCAL_COORD_OF_ELEM(theElement,l));
              bptr->data[t] += value/(DOUBLE)(bptr->nElem);
            }
          }
        }
        t++;
      }
      else if (theEVec[k]!=NULL)
      {
        /* save data from evec procs */
        for (i=0; i<=TOPLEVEL(theMG); i++)
        {
          theGrid = GRID_ON_LEVEL(theMG,i);
          for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
          {
            tag = TAG(theElement);
            coe = CORNERS_OF_ELEM(theElement);
            for (l=0; l<coe; l++)
              x[l] = CVECT(MYVERTEX(CORNER(theElement,l)));
            for (l=0; l<coe; l++)
            {
              id = ID(CORNER(theElement,l));
              if (id<fb || id>=lb) continue;
              bptr = DTIO_BLP(block,blocksize,id-fb);
              (*(theEVec[k]->EvalProc))(theElement,x,LOCAL_COORD_OF_ELEM(theElement,l),vector);
              for (q=0; q<DIM; q++)
                bptr->data[t+q] += vector[q]/(DOUBLE)(bptr->nElem);
            }
          }
        }
        t += DIM;
      }
      else                                                                                                                                                                    {CloseDTFile(); return (1);}
    }

    /* save data to file */
    for (i=0; i<lb-fb; i++)
    {
      bptr = DTIO_BLP(block,blocksize,i);
      if (Bio_Write_mdouble(ncomp,bptr->data))                                                                                                {CloseDTFile(); return (1);}
    }
  }
  ReleaseTmpMem(theHeap,MarkKey);

  /* save the data field of the elements */
  if (EDATA_DEF_IN_MG(theMG))
    if (SaveElementData(theMG))
    {CloseDTFile(); UserWrite("ERROR: save element data failed\n"); return (1);}

  /* save edge vector data (if exists) */
  for (i=0; i<n; i++)
  {
    if (theVDList[i]==NULL) continue;
    if (VD_NCMPS_IN_TYPE(theVDList[i],EDGEVEC) > 0)
      if (SaveEdgeData(theMG,theVDList[i])) {
        CloseDTFile();
        UserWrite("ERROR: save element data failed\n");
        return (1);
      }
  }

  /* close file */
  if (CloseDTFile()) return (1);

  return (0);
}
