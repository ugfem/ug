// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugpfm.c														*/
/*																			*/
/* Purpose:   parallel file merger											*/
/*																			*/
/* Author:	  Klaus Johannsen/Stefan Lnag                                                                   */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   12.05.98 begin, version 1.0									*/
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "compiler.h"
#include "fileopen.h"
#include "heaps.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"
#include "bio.h"
#include "ugstruct.h"

#include "devices.h"

#include "gm.h"
#include "algebra.h"
#include "misc.h"
#include "ugm.h"
#include "ugio.h"
#include "elements.h"
#include "shapes.h"
#include "mgio.h"
#include "dio.h"


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define HASH_MULTIPLICATOR                      172742983
#define ELEMPROCLISTSIZE                2000
#define PROCLISTSIZE                    (ELEMPROCLISTSIZE*MAX_SONS * MAX(5,(int)(2.0+log((double)nparfiles))))

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {

  void *obj;                                            /* object                                                                       */
  int lid;                                              /* local id created automatically			*/
  int key[1];                       /* multiple key								*/
} ENTRY;

typedef struct {
  int n_obj;                                    /* current number of objects				*/
  int table_len;                                /* max number of objects					*/
  int lock;                                             /* lock for push, pop						*/
  int next_get;                                 /* next entry returned by 'get'				*/
  int next_lid;
  ENTRY *entry[1];                      /* table used								*/
} HASH_TABLE;


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

static int nparfiles;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

int Get_h(double *in, double *out) {
  return (0);
}
int AllMemElements(int nelements){
  return (0);
}
int AllMemInnerPoints(int npoints){
  return (0);
}
int AddInnerNode(double x, double y, double z){
  return (0);
}
int AddTetrahedron (int node0, int node1, int node2, int node3){
  return (0);
}
int AddElement (int node0, int node1, int node2, int node3) {
  return (0);
}

/****************************************************************************/
/****************************** ugpfm-low ***********************************/
/****************************************************************************/

HASH_TABLE *CreateHashTable (int max_obj)
{
  int i,j,table_len,max_factor;
  HASH_TABLE *ht;

  /* determine table_len */
  max_factor=1;
  for (i=2*max_obj-1; i<2*max_obj+200; i+=2)
    for (j=1; j<sqrt(i); j+=2)
      if (i%j==0)
      {
        if (j>max_factor)
        {
          max_factor=j;
          table_len=i;
        }
        break;
      }
  ht = (HASH_TABLE *)malloc(sizeof(HASH_TABLE)+(table_len-1)*sizeof(ENTRY*));
  if (ht==NULL) return (NULL);
  for (i=0; i<table_len; i++) ht->entry[i]=NULL;
  ht->table_len=table_len;
  ht->n_obj=0;
  ht->lock=0;
  ht->next_get=0;
  ht->next_lid=0;

  return (ht);
}

int DisposeHashTable (HASH_TABLE *ht)
{
  free((void*)ht);

  return (0);
}

int HashFunc (HASH_TABLE *ht, int *key)
{
  int i,skey;

  skey=HASH_MULTIPLICATOR;
  for (i=1; i<=key[0]; i++)
    skey *= key[i];
  skey %= ht->table_len;
  return (skey);
}

int PushHashEntry (HASH_TABLE *ht, int *key, void *obj)
{
  int i,skey;

  if (ht->lock)
  {
    printf("ERROR in 'PopHashEntry': cannot push locked hashtable\n");
    return (1);
  }
  if (ht->n_obj>=ht->table_len)
  {
    printf("ERROR in 'PushHashEntry': overflow\n");
    return (1);
  }

  skey = HashFunc(ht,key);
  while (ht->entry[skey]!=NULL) skey=(skey+1)%ht->table_len;
  ht->entry[skey]=(ENTRY *)malloc(sizeof(ENTRY)+key[0]*sizeof(int));
  if (ht->entry[skey]==NULL)
  {
    printf("ERROR in 'PushHashEntry': cannot allocate entry\n");
    return (1);
  }
  ht->entry[skey]->obj=obj;
  for (i=0; i<=key[0]; i++)
    ht->entry[skey]->key[i]=key[i];
  ht->n_obj++;
  ht->entry[skey]->lid=ht->next_lid;
  ht->next_lid++;

  return (0);
}

void *PopHashEntry (HASH_TABLE *ht, int *key)
{
  int skey,i,rkey;
  void *obj;

  if (ht->lock)
  {
    printf("ERROR in 'PopHashEntry': cannot pop locked hashtable\n");
    return (NULL);
  }

  skey = HashFunc(ht,key);
  for (rkey=skey; rkey!=(skey-1+ht->table_len)%ht->table_len; rkey=(rkey+1)%ht->table_len)
  {
    if (ht->entry[rkey]==NULL || ht->entry[rkey]->key[0]!=key[0]) continue;
    for (i=1; i<=key[0]; i++)
      if (ht->entry[rkey]->key[i]!=key[i])
        break;
    if (i<=key[0]) continue;
    obj=ht->entry[rkey]->obj;
    free((void*)ht->entry[rkey]);
    ht->n_obj--;
    return(obj);
  }

  return (NULL);
}

int BeginHashGet (HASH_TABLE *ht)
{
  if (ht->lock)
  {
    printf("ERROR in 'BeginHashGet': hashtable already locked\n");
    return (1);
  }
  ht->lock=1;
  ht->next_get=0;

  return (0);
}

int HashGet (HASH_TABLE *ht, void **obj, int *lid)
{
  int i;

  if (ht->next_get>=ht->table_len)
  {
    *obj=NULL;
    return (0);
  }
  for (i=ht->next_get; i<ht->table_len; i++)
    if (ht->entry[i]!=NULL)
      break;
  ht->next_get=i+1;
  if (i==ht->table_len) *obj=NULL;
  else
  {
    *obj=ht->entry[i]->obj;
    *lid=ht->entry[i]->lid;
  }

  return (0);
}

int EndHashGet (HASH_TABLE *ht)
{
  if (!ht->lock)
  {
    printf("ERROR in 'EndHashGet': hashtable not locked\n");
    return (1);
  }
  ht->lock=0;
  ht->next_get=0;

  return (0);
}

int LocalIndexHash (HASH_TABLE *ht, int *key, int *lid)
{
  int skey,i,rkey;
  void *obj;

  *lid=-1;
  skey = HashFunc(ht,key);
  for (rkey=skey; rkey!=(skey-1+ht->table_len)%ht->table_len; rkey=(rkey+1)%ht->table_len)
  {
    if (ht->entry[rkey]==NULL || ht->entry[rkey]->key[0]!=key[0]) continue;
    for (i=1; i<=key[0]; i++)
      if (ht->entry[rkey]->key[i]!=key[i])
        break;
    if (i<=key[0]) continue;
    *lid=ht->entry[rkey]->lid;
    break;
  }

  return (0);
}



/****************************************************************************/
/******************************* ugpfm-mg ***********************************/
/****************************************************************************/

int MergeMultigrid (char *in, int rename)
{
  MGIO_MG_GENERAL mg_general,*mg_general_list,mg_general_dummy;
  MGIO_GE_GENERAL ge_general;
  MGIO_GE_ELEMENT ge_element[TAGS];
  MGIO_RR_GENERAL rr_general;
  MGIO_RR_RULE *rr_rules;
  MGIO_CG_GENERAL cg_general;
  MGIO_CG_POINT *cg_point;
  MGIO_CG_ELEMENT *cg_element;
  MGIO_BD_GENERAL bd_general;
  MGIO_PARINFO cg_pinfo;
  MGIO_REFINEMENT refinement;
  unsigned short *ProcList;
  BNDP **BndPList;
  char prefix[128],appdix[128],outname[128],tmp[128],tmp2[28],*p;
  int i,j,k,non,foid,tag,*vidlist;

  /*************************************************************************/
  /************************ read input file ********************************/
  /*************************************************************************/

  /* open proc 0  file */
  sprintf(tmp,"%s/mg.0000",in);
  if (Read_OpenMGFile(tmp))                                               {printf("ERROR in 'MergeMultigrid': cannot open proc 0 file\n");return (1);}
  if (Read_MG_General(&mg_general))                               {printf("ERROR in 'MergeMultigrid': cannot read mg_general 0 file\n");return (1);}
  if (strcmp(mg_general.version,MGIO_VERSION)!=0) {printf("ERROR in 'MergeMultigrid': version mismatch\n");return (1);}
  nparfiles=mg_general.nparfiles;
  if (nparfiles<=1)                                                               {printf("ERROR in 'MergeMultigrid': cannot merge %d parfile\n",nparfiles);return (1);}
  if (CloseMGFile())                                                              {printf("ERROR in 'MergeMultigrid': cannot close proc 0 file\n");return (1);}

  /* allocate dynamic lists */
  mg_general_list=(MGIO_MG_GENERAL*)malloc(nparfiles*sizeof(MGIO_MG_GENERAL));
  if (mg_general_list==NULL)                                              {printf("ERROR in 'MergeMultigrid': cannot allocate mg_general_list\n");return (1);}

  /* read all mg_generals */
  for (i=0; i<nparfiles; i++)
  {
    sprintf(tmp,"%s/mg.%04d",in,i);
    if (Read_OpenMGFile(tmp))                                       {printf("ERROR in 'MergeMultigrid': cannot open proc %d file\n",i);return (1);}
    if (Read_MG_General(mg_general_list+i))         {printf("ERROR in 'MergeMultigrid': cannot read mg_general %d file\n",i);return (1);}
    if (CloseMGFile())                                                      {printf("ERROR in 'MergeMultigrid': cannot close proc %d file\n",i);return (1);}
  }

  /* create mg_general for out-file */
  mg_general.mode=mg_general_list[0].mode;
  strcpy(mg_general.version,mg_general_list[0].version);
  mg_general.magic_cookie=mg_general_list[0].magic_cookie;
  strcpy(mg_general.ident,mg_general_list[0].ident);
  mg_general.nparfiles=1;
  mg_general.me=0;
  mg_general.nLevel=mg_general_list[0].nLevel;
  mg_general.nNode=mg_general.nPoint=mg_general.nElement=0;
  for (i=0; i<nparfiles; i++)
  {
    mg_general.nNode+=mg_general_list[i].nNode;
    mg_general.nPoint+=mg_general_list[i].nPoint;
    mg_general.nElement+=mg_general_list[i].nElement;
  }
  mg_general.dim=mg_general_list[0].dim;
  strcpy(mg_general.DomainName,mg_general_list[0].DomainName);
  strcpy(mg_general.MultiGridName,mg_general_list[0].MultiGridName);
  strcpy(mg_general.Formatname,mg_general_list[0].Formatname);
  mg_general.heapsize=nparfiles*mg_general_list[0].heapsize;
  mg_general.VectorTypes=mg_general_list[0].VectorTypes;

  /* create hashes */

  /* scan each input file */
  for (i=0; i<nparfiles; i++)
  {
    sprintf(tmp,"%s/mg.%04d",in,i);
    if (Read_OpenMGFile(tmp))                                                       {printf("ERROR in 'MergeMultigrid': cannot open proc %d file\n",i);return (1);}
    if (Read_MG_General(&mg_general_dummy))                         {printf("ERROR in 'MergeMultigrid': cannot read mg_general of proc %d file\n",i);return (1);}
    if (Read_GE_General(&ge_general))                       {printf("ERROR in 'MergeMultigrid': cannot read 'ge_general' of proc %d file\n",i);return (1);}
    if (Read_GE_Elements(TAGS,ge_element))                  {printf("ERROR in 'MergeMultigrid': cannot read 'ge_element' of proc %d file\n",i);return (1);}
    if (Read_RR_General(&rr_general))                       {printf("ERROR in 'MergeMultigrid': cannot read 'rr_general' of proc %d file\n",i);return (1);}
    rr_rules = (MGIO_RR_RULE *)malloc(rr_general.nRules*sizeof(MGIO_RR_RULE));
    if (rr_rules==NULL)                                                             {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'rr_rules' of proc %d file\n",i);return (1);}
    if (Read_RR_Rules(rr_general.nRules,rr_rules))          {printf("ERROR in 'MergeMultigrid': cannot read 'rr_rules' of proc %d file\n",i);return (1);}
    if (Read_CG_General(&cg_general))                                       {printf("ERROR in 'MergeMultigrid': cannot read 'cg_general' of proc %d file\n",i);return (1);}
    if (cg_general.nPoint>0)
    {
      cg_point=(MGIO_CG_POINT*)malloc(cg_general.nPoint*sizeof(MGIO_CG_POINT));
      if (cg_point==NULL)                                                                     {printf("ERROR in 'MergeMultigrid': cannot allocate 'cg_point' of proc %d file\n",i);return (1);}
      if (Read_CG_Points(cg_general.nPoint,cg_point))         {printf("ERROR in 'MergeMultigrid': cannot read 'cg_point' of proc %d file\n",i);return (1);}
    }
    if (Bio_Read_mint(1,&non))                              {printf("ERROR in 'MergeMultigrid': cannot read 'non' of proc %d file\n",i);return (1);}
    if (Bio_Read_mint(1,&foid))                             {printf("ERROR in 'MergeMultigrid': cannot read 'foid' of proc %d file\n",i);return (1);}
    vidlist = (int*)malloc(non*sizeof(int));
    if (Bio_Read_mint(non,vidlist))                         {printf("ERROR in 'MergeMultigrid': cannot read 'vidlist' of proc %d file\n",i);return (1);}
    if (cg_general.nElement>0)
    {
      cg_element=(MGIO_CG_ELEMENT*)malloc(cg_general.nElement*sizeof(MGIO_CG_ELEMENT));
      if (cg_element==NULL)                                                                   {printf("ERROR in 'MergeMultigrid': cannot allocate 'cg_element' of proc %d file\n",i);return (1);}
      if (Read_CG_Elements(cg_general.nElement,cg_element))   {printf("ERROR in 'MergeMultigrid': cannot read 'cg_element' of proc %d file\n",i);return (1);}
    }
    if (Bio_Jump (0))                                                                       {printf("ERROR in 'MergeMultigrid': cannot 'Bio_Jump (0)' in proc %d file\n",i);return (1);}
    if (Read_BD_General (&bd_general))                                      {printf("ERROR in 'MergeMultigrid': cannot read 'bd_general' in proc %d file\n",i);return (1);}
    if (bd_general.nBndP>0)
    {
      BndPList = (BNDP**)malloc(bd_general.nBndP*sizeof(BNDP*));
      if (BndPList==NULL)                                                     {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'BndPList' of proc %d file\n",i);return (1);}
      if (Read_PBndDesc (NULL,NULL,bd_general.nBndP,BndPList))
      {printf("ERROR in 'MergeMultigrid': cannot read 'BndPList' of proc %d file\n",i);return (1);}
    }
    if (cg_general.nElement>0)
    {
      ProcList = (unsigned short*)malloc(PROCLISTSIZE*sizeof(unsigned short));
      if (ProcList==NULL)                                                             {printf("ERROR in 'MergeMultigrid': cannot allocate 'ProcList' of proc %d file\n",i);return (1);}
      cg_pinfo.proclist = ProcList;
      for (j=0; j<cg_general.nElement; j++)
      {
        tag=ge_element[cg_element[j].ge].tag;
        if (Read_pinfo(tag,&cg_pinfo))                                  {printf("ERROR in 'MergeMultigrid': cannot read 'cg_pinfo' of proc %d file\n",i);return (1);}
      }
      for (j=0; j<cg_general.nElement; j++)
        for (k=0; k<cg_element[j].nref; k++)
        {
          if (Read_Refinement(&refinement,rr_rules))  {printf("ERROR in 'MergeMultigrid': cannot read 'refinement' of proc %d file\n",i);return (1);}
        }
    }

    /* free memory */
    if (cg_general.nElement>0) free((void*)ProcList);
    if (bd_general.nBndP>0) free((void*)BndPList);
    free((void*)vidlist);
    if (cg_general.nElement>0) free((void*)cg_element);
    if (cg_general.nPoint>0) free((void*)cg_point);


    if (CloseMGFile())                                                              {printf("ERROR in 'MergeMultigrid': cannot close proc 0 file\n");return (1);}

  }

  /*************************************************************************/
  /********************** write output file ********************************/
  /*************************************************************************/

  /* create outname */
  strcpy(tmp,in);p=strtok(tmp,".");
  if (p==NULL)                                                                    {printf("ERROR in 'MergeMultigrid': cannot create outfilename\n");return (1);}
  p+=strlen(p)+1;
  strcpy(outname,tmp);sprintf(tmp2,"_%d.",(int)nparfiles);strcat(outname,tmp2);strcat(outname,p);

  /* write part 0 to output file */
  if (Write_OpenMGFile(outname,rename))                   {printf("ERROR in 'MergeMultigrid': cannot open output file\n");return (1);}
  if (Write_MG_General(&mg_general))                              {printf("ERROR in 'MergeMultigrid': cannot 'mg_general' to output file\n");return (1);}
  if (Write_GE_General(&ge_general))                          {printf("ERROR in 'MergeMultigrid': cannot 'ge_general' to output file\n");return (1);}
  if (Write_GE_Elements(TAGS,ge_element))             {printf("ERROR in 'MergeMultigrid': cannot 'ge_element' to output file\n");return (1);}
  if (Write_RR_General(&rr_general))                          {printf("ERROR in 'MergeMultigrid': cannot 'rr_general' to output file\n");return (1);}
  if (Write_RR_Rules(rr_general.nRules,rr_rules)) {printf("ERROR in 'MergeMultigrid': cannot 'rr_rules' to output file\n");return (1);}
  if (CloseMGFile())                                                              {printf("ERROR in 'MergeMultigrid': cannot close output file\n");return (1);}

  return (0);
}




int main (int argc, char **argv)
{
  char in[128];

  if (argc<2)                                                                     {printf("filename required\n"); return (0);}
  strcpy(in,argv[1]);
  if (MergeMultigrid (in,0)) printf("some error\n");

  return (0);
}
