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

#define ModelP

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
#include "parallel.h"


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAX_KEY_LEN                                     10
#define HASH_MULTIPLICATOR                  43
#define HASH_STARTSIZE                          1000
#define HASH_RESIZEFRAC             0.5
#define ELEMPROCLISTSIZE                2000
#define PROCLISTSIZE                    (ELEMPROCLISTSIZE*MAX_SONS * MAX(5,(int)(2.0+log((double)nparfiles))))
#define VERBOSE
#define MERGE_DEBUG

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
  int ncol;                         /* nb. of collisions						*/
  int nmerge;                       /* nb. of merges							*/
  int mem;                                              /* memory used in bytes						*/
  ENTRY *entry[1];                      /* table used								*/
} HASH_TABLE;

typedef struct {
  int table_len;                                /* max number of objects					*/
  int n_obj;                        /* current number of objects                */
  int ncol;                         /* nb. of collisions                        */
  int nmerge;                       /* nb. of merges                            */
} HASH_STAT;

typedef struct {
  struct mgio_refinement_seq ref;
  int sonex;
} MERGE_REFINEMENT;

typedef int (*HashEntryProc)(void **object);
typedef int (*RefinementProcPtr)(MERGE_REFINEMENT *ref);

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
static int debug=0;

#ifdef MERGE_DEBUG
static HASH_TABLE *ht_nodes=NULL;
#endif
HASH_TABLE *ht_mem=NULL;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

void *PopHashEntry (HASH_TABLE *ht, int *key);
int PushHashEntry (HASH_TABLE **hth, int *key, void *obj);

/****************************************************************************/
/****************************** ugpfm-low ***********************************/
/****************************************************************************/

void *ht_malloc (size_t size, char *ident)
{
  int i,key[MAX_KEY_LEN+1];
  void *mem,*obj;

  mem=malloc(size);
  if (mem!=NULL)
  {
    assert(strlen(ident)<=MAX_KEY_LEN);
    for (i=0; i<strlen(ident); i++)
      key[i+1]=(int)(ident[i]);
    key[0]=i;
    obj=PopHashEntry(ht_mem,key);
    if (obj!=NULL) size+=(int)obj;
    if (PushHashEntry(&ht_mem,key,(void*)size)) return (NULL);
  }
  return (mem);
}

HASH_TABLE *CreateHashTable (HASH_TABLE **hth, int size)
{
  int i,j,table_len,max_factor;
  HASH_TABLE *ht;

  /* determine table_len */
  max_factor=1;
  for (i=2*size-1; i<2*size+200; i+=2)
  {
    for (j=3; j<sqrt(i); j+=2)
      if (i%j==0)
      {
        if (j>max_factor)
        {
          max_factor=j;
          table_len=i;
        }
        break;
      }
    if (j>=sqrt(i))
    {
      table_len=i;
      break;
    }
  }
  if (hth!=&ht_mem) ht = (HASH_TABLE *)ht_malloc(sizeof(HASH_TABLE)+(table_len-1)*sizeof(ENTRY*),"hashtable");
  else ht = (HASH_TABLE *)malloc(sizeof(HASH_TABLE)+(table_len-1)*sizeof(ENTRY*));
  if (ht==NULL) return (NULL);
  for (i=0; i<table_len; i++) ht->entry[i]=NULL;
  ht->table_len=table_len;
  ht->n_obj=0;
  ht->ncol=0;
  ht->lock=0;
  ht->next_get=0;
  ht->next_lid=0;
  ht->nmerge=0;
  ht->mem=sizeof(HASH_TABLE)+(table_len-1)*sizeof(ENTRY*);

  return (ht);
}

HASH_TABLE *ResizeHashTable (HASH_TABLE *ht, int add_obj)
{
  HASH_TABLE *nht;
  int i,ncol,nmerge;

  nht=CreateHashTable(&ht,ht->n_obj+add_obj);
  if (nht==NULL) return (NULL);
  for (i=ht->next_get; i<ht->table_len; i++)
  {
    if (ht->entry[i]==NULL) continue;
    if (PushHashEntry (&nht,ht->entry[i]->key,ht->entry[i]->obj))
      free(ht->entry[i]);
  }
  nht->ncol=ht->ncol;
  nht->nmerge=ht->nmerge;
  free((void*)ht);

  return (nht);
}

int HashFunc (HASH_TABLE *ht, int *key)
{
  int i,skey;

  skey=0;
  for (i=1; i<=key[0]; i++)
    skey += key[i];
  skey*=HASH_MULTIPLICATOR;
  skey %= ht->table_len;
  if (skey<0) skey=-skey;
  return (skey);
}

int key_eq (int* k1, int *k2)
{
  int i;
  if (k1[0]!=k2[0]) return (0);
  for (i=1; i<=k1[0]; i++)
    if (k1[i]!=k2[i])
      return (0);
  return (1);
}

int key_eq2 (int* k1, int *k2)
{
  int i,j,jfound,ifound;

  if (k1[0]!=k2[0]) return (0);
  ifound=0;
  for (i=1; i<=k1[0]; i++)
  {
    jfound=0;
    for (j=1; j<=k1[0]; j++)
      if (k2[j]==k1[i])
        jfound++;
    assert(jfound<=1);
    if (jfound==1) ifound++;
  }
  if (ifound==k1[0]) return (1);
  return (0);
}

int PushHashEntry (HASH_TABLE **hth, int *key, void *obj)
{
  HASH_TABLE *ht;
  int i,skey;

  if (*hth==NULL) *hth=CreateHashTable(hth,HASH_STARTSIZE);
  if (*hth==NULL)
  {
    printf("ERROR in 'PushHashEntry': cannot resize hashtable\n");
    return (1);
  }

  ht=*hth;
  if (ht->lock)
  {
    printf("ERROR in 'PushHashEntry': cannot push locked hashtable\n");
    return (1);
  }
  if (ht->n_obj>=HASH_RESIZEFRAC*ht->table_len) ht=*hth=ResizeHashTable(ht,ht->table_len);
  if (ht==NULL)
  {
    printf("ERROR in 'PushHashEntry': cannot resize hashtable\n");
    return (1);
  }

  skey = HashFunc(ht,key);
  while (ht->entry[skey]!=NULL && !key_eq(key,ht->entry[skey]->key)) {skey=(skey+1)%ht->table_len; ht->ncol++;}
  if (ht->entry[skey]!=NULL)
  {
    ht->entry[skey]->obj=obj;
    return (0);
  }
  if (hth!=&ht_mem) ht->entry[skey]=(ENTRY *)ht_malloc(sizeof(ENTRY)+key[0]*sizeof(int),"hashtable");
  else ht->entry[skey]=(ENTRY *)malloc(sizeof(ENTRY)+key[0]*sizeof(int));
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
  ht->mem+=sizeof(ENTRY)+key[0]*sizeof(int);

  return (0);
}

void *PopHashEntry (HASH_TABLE *ht, int *key)
{
  int skey,i,rkey;

  if (ht==NULL) return (NULL);
  skey = HashFunc(ht,key);
  for (rkey=skey; rkey!=(skey-1+ht->table_len)%ht->table_len; rkey=(rkey+1)%ht->table_len)
  {
    if (ht->entry[rkey]==NULL) return (NULL);
    if (!key_eq (key,ht->entry[rkey]->key)) continue;
    return(ht->entry[rkey]->obj);
  }
  return (NULL);
}

int BeginHashGet (HASH_TABLE *ht)
{
  if (ht==NULL) return (0);
  if (ht->lock)
  {
    printf("ERROR in 'BeginHashGet': hashtable already locked\n");
    return (1);
  }
  ht->lock=1;
  ht->next_get=0;

  return (0);
}

int HashGet (HASH_TABLE *ht, void **obj, int *key, int *lid, int *found)
{
  int i,j;

  if (ht==NULL) {*found=0; return (0);}
  if (ht->next_get>=ht->table_len)
  {
    *found=0;
    return (0);
  }
  for (i=ht->next_get; i<ht->table_len; i++)
    if (ht->entry[i]!=NULL)
      break;
  ht->next_get=i+1;
  if (i==ht->table_len) *found=0;
  else
  {
    *obj=ht->entry[i]->obj;
    *lid=ht->entry[i]->lid;
    for (j=0; j<=ht->entry[i]->key[0]; j++)
      key[j]=ht->entry[i]->key[j];
    *found=1;
  }

  return (0);
}

int EndHashGet (HASH_TABLE *ht)
{
  if (ht==NULL) return (0);
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
  if (ht==NULL) return (0);
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

int HashTableProc (HASH_TABLE *ht, HashEntryProc Proc, int *result)
{
  int i;

  *result=0;
  if (ht==NULL) return (0);
  for (i=0; i<ht->table_len; i++)
    if (ht->entry[i]!=NULL)
      *result += (*Proc)(&(ht->entry[i]->obj));

  return (0);
}

int HashTableStat (HASH_TABLE *ht, HASH_STAT *hst)
{
  if (ht!=NULL)
  {
    hst->table_len=ht->table_len;
    hst->n_obj=ht->n_obj;
    hst->ncol=ht->ncol;
    hst->nmerge=ht->nmerge;
  }
  else
  {
    hst->table_len=0;
    hst->n_obj=0;
    hst->ncol=0;
    hst->nmerge=0;
  }

  return (0);
}

int HashTablePrint (HASH_TABLE *ht, char *name)
{
  printf("######### %s #########\n");
  if (ht!=NULL)
  {
    printf("n_obj: %d\n",ht->n_obj);
    printf("n_col: %d\n",ht->ncol);
    printf("n_mrg: %d\n",ht->nmerge);
    printf("mem:   %d\n",ht->mem);
    printf("\n");
  }
  else
  {
    printf("n_obj: 0\n");
    printf("n_col: 0\n");
    printf("n_mrg: 0\n");
    printf("mem:   0\n");
    printf("\n");
  }

  return (0);
}

int HashTableInsertAtBegin (HASH_TABLE **ht, HASH_TABLE *insert)
{
  int i,key[MAX_KEY_LEN],lid,found,n_obj_ht,n_obj_in;
  void *object;

  if (*ht==NULL) n_obj_ht=0;
  else n_obj_ht=(*ht)->n_obj;
  if (insert==NULL) n_obj_in=0;
  else n_obj_in=insert->n_obj;
  if (BeginHashGet(insert)) return (1);
  while(1)
  {
    if (HashGet (insert,&object,key,&lid,&found)) return (1);
    if (!found) break;
    if (PushHashEntry (ht,key,object)) return (1);
  }
  if (EndHashGet(insert)) return (1);
  assert((*ht)->n_obj==n_obj_ht+n_obj_in);
  for (i=0; i<(*ht)->table_len; i++)
    if ((*ht)->entry[i]!=NULL)
      if ((*ht)->entry[i]->lid<n_obj_ht)
      {
        (*ht)->entry[i]->lid+=n_obj_in;
        assert((*ht)->entry[i]->lid>=n_obj_in);
        assert((*ht)->entry[i]->lid<(*ht)->table_len);
      }
      else
      {
        (*ht)->entry[i]->lid-=n_obj_ht;
        assert((*ht)->entry[i]->lid>=0);
        assert((*ht)->entry[i]->lid<n_obj_in);
      }

  return (0);
}

int ht_malloc_display (void)
{
  int total,i,lid,found,key[MAX_KEY_LEN+1];
  char name[MAX_KEY_LEN+1];
  void *obj;

  printf("%%%%%%%%%% ht_mem display %%%%%%%%%%\n");
  if (ht_mem==NULL) return (0);
  total=0;
  if (BeginHashGet(ht_mem)) return (1);
  while(1)
  {
    if (HashGet(ht_mem,&obj,key,&lid,&found)) return(1);
    if (!found) break;
    for (i=1; i<=key[0]; i++)
      name[i-1]=(char)key[i];
    name[key[0]]='\0';
    printf("%12s: %d bytes\n",name,(int)obj);
    total+=(int)obj;
  }
  if (EndHashGet(ht_mem)) return (1);
  printf("total: %d bytes\n",total);
  printf("\n");

  return (0);
}

/****************************************************************************/
/******************************* ugpfm-mg ***********************************/
/****************************************************************************/

int ReconstructNodeContext (MERGE_REFINEMENT *mref, MGIO_GE_ELEMENT *ge_element, MGIO_RR_RULE *rr_rules, int *nc)
{
  int i,j,sonex;
  struct mgio_refinement_seq *ref;

  ref=&(mref->ref);
  sonex=mref->sonex;
  for (i=0; i<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS; i++) nc[i]=-1;
  for (i=0; i<rr_rules[ref->refrule].nsons; i++)
    if (sonex & (1<<i))
      for (j=0; j<ge_element[rr_rules[ref->refrule].sons[i].tag].nCorner; j++)
      {
        assert(rr_rules[ref->refrule].sons[i].corners[j]>=0);
        assert(rr_rules[ref->refrule].sons[i].corners[j]<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS);
        nc[rr_rules[ref->refrule].sons[i].corners[j]]=1;
      }
  for (i=j=0; i<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS; i++) if (nc[i]!=-1) nc[i]=ref->newcornerid[j++];
  assert(ref->nnewcorners==j);

  return (0);
}

int MergeRefinement (HASH_TABLE **hth, int *key, MERGE_REFINEMENT *refinement, MGIO_GE_ELEMENT *ge_element, MGIO_RR_RULE *rr_rules)
{
  MERGE_REFINEMENT *popref;
  int i,ncold[MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS],ncnew[MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS];

  popref=(MERGE_REFINEMENT*)PopHashEntry(*hth,key);
  if (popref==NULL) return(PushHashEntry(hth,key,(void*)refinement));
  if (ReconstructNodeContext(popref,ge_element,rr_rules,ncold)) return (1);
  if (ReconstructNodeContext(refinement,ge_element,rr_rules,ncnew)) return (1);
  for (i=0; i<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS; i++)
  {
    if (ncnew[i]!=-1 && ncold[i]!=-1)       {assert(ncnew[i]==ncold[i]); continue;}
    if (ncnew[i]!=-1) continue;
    if (ncold[i]!=-1) ncnew[i]=ncold[i];
  }
  for (i=popref->ref.nnewcorners=0; i<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS; i++)
    if (ncnew[i]!=-1)
      popref->ref.newcornerid[popref->ref.nnewcorners++]=ncnew[i];
  popref->ref.sonref|=refinement->ref.sonref;
  popref->sonex|=refinement->sonex;

  return (0);
}

int InsertRefinement (HASH_TABLE **hth, int *key, MERGE_REFINEMENT *refinement, MGIO_GE_ELEMENT *ge_element, MGIO_RR_RULE *rr_rules, int *nref)
{
  int i,j,nc[MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS];
  static MERGE_REFINEMENT *next_ref;
  MERGE_REFINEMENT *ref;

  if (refinement!=NULL) next_ref=refinement;
  else next_ref++;
  ref=next_ref;
  if (MergeRefinement(hth,key,ref,ge_element,rr_rules))                   {printf("ERROR in 'InsertRefinement': cannot 'PopRefinement'\n");return (1);}
  if (ReconstructNodeContext(ref,ge_element,rr_rules,nc)) return (1);
  for (i=0; i<rr_rules[ref->ref.refrule].nsons; i++)
  {
    if (((ref->sonex&(1<<i))==0)||((ref->ref.sonref&(1<<i))==0)) continue;
    for (j=0; j<ge_element[rr_rules[ref->ref.refrule].sons[i].tag].nCorner; j++)
    {
      assert(nc[rr_rules[ref->ref.refrule].sons[i].corners[j]]!=-1);
      key[j+1]=nc[rr_rules[ref->ref.refrule].sons[i].corners[j]];
    }
    key[0]=j;
    if (InsertRefinement(hth,key,NULL,ge_element,rr_rules,nref))
    {printf("ERROR in 'InsertRefinement': cannot recursivly call myself\n");return (1);}
  }
  (*nref)++;

  return (0);
}


static int HT_ShiftEntry_Shift;

int HT_ShiftEntry (void **obj)
{
  int val;
  val=(int)(*obj);
  val+=HT_ShiftEntry_Shift;
  *obj=(void*)val;

  return (0);
}

void write_cge_gnu(FILE *stream, MGIO_CG_ELEMENT *elem, struct mgio_cg_point_seq *cg_point_out)
{
  int i,id;

  for (i=0; i<3; i++)
  {
    id=elem->cornerid[i];
    fprintf(stream,"%f %f;\n",cg_point_out[id].position[0],cg_point_out[id].position[1]);
    id=elem->cornerid[(i+1)%3];
    fprintf(stream,"%f %f;\n",cg_point_out[id].position[0],cg_point_out[id].position[1]);
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n");

  return;
}

#ifdef MERGE_DEBUG
static HASH_TABLE *ht_ref_crosscheck=NULL;
static int ref_crosscheck=0;
#endif
int ClimbRefinementTree(int *key, HASH_TABLE *ht_ref, RefinementProcPtr Proc, MGIO_GE_ELEMENT *ge_element, MGIO_RR_RULE *rr_rules)
{
  int i,j,sonkey[MGIO_MAX_CORNERS_OF_ELEM+1],nc[MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS];
  MERGE_REFINEMENT *ref;

  ref=(MERGE_REFINEMENT*)PopHashEntry(ht_ref,key);
  if (ref==NULL) return (0);

#ifdef MERGE_DEBUG
  if (ref_crosscheck)
    if (PushHashEntry(&ht_ref_crosscheck,key,NULL)) return (1);
#endif

  if (ReconstructNodeContext (ref,ge_element,rr_rules,nc)) return (1);
  if (Proc!=NULL) if ((*Proc)(ref)) return (1);
  for (i=0; i<rr_rules[ref->ref.refrule].nsons; i++)
  {
    if (((ref->sonex&(1<<i))==0)||((ref->ref.sonref&(1<<i))==0)) continue;
    for (j=0; j<ge_element[rr_rules[ref->ref.refrule].sons[i].tag].nCorner; j++)
    {
      sonkey[j+1]=nc[rr_rules[ref->ref.refrule].sons[i].corners[j]];
      assert(sonkey[j+1]!=-1);
    }
    sonkey[0]=j;
    if (ClimbRefinementTree(sonkey,ht_ref,Proc,ge_element,rr_rules)) return (1);
  }

  return (0);
}

static HASH_TABLE *WR_GOL;
static MGIO_RR_RULE *WR_rr_rules;
int WriteRefinement (MERGE_REFINEMENT *ref)
{
  int i,lid,key[2];

  for (i=0; i<ref->ref.nnewcorners; i++)
  {
    key[0]=1; key[1]=ref->ref.newcornerid[i];
    if (LocalIndexHash(WR_GOL,key,&lid)) return (1);
    assert(lid!=-1);
    ref->ref.newcornerid[i]=lid;
  }
  if (Write_Refinement((MGIO_REFINEMENT*)ref,WR_rr_rules)) return (1);

  return (0);
}

static int CR_nref;
int CountRefinements (MERGE_REFINEMENT *ref)
{
  CR_nref++;
  return (0);
}

int MergeMultigrid (char *in, int rename)
{
  HASH_TABLE *ht_cgv,*ht_ref,*ht_bnp, *ht_cgvlid;
  HASH_TABLE *ht_bn_l0, *ht_in_l0, *ht_gol,*ht_cge;
  HASH_STAT hst;
  MGIO_MG_GENERAL mg_general,*mg_general_list,mg_general_dummy;
  MGIO_GE_GENERAL ge_general;
  MGIO_GE_ELEMENT ge_element[TAGS];
  MGIO_RR_GENERAL rr_general;
  MGIO_RR_RULE *rr_rules;
  MGIO_CG_GENERAL *cg_general,cg_general_out;
  MGIO_CG_POINT **cg_point;
  struct mgio_cg_point_seq *cg_point_out;
  MGIO_CG_ELEMENT **cg_element,*o_element,*elem;
  MGIO_BD_GENERAL *bd_general,bd_general_out;
  MGIO_PARINFO cg_pinfo;
  MERGE_REFINEMENT ***refinement;
  MGIO_REFINEMENT loc_ref,*ref;
  unsigned short *ProcList;
  BNDP **BndPList;
  char prefix[128],appdix[128],outname[128],tmp[128],tmp2[28],*p;
  int i,j,k,l,s,t,non,foid,tag,*vidlist,key[MGIO_MAX_CORNERS_OF_ELEM+1],*ncge,n_ref_tot,nref_read,level,*in_lid2gid;
  int nc[MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS],*o_element_im,out_blid_offset,out_ilid_offset,nid_l0_max,vid_l0_max,vid_bl0_max,lid,error;
  int found,n_bn_l0,n_in_l0,id,gecid[MGIO_MAX_CORNERS_OF_ELEM],crosscheck_nref_tot;
  void *object;
  FILE *stream;

  stream=fopen("foo","w");

  /*************************************************************************/
  /************************ read input file ********************************/
  /*************************************************************************/

  /* open proc 0  file */
  sprintf(tmp,"%s/mg.0000",in);
  if (Read_OpenMGFile(tmp))                                               {printf("ERROR in 'MergeMultigrid': cannot open proc 0 file\n");return (1);}
  if (Read_MG_General(&mg_general))                               {printf("ERROR in 'MergeMultigrid': cannot read mg_general 0 file\n");return (1);}
  if (strcmp(mg_general.version,MGIO_VERSION)!=0) {printf("ERROR in 'MergeMultigrid': version mismatch\n");return (1);}
  nparfiles=mg_general.nparfiles;
  if (nparfiles<=0)                                                               {printf("ERROR in 'MergeMultigrid': cannot merge %d parfile\n",nparfiles);return (1);}
  if (CloseMGFile())                                                              {printf("ERROR in 'MergeMultigrid': cannot close proc 0 file\n");return (1);}

  /* allocate dynamic lists */
  mg_general_list=(MGIO_MG_GENERAL*)ht_malloc(nparfiles*sizeof(MGIO_MG_GENERAL),"const");
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

  /* prepare hashes */
  ht_cgv = NULL;
  ht_ref = NULL;
  ht_bnp = NULL;
  ht_bn_l0 = NULL;
  ht_in_l0 = NULL;
  ht_gol = NULL;

  /* scan each input file */
  cg_general=(MGIO_CG_GENERAL*)ht_malloc(nparfiles*sizeof(MGIO_CG_GENERAL),"const");
  if (cg_general==NULL)                                                                   {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'cg_general' \n");return (1);}
  cg_point=(MGIO_CG_POINT**)ht_malloc(nparfiles*sizeof(MGIO_CG_POINT*),"const");
  if (cg_point==NULL)                                                                             {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'cg_point' \n");return (1);}
  cg_element=(MGIO_CG_ELEMENT**)ht_malloc(nparfiles*sizeof(MGIO_CG_ELEMENT*),"const");
  if (cg_element==NULL)                                                                   {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'cg_element' \n");return (1);}
  bd_general=(MGIO_BD_GENERAL*)ht_malloc(nparfiles*sizeof(MGIO_BD_GENERAL),"const");
  if (bd_general==NULL)                                                                   {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'bd_general' \n");return (1);}
  refinement=(MERGE_REFINEMENT***)ht_malloc(nparfiles*sizeof(MERGE_REFINEMENT**),"const");
  if (refinement==NULL)                                                                   {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'refinement' \n");return (1);}
  ncge=(int*)ht_malloc(nparfiles*sizeof(int),"const");
  if (ncge==NULL)                                                                                 {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'ncge' \n");return (1);}
  printf("merging '%s': ",in);
  fflush(stdout);
  for (i=0; i<nparfiles; i++)
  {
    sprintf(tmp,"%s/mg.%04d",in,i);
    if (i<nparfiles-1) printf("[%d]",i);
    else printf("[%d]\n",i);
    fflush(stdout);
    if (Read_OpenMGFile(tmp))                                                       {printf("ERROR in 'MergeMultigrid': cannot open proc %d file\n",i);return (1);}
    if (Read_MG_General(&mg_general_dummy))                         {printf("ERROR in 'MergeMultigrid': cannot read mg_general of proc %d file\n",i);return (1);}
    in_lid2gid=(int*)ht_malloc(mg_general_dummy.nNode*sizeof(int),"procloc");
    if (in_lid2gid==NULL)                                                                   {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'in_lid2gid' of proc %d file\n",i);return (1);}
    for (j=0; j<mg_general_dummy.nNode; j++) in_lid2gid[j]=-1;
    nid_l0_max=vid_l0_max=vid_bl0_max=-1;
    if (Read_GE_General(&ge_general))                       {printf("ERROR in 'MergeMultigrid': cannot read 'ge_general' of proc %d file\n",i);return (1);}
    if (Read_GE_Elements(TAGS,ge_element))                  {printf("ERROR in 'MergeMultigrid': cannot read 'ge_element' of proc %d file\n",i);return (1);}
    if (Read_RR_General(&rr_general))                       {printf("ERROR in 'MergeMultigrid': cannot read 'rr_general' of proc %d file\n",i);return (1);}
    rr_rules = (MGIO_RR_RULE *)ht_malloc(rr_general.nRules*sizeof(MGIO_RR_RULE),"const");
    if (rr_rules==NULL)                                                             {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'rr_rules' of proc %d file\n",i);return (1);}
    if (Read_RR_Rules(rr_general.nRules,rr_rules))          {printf("ERROR in 'MergeMultigrid': cannot read 'rr_rules' of proc %d file\n",i);return (1);}
    if (Read_CG_General(cg_general+i))                                      {printf("ERROR in 'MergeMultigrid': cannot read 'cg_general' of proc %d file\n",i);return (1);}
    if (cg_general[i].nPoint>0)
    {
      cg_point[i]=(MGIO_CG_POINT*)ht_malloc(cg_general[i].nPoint*sizeof(MGIO_CG_POINT),"procglob");
      if (cg_point[i]==NULL)                                                                  {printf("ERROR in 'MergeMultigrid': cannot allocate 'cg_point' of proc %d file\n",i);return (1);}
      if (Read_CG_Points(cg_general[i].nPoint,cg_point[i]))   {printf("ERROR in 'MergeMultigrid': cannot read 'cg_point' of proc %d file\n",i);return (1);}
    }
    if (Bio_Read_mint(1,&non))                              {printf("ERROR in 'MergeMultigrid': cannot read 'non' of proc %d file\n",i);return (1);}
    if (Bio_Read_mint(1,&foid))                             {printf("ERROR in 'MergeMultigrid': cannot read 'foid' of proc %d file\n",i);return (1);}
    vidlist = (int*)ht_malloc(non*sizeof(int),"procloc");
    if (Bio_Read_mint(non,vidlist))                         {printf("ERROR in 'MergeMultigrid': cannot read 'vidlist' of proc %d file\n",i);return (1);}
    if (cg_general[i].nElement>0)
    {
      o_element=(MGIO_CG_ELEMENT*)ht_malloc(cg_general[i].nElement*sizeof(MGIO_CG_ELEMENT),"procloc");
      o_element_im=(int*)ht_malloc(cg_general[i].nElement*sizeof(int),"procloc");
      if (o_element==NULL)                                                    {printf("ERROR in 'MergeMultigrid': cannot allocate 'o_element' of proc %d file\n",i);return (1);}
      if (Read_CG_Elements(cg_general[i].nElement,o_element)) {printf("ERROR in 'MergeMultigrid': cannot read 'o_element' of proc %d file\n",i);return (1);}
      for (ncge[i]=0; ncge[i]<cg_general[i].nElement; ncge[i]++) if (o_element[ncge[i]].level>0) break;
      cg_element[i]=NULL;
      /* we use: 1.: level for global id */
    }
    if (Bio_Jump (0))                                                                       {printf("ERROR in 'MergeMultigrid': cannot 'Bio_Jump (0)' in proc %d file\n",i);return (1);}
    if (Read_BD_General (bd_general+i))                                     {printf("ERROR in 'MergeMultigrid': cannot read 'bd_general' in proc %d file\n",i);return (1);}
    if (bd_general[i].nBndP>0)
    {
      BndPList = (BNDP**)ht_malloc(bd_general[i].nBndP*sizeof(BNDP*),"const");
      if (BndPList==NULL)                                                     {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'BndPList' of proc %d file\n",i);return (1);}
      if (Read_PBndDesc (NULL,NULL,bd_general[i].nBndP,BndPList))
      {printf("ERROR in 'MergeMultigrid': cannot read 'BndPList' of proc %d file\n",i);return (1);}
    }
    if (cg_general[i].nElement>0)
    {
      ProcList = (unsigned short*)ht_malloc(PROCLISTSIZE*sizeof(unsigned short),"const");
      if (ProcList==NULL)                                                             {printf("ERROR in 'MergeMultigrid': cannot allocate 'ProcList' of proc %d file\n",i);return (1);}
      cg_pinfo.proclist = ProcList;
      for (j=0; j<cg_general[i].nElement; j++)
      {
        tag=ge_element[o_element[j].ge].tag;
        if (Read_pinfo(tag,&cg_pinfo))                                  {printf("ERROR in 'MergeMultigrid': cannot read 'cg_pinfo' of proc %d file\n",i);return (1);}
        if (EGHOSTPRIO(cg_pinfo.prio_elem)) o_element_im[j]=0;
        else o_element_im[j]=1;
        level=o_element[j].level; o_element[j].level=cg_pinfo.e_ident;
        for (k=0; k<ge_element[o_element[j].ge].nCorner; k++)
          if (cg_pinfo.prio_node[k]==PrioMaster)
          {
            assert(o_element[j].cornerid[k]<mg_general_dummy.nNode);
            assert(o_element[j].cornerid[k]>=foid && o_element[j].cornerid[k]<foid+non);
            l=o_element[j].cornerid[k];
            assert(in_lid2gid[l]==-1 || in_lid2gid[l]==cg_pinfo.n_ident[k]);
            in_lid2gid[l]=cg_pinfo.n_ident[k];
#ifdef MERGE_DEBUG
            key[0]=1; key[1]=cg_pinfo.n_ident[k];
            if (PushHashEntry(&ht_nodes,key,NULL))                          {printf("ERROR in 'MergeMultigrid': cannot insert vertex in 'ht_nodes' of proc %d file\n",i);return (1);}
#endif
          }
        if (level==0 && o_element_im[j])
          for (k=0; k<ge_element[o_element[j].ge].nCorner; k++)
            if (MASTERPRIO(cg_pinfo.prio_node[k]))
            {
              assert(o_element[j].cornerid[k]>=foid);
              assert(o_element[j].cornerid[k]<non+foid);
              id=vidlist[o_element[j].cornerid[k]-foid];
              key[0]=1; key[1]=cg_pinfo.n_ident[k];
              if (PushHashEntry(&ht_cgv,key,(void*)(cg_point[i]+id)))
              {printf("ERROR in 'MergeMultigrid': cannot insert vertex in 'ht_cgv' of proc %d file\n",i);return (1);}
              if (nid_l0_max<o_element[j].cornerid[k]) nid_l0_max=o_element[j].cornerid[k];
              if (vid_l0_max<id) vid_l0_max=id;
              if (id<bd_general[i].nBndP)
              {
                if (PushHashEntry(&ht_bnp,key,(void*)(BndPList[id])))
                {printf("ERROR in 'MergeMultigrid': cannot insert vertex in 'ht_bnp' of proc %d file\n",i);return (1);}
                if (vid_bl0_max<id) vid_bl0_max=id;
              }
            }
        for (k=0; k<ge_element[o_element[j].ge].nCorner; k++)
          o_element[j].cornerid[k]=cg_pinfo.n_ident[k];
      }
      refinement[i]=(MERGE_REFINEMENT**)ht_malloc(cg_general[i].nElement*sizeof(MERGE_REFINEMENT*),"refinement");
      if (refinement[i]==NULL)                                                        {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'refinement' \n");return (1);}
      for (j=0; j<cg_general[i].nElement; j++)
      {
        if (o_element[j].nref==0) continue;
        refinement[i][j]=(MERGE_REFINEMENT*)ht_malloc(o_element[j].nref*sizeof(MERGE_REFINEMENT),"refinement");
        ref=&loc_ref;
        for (l=0; l<MGIO_MAX_SONS_OF_ELEM; l++)
          ref->pinfo[l].proclist=ProcList;
        for (k=0; k<o_element[j].nref; k++)
        {
          if (Read_Refinement(ref,rr_rules))  {printf("ERROR in 'MergeMultigrid': cannot read 'refinement' of proc %d file\n",i);return (1);}
          if (ref->sonex==0) assert(ref->nnewcorners==0);
          for (l=0; l<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS; l++) nc[l]=-1;
          for (l=0; l<rr_rules[ref->refrule].nsons; l++)
            if (ref->sonex & (1<<l))
              for (s=0; s<ge_element[rr_rules[ref->refrule].sons[l].tag].nCorner; s++)
              {
                assert(rr_rules[ref->refrule].sons[l].corners[s]>=0);
                assert(rr_rules[ref->refrule].sons[l].corners[s]<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS);
                nc[rr_rules[ref->refrule].sons[l].corners[s]]=1;
              }
          for (l=s=0; l<MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS; l++) if (nc[l]!=-1) nc[l]=s++;
          assert(ref->nnewcorners==s);
          for (l=0; l<rr_rules[ref->refrule].nsons; l++)
            if (ref->sonex & (1<<l))
              for (s=0; s<ge_element[rr_rules[ref->refrule].sons[l].tag].nCorner; s++)
              {
                if (ref->pinfo[l].prio_node[s]!=PrioMaster) continue;
                t=nc[rr_rules[ref->refrule].sons[l].corners[s]];
                if (t==-1) continue;
                assert(ref->newcornerid[t]>=foid);
                assert(ref->newcornerid[t]<mg_general_dummy.nNode);
                t=ref->newcornerid[t];
                in_lid2gid[t]=ref->pinfo[l].n_ident[s];
                assert(in_lid2gid[t]==-1 || in_lid2gid[t]==ref->pinfo[l].n_ident[s]);
#ifdef MERGE_DEBUG
                key[0]=1; key[1]=ref->pinfo[l].n_ident[s];
                if (PushHashEntry(&ht_nodes,key,NULL))                          {printf("ERROR in 'MergeMultigrid': cannot insert vertex in 'ht_nodes' of proc %d file\n",i);return (1);}
#endif
              }
          for (l=0; l<rr_rules[ref->refrule].nsons; l++)
            if (ref->sonex & (1<<l))
              for (s=0; s<ge_element[rr_rules[ref->refrule].sons[l].tag].nCorner; s++)
              {
                t=nc[rr_rules[ref->refrule].sons[l].corners[s]];
                if (t==-1) continue;
                assert(ref->newcornerid[t]>=0);
                assert(ref->newcornerid[t]<mg_general_dummy.nNode);
                ref->newcornerid[t]=ref->pinfo[l].n_ident[s];
                nc[rr_rules[ref->refrule].sons[l].corners[s]]=-1;
              }
          memcpy((void*)(refinement[i][j]+k),(const void*)ref,sizeof(struct mgio_refinement_seq));
          refinement[i][j][k].sonex=ref->sonex;
        }
      }
      for (j=0; j<cg_general[i].nElement; j++)
      {
        if (o_element[j].nref==0) continue;
        for (k=0; k<ge_element[o_element[j].ge].nCorner; k++)
          key[k+1]=o_element[j].cornerid[k];
        key[0]=k; nref_read=0;
        if (InsertRefinement(&ht_ref,key,refinement[i][j],ge_element,rr_rules,&nref_read))
        {printf("ERROR in 'MergeMultigrid': cannot insert local refinement tree of proc %d file\n",i);return (1);}
        if (nref_read!=o_element[j].nref)
        {printf("ERROR in 'MergeMultigrid': nref mismatch in  refinement tree of proc %d file\n",i);return (1);}
      }
      cg_element[i]=(MGIO_CG_ELEMENT*)ht_malloc(ncge[i]*sizeof(MGIO_CG_ELEMENT),"procglob");
      if (cg_element[i]==NULL)                                                {printf("ERROR in 'MergeMultigrid': cannot allocate 'cg_element[i]' of proc %d file\n",i);return (1);}
      memcpy((void*)cg_element[i],(const void*)o_element,ncge[i]*sizeof(MGIO_CG_ELEMENT));
      for (j=0; j<ncge[i]; j++)                   /* globalize neighbors */
      {
        if (o_element_im[j]==0) {cg_element[i][j].ge=-1; continue;}
        for (k=0; k<ge_element[cg_element[i][j].ge].nSide; k++)
          if (cg_element[i][j].nbid[k]!=-1)
            cg_element[i][j].nbid[k]=cg_element[i][cg_element[i][j].nbid[k]].level;
      }
    }

    /* insert in hash tables */
    key[0]=1;
    for (j=foid; j<=nid_l0_max; j++)
    {
      if (in_lid2gid[j]==-1) continue;
      key[1]=in_lid2gid[j];
      k=vidlist[j-foid];
      if (k<=vid_bl0_max)
      {
        if (PushHashEntry(&(ht_bn_l0),key,NULL))  {printf("ERROR in 'MergeMultigrid': cannot insert bnd-node in 'ht_bn_l0' of proc %d file\n",i);return (1);}
      }
      else
      {
        assert(k<=vid_l0_max);
        if (PushHashEntry(&(ht_in_l0),key,NULL))  {printf("ERROR in 'MergeMultigrid': cannot insert bnd-node in 'ht_in_l0' of proc %d file\n",i);return (1);}
      }
    }
    for (; j<mg_general_dummy.nNode; j++)
    {
      if (in_lid2gid[j]==-1) continue;
      key[1]=in_lid2gid[j];
      if (PushHashEntry(&(ht_gol),key,NULL))  {printf("ERROR in 'MergeMultigrid': cannot insert bnd-node in 'ht_hn' of proc %d file\n",i);return (1);}
    }

    /* free memory */
    if (cg_general[i].nElement>0) free((void*)o_element);
    if (cg_general[i].nElement>0) free((void*)o_element_im);
    free((void*)vidlist);
    free((void*)in_lid2gid);


    if (CloseMGFile())                                                              {printf("ERROR in 'MergeMultigrid': cannot close proc 0 file\n");return (1);}

  }

  HashTableStat(ht_bn_l0,&hst); n_bn_l0=hst.n_obj;
  HashTableStat(ht_in_l0,&hst); n_in_l0=hst.n_obj;
  if (HashTableInsertAtBegin (&ht_gol,ht_in_l0)) {printf("ERROR in 'MergeMultigrid': cannot insert hash table\n");return (1);}
  if (HashTableInsertAtBegin (&ht_gol,ht_bn_l0)) {printf("ERROR in 'MergeMultigrid': cannot insert hash table\n");return (1);}

#ifdef VERBOSE
  HashTablePrint(ht_cgv,"vtx");
  HashTablePrint(ht_bnp,"bnp");
  HashTablePrint(ht_ref,"ref");
#ifdef MERGE_DEBUG
  HashTableStat(ht_ref,&hst);
  crosscheck_nref_tot=hst.n_obj;
#endif
  HashTablePrint(ht_gol,"gol");
#ifdef MERGE_DEBUG
  HashTablePrint(ht_nodes,"nds");
#endif
#endif

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

  /* write part 1 to output file */
  HashTableStat(ht_cgv,&hst);
  cg_general_out.nPoint=hst.n_obj;
  HashTableStat(ht_bnp,&hst);
  cg_general_out.nBndPoint=hst.n_obj;
  cg_general_out.nInnerPoint=cg_general_out.nPoint-cg_general_out.nBndPoint;
  cg_general_out.nElement=cg_general_out.nBndElement=0;
  for (i=0; i<nparfiles; i++)
    for (j=0; j<ncge[i]; j++)
    {
      if (cg_element[i][j].ge < 0) continue;
      cg_general_out.nElement++;
      for (k=0; k<ge_element[cg_element[i][j].ge].nSide; k++)
        if (cg_element[i][j].se_on_bnd & (1<<k))
        {
          cg_general_out.nBndElement++;
          break;
        }
    }
  cg_general_out.nInnerElement=cg_general_out.nElement-cg_general_out.nBndElement;
#ifdef VERBOSE
  printf("coarse grid info: nPoints:     %d\n"
         "                  nInnerPoints:%d\n"
         "                  nBnbPoints:  %d\n"
         "                  nElem:       %d\n"
         "                  nInnerElem:  %d\n"
         "                  nBndElem:    %d\n"
         ,cg_general_out.nPoint,cg_general_out.nInnerPoint,cg_general_out.nBndPoint,cg_general_out.nElement,cg_general_out.nInnerElement,cg_general_out.nBndElement);
#endif
  if (Write_CG_General(&cg_general_out))                                              {printf("ERROR in 'MergeMultigrid': cannot write 'cg_general' to output file\n");return (1);}
  cg_point_out=(struct mgio_cg_point_seq*)ht_malloc(cg_general_out.nPoint*sizeof(struct mgio_cg_point_seq),"out");
  if (cg_point_out==NULL)                                                                         {printf("ERROR in 'MergeMultigrid': cannot allocate array for 'cg_point for output file\n");return (1);}
  if (BeginHashGet(ht_cgv))                                                                       {printf("ERROR in 'MergeMultigrid': cannot 'BeginHashGet' of 'ht_cgv' for output file\n");return (1);}
  while (1)
  {
    if (HashGet(ht_cgv,&object,key,&lid,&found))                    {printf("ERROR in 'MergeMultigrid': cannot 'HashGet' of 'ht_cgv' for output file\n");return (1);}
    if (!found) break;
    if (LocalIndexHash(ht_gol,key,&lid))                    {printf("ERROR in 'MergeMultigrid': cannot 'LocalIndexHash'of 'ht_gol' for output file\n");return (1);}
    assert(lid>=0);
    assert(lid<cg_general_out.nPoint);
    for (i=0; i<mg_general.dim; i++)
      cg_point_out[lid].position[i]=((double*)object)[i];
  }
  if (EndHashGet(ht_cgv))                                                                         {printf("ERROR in 'MergeMultigrid': cannot 'EndHashGet' of 'ht_cgv' for output file\n");return (1);}
  if (Write_CG_Points(cg_general_out.nPoint,(MGIO_CG_POINT*)cg_point_out))
  {printf("ERROR in 'MergeMultigrid': cannot write 'cg_point_out' to output file\n");return (1);}
  ht_cge=NULL;
  for (i=0; i<nparfiles; i++)
    for (j=0; j<ncge[i]; j++)
      if (cg_element[i][j].ge!=-1)
      {
        key[0]=1; key[1]=cg_element[i][j].level;
        if (PushHashEntry(&ht_cge,key,(void*)(&(cg_element[i][j]))))
        {printf("ERROR in 'MergeMultigrid': cannot push 'cg_element' to 'ht_cge'\n");return (1);}
      }
#ifdef MERGE_DEBUG
  HashTablePrint(ht_cge,"cge");
#endif
  for (i=0; i<nparfiles; i++)
    for (j=0; j<ncge[i]; j++)
      if (cg_element[i][j].ge!=-1)
      {
        elem=&(cg_element[i][j]);
        for (k=0; k<ge_element[elem->ge].nSide; k++)
        {
          if (elem->nbid[k]==-1) continue;
          key[0]=1; key[1]=elem->nbid[k];
          if (LocalIndexHash(ht_cge,key,&lid))        {printf("ERROR in 'MergeMultigrid': cannot 'LocalIndexHash'of 'ht_cge' for output file\n");return (1);}
          assert(lid>=0);
          assert(lid<cg_general_out.nElement);
          elem->nbid[k]=lid;
        }
        for (k=0; k<ge_element[elem->ge].nCorner; k++)
          key[k+1]=elem->cornerid[k];
        key[0]=k;
        CR_nref=0;
        if (ClimbRefinementTree(key,ht_ref,CountRefinements,ge_element,rr_rules))
        {printf("ERROR in 'MergeMultigrid': cannot 'ClimbRefinementTree' for output file\n");return (1);}
        elem->nref=CR_nref;
        for (k=0; k<ge_element[elem->ge].nCorner; k++)
        {
          key[0]=1; key[1]=gecid[k]=elem->cornerid[k];
          if (LocalIndexHash(ht_gol,key,&lid))        {printf("ERROR in 'MergeMultigrid': cannot 'LocalIndexHash'of 'ht_gol' for output file\n");return (1);}
          assert(lid!=-1);
          elem->cornerid[k]=lid;
        }
        if (Write_CG_Elements(1,elem))                                  {printf("ERROR in 'MergeMultigrid': cannot 'Write_CG_Elements' for output file\n");return (1);}
        for (k=0; k<ge_element[elem->ge].nCorner; k++)
          elem->cornerid[k]=gecid[k];
      }
  if (Bio_Jump_From ())                                                                           {printf("ERROR in 'MergeMultigrid': cannot 'Bio_Jump_From' for output file\n");return (1);}
  bd_general_out.nBndP=cg_general_out.nBndPoint;
  if (Write_BD_General(&bd_general_out))                                          {printf("ERROR in 'MergeMultigrid': cannot 'Write_BD_General' for output file\n");return (1);}
  BndPList=(BNDP**)ht_malloc(bd_general_out.nBndP*sizeof(BNDP*),"out");
  if (BndPList==NULL)                                                                                     {printf("ERROR in 'MergeMultigrid': cannot allocate 'BndPList' for output file\n");return (1);}
#ifdef MERGE_DEBUG
  for (i=0; i<bd_general_out.nBndP; i++)
    BndPList[i]=NULL;
#endif
  if (BeginHashGet(ht_bnp))                                                                       {printf("ERROR in 'MergeMultigrid': cannot 'BeginHashGet' of 'ht_bnp' for output file\n");return (1);}
  while (1)
  {
    if (HashGet(ht_bnp,&object,key,&lid,&found))            {printf("ERROR in 'MergeMultigrid': cannot 'HashGet' of 'ht_cge' for output file\n");return (1);}
    if (!found) break;
    if (LocalIndexHash(ht_gol,key,&lid))                    {printf("ERROR in 'MergeMultigrid': cannot 'LocalIndexHash'of 'ht_gol' for output file\n");return (1);}
    assert(lid>=0);
    assert(lid<bd_general_out.nBndP);
#ifdef MERGE_DEBUG
    assert(BndPList[lid]==NULL);
#endif
    BndPList[lid]=(BNDP*)object;
  }
  if (EndHashGet(ht_bnp))                                     {printf("ERROR in 'MergeMultigrid': cannot 'EndHashGet' of 'ht_bnp' for output file\n");return (1);}
#ifdef MERGE_DEBUG
  for (i=0; i<bd_general_out.nBndP; i++)
    assert(BndPList[i]!=NULL);
#endif
  if (Write_PBndDesc(-bd_general_out.nBndP,BndPList))                     {printf("ERROR in 'MergeMultigrid': cannot 'Write_PBndDesc' of 'BndPList' for output file\n");return (1);}
  if (Bio_Jump_To ())                                                                             {printf("ERROR in 'MergeMultigrid': cannot 'Bio_Jump_To' for output file\n");return (1);}
  WR_GOL=ht_gol;
  WR_rr_rules=rr_rules;
#ifdef MERGE_DEBUG
  ref_crosscheck=1;
#endif
  for (i=0; i<nparfiles; i++)
    for (j=0; j<ncge[i]; j++)
      if (cg_element[i][j].ge!=-1)
      {
        elem=&(cg_element[i][j]);
        for (k=0; k<ge_element[elem->ge].nCorner; k++)
          key[k+1]=elem->cornerid[k];
        key[0]=k;
        if (ClimbRefinementTree(key,ht_ref,WriteRefinement,ge_element,rr_rules))
        {printf("ERROR in 'MergeMultigrid': cannot 'ClimbRefinementTree' for output file\n");return (1);}
      }

  /* cross-check */
  assert(n_bn_l0==cg_general_out.nBndPoint);
  assert(n_in_l0==cg_general_out.nInnerPoint);
#ifdef MERGE_DEBUG
  HashTablePrint(ht_ref_crosscheck,"refcrosscheck");
  HashTableStat(ht_ref_crosscheck,&hst);
  assert(hst.nmerge==0);
  assert(hst.n_obj==crosscheck_nref_tot);
#endif

  /* we are done */
  if (CloseMGFile())                                                              {printf("ERROR in 'MergeMultigrid': cannot close output file\n");return (1);}

#ifdef VERBOSE
  printf("\n");
#endif

  fclose(stream);

  /* display memory requirements */
  if (ht_malloc_display()) return (1);

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
