// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  er.c															*/
/*																			*/
/* Purpose:   extract rules realized in a multigrid                                             */
/*																			*/
/* Author:	  Henrik Rentz-Reichert			                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   06.12.97 begin, ug version 3.7								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#undef __DEBUG_ER__
#undef __PATCH_UG_RULES__

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>

/* low */
#include "general.h"

/* gm */
#include "rm.h"
#include "mgio.h"
#include "ugm.h"
#include "GenerateRules.h"

/* own header */
#include "er.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* TODO (HRR 971207): increase REFINE_LEN by shifting REFINE-ce to flag cw */
#define MAX_HRID                                        (1<<REFINE_LEN)
#define MAX_IFID                                        MAX_HRID

#ifdef __PATCH_UG_RULES__
/* CAUTION: save/load ok but REFINE(e) will not be consistent with rm anymore */
        #define UGMAXRULE(tag)                  ((MaxRules[tag]>0) ? 1 : 0)       /* include NO_REFINEMENT from ug */
        #define HAS_NO_RULE(e)                  (REFINE(e)>=UGMAXRULE(TAG(e)))
#else
        #define UGMAXRULE(tag)                  (MaxRules[tag])
        #define HAS_NO_RULE(e)                  (REFINE(e)==COPY && REFINECLASS(e)==GREEN_CLASS)
#endif

#define BEYOND_UG_RULES(e)                      (REFINE(e)>=UGMAXRULE(TAG(e)))

#define NOTDONE         -1

/* for computing son paths */
enum NB_STATUS {

  NB_DONE,
  NB_NOTDONE,
  NB_TOUCHED
};

#define HASH_SIZE                               1000
#define HASH_FACTOR                             .618033988              /* 0.5*(sqrt(5)-1)			*/
#define HASH_ADDRESS(k)                 floor(HASH_SIZE*(k*HASH_FACTOR-floor(k*HASH_FACTOR)))

#define ER_NSONS(p)                             ((p)->nsons)
#define ER_NCO(p,i)                             ((p)->nco[i])
#define ER_DCO(p,i)                             ((p)->dco[i])

#define HR_ID(p)                                ((p)->id)
#define HR_KEY(p)                               ((p)->key)
#define HR_TAG(p)                               ((p)->tag)
#define HR_NEXT(p)                              ((p)->next)
#define HR_ERULE(p)                             (&((p)->erule))
#define HR_NSONS(p)                             ((p)->erule.nsons)
#define HR_NCO(p,i)                             ((p)->erule.nco[i])
#define HR_DCO(p,i)                             ((p)->erule.dco[i])
#define HR_OCOPTR(p)                    ((p)->erule.dco+HR_NSONS(p))
#define HR_OCO(p,i)                             ((p)->erule.dco[HR_NSONS(p)+i])

/* debug macros */
#define ELEM_INFO(f,c,e,ns)             f ": %sconsiderd elem%ld tag%d refine%ld nsons%d (actually %d)\n",                      \
  (c) ? "not " : "",                                                                                                                         \
  (long)ID(e),                                                                                                                            \
  (int)TAG(e),                                                                                                                            \
  (long)REFINE(e),                                                                                                                        \
  (int)NSONS(e),                                                                                                                          \
  (int)ns

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {

  SHORT nsons;
  SHORT nco[MAX_SONS];
  DOUBLE dco[MAX_SONS];

} ERULE;

typedef INT HRID;

struct hashed_rule {

  HRID id;
  DOUBLE key;
  SHORT tag;
  struct hashed_rule *next;

  ERULE erule;
  /* will be allocated with additional doubles for ordered sons			*/
  /* CAUTION: storage occupied may be < sizeof(HRULE)						*/
};

typedef struct hashed_rule HRULE;
typedef REFRULE URULE;

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

static struct {

  HRULE **hash_table;
  HRULE **hrule[TAGS];
  long maxrule[TAGS];
  long maxrules;
  long nelem_inspected[TAGS];
  long nelems_inspected;
  long nelem_not_inspected[TAGS];
  long nelems_not_inspected;
  HEAP *heap;

        #ifdef ModelP
  long if_elems;
  ERULE *interface_rules;
        #endif

} global;

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*
        storage allocation stages:

        Notation:	M - Mark(Tmp)Mem
                                G - Get(Tmp)Mem
                                R - Release(Tmp)Mem


        ------------------------+-+-+-+-+---+-+-+-+-+
   |BOTTOM	|	|  TOP	|
        ------------------------+-+-+-+-+---+-+-+-+-+
                stack level			|0|1|2|3|	|0|1|2|3|
        ------------------------+-+-+-+-+---+-+-+-+-+
        (1) for hrules			|M| | | |	| | | | |
        (2) hash table			| | | | |	| |M| | |
 | | | | |	| |G| | |
        (3) interface rules		| | | | |	| | |M| |
 | | | | |	| | |G| |
        (4) hrules				|G| | | |	| | | | |
        (5) free (3)			| | | | |	| | |R| |
        (6) further hrules		|G| | | |	| | | | |
        (6) hrule table			|G| | | |	| | | | |
        (7) free (2)			| | | | |	| |R| | |
        (7) mrule table (rm+er)	| | | | |	|G| | | |
        (9) free hrules			|R| | | |	| | | | |
        ------------------------+-+-+-+-+---+-+-+-+-+
 */
/****************************************************************************/


/****************************************************************************/
/*D
   name - short_description

   SYNOPSIS:

   PARAMETERS:
   .  par - meaning

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

static DOUBLE Corner2DCorners (INT n, SHORT corners[])
{
  DOUBLE dco = corners[0];
  int i;

  for (i=1; i<n; i++)
  {
    dco *= MAX_REFINED_CORNERS_DIM;
    dco += corners[i];
  }
  return (dco);
}

static void DCorners2Corners (INT n, DOUBLE dco, SHORT corners[])
{
  int i;

  for (i=n-1; i>=0; i--)
  {
    DOUBLE x = floor(dco/((DOUBLE)MAX_REFINED_CORNERS_DIM));
    corners[i] = (SHORT)(dco-x*MAX_REFINED_CORNERS_DIM);
    dco = x;
  }
  return;
}

static int CornerCompare (const SHORT *c0, const SHORT *c1)
{
  if (*c0>*c1)
    return ( 1);
  else
    return (-1);
}

static int SonCompare (const DOUBLE *s0, const DOUBLE *s1)
{
  if (*s0>*s1)
    return ( 1);
  else
    return (-1);
}

static void FillOrderedSons (const ERULE *er, DOUBLE oco[])
{
  SHORT corners[MAX_CORNERS_OF_ELEM_DIM];
  int s;

  /* sort corners of each son */
  for (s=0; s<ER_NSONS(er); s++)
  {
    DCorners2Corners(ER_NCO(er,s),ER_DCO(er,s),corners);
    qsort(corners,ER_NCO(er,s),sizeof(*corners),
          (int (*)(const void *, const void *))CornerCompare);
    oco[s] = Corner2DCorners(ER_NCO(er,s),corners);
  }

  /* sort sons */
  qsort(oco,ER_NSONS(er),sizeof(*oco),
        (int (*)(const void *, const void *))SonCompare);

  return;
}

static INT Hash_Init (int MarkKey)
{
  int i;

  global.hash_table = (HRULE**) GetTmpMem(global.heap,HASH_SIZE*sizeof(HRULE*),MarkKey);
  if (global.hash_table==NULL)
    REP_ERR_RETURN(1);

  for (i=0; i<HASH_SIZE; i++)
    global.hash_table[i] = NULL;

  return (0);
}

static HRID Hash_InsertRule (INT etag, INT key, const ERULE *er, const DOUBLE oco[], HRULE **next_handle)
{
  size_t size = sizeof(HRULE)                           /* full HRULE				*/
                +sizeof(DOUBLE)*
                (2*ER_NSONS(er)                                 /* #DOUBLEs needed			*/
                 -MAX_SONS);                                            /* #DOUBLEs at end of HRULE	*/
  HRULE *hr       = (HRULE*) GetMem(global.heap,size,FROM_BOTTOM);
  HRID id         = global.maxrule[etag]++;

  if (hr==NULL)
    REP_ERR_RETURN (-1);

  /* insert in list */
  HR_NEXT(hr) = *next_handle;
  *next_handle = hr;

  /* init */
  HR_ID(hr)               = id;
  HR_KEY(hr)              = key;
  HR_TAG(hr)              = etag;
  memcpy(HR_ERULE(hr),    er,             sizeof(ERULE)+sizeof(DOUBLE)*(ER_NSONS(er)-MAX_SONS));
  memcpy(HR_OCOPTR(hr),   oco,    sizeof(DOUBLE)*ER_NSONS(er));

  return (id);
}

static INT SonsAreEqual (INT nsons, const DOUBLE oco[], const HRULE *hr)
{
  if (nsons!=HR_NSONS(hr))
    return (NO);
  else
  {
    int s;
    for (s=0; s<nsons; s++)
      if (oco[s]!=HR_OCO(hr,s))
        return (NO);
    return (YES);
  }
}

static HRID GetRuleID
(
        #ifdef Debug
  ELEMENT *elem,
        #endif
  INT etag,
  const ERULE *er
)
{
  DOUBLE key = 0;
  HRULE *hr;
  DOUBLE oco[MAX_SONS];
  int s,h;

  global.nelem_inspected[etag]++;
  global.nelems_inspected++;

  PRINTDEBUG(gm,2,(ELEM_INFO("GetRuleID",YES,elem,ER_NSONS(er))));

  FillOrderedSons(er,oco);

  for (s=0; s<ER_NSONS(er); s++)
    key += oco[s];

  /* TODO (HRR 971211): use tag also for key? */

  h = HASH_ADDRESS(key);

  hr = global.hash_table[h];
  if (hr==NULL)
    return (Hash_InsertRule(etag,key,er,oco,&(global.hash_table[h])));

  for (; hr!=NULL; hr=HR_NEXT(hr))
    if (key==HR_KEY(hr) && etag==HR_TAG(hr) && SonsAreEqual(ER_NSONS(er),oco,hr))
      return (HR_ID(hr));
    else if ((HR_NEXT(hr)==NULL) || (key<HR_KEY(HR_NEXT(hr))))
      break;

  return (Hash_InsertRule(etag,key,er,oco,&HR_NEXT(hr)));
}

static INT RuleCompare (int id, const URULE *ur, const ERULE *er)
{
  const int ns0   = NSONS_OF_RULE(ur);
  const int ns1   = ER_NSONS(er);
  int s0,s1;

  if (ns0!=ns1)
  {
    PrintDebug("CompareRules: rules not equal: elem%d (ns0=%d != ns1=%d)\n",id,ns0,ns1);
    return (NO);
  }

  /* compare sons */
  for (s0=0; s0<ns0; s0++)
  {
    const SONDATA *son0     = SON_OF_RULE(ur,s0);
    int nco0                        = CORNERS_OF_TAG(SON_TAG(son0));

    for (s1=0; s1<ns1; s1++)
    {
      int nco1 = ER_NCO(er,s1);

      if (nco0==nco1)
      {
        int corners_match = TRUE;
        int i;

        for (i=0; i<nco0; i++)
        {
          int co0 = SON_CORNER(son0,i);
          SHORT corners[MAX_CORNERS_OF_ELEM_DIM];
          int j;

          DCorners2Corners(ER_NCO(er,s1),ER_DCO(er,s1),corners);

          for (j=0; j<nco1; j++)
            if (co0==corners[j])
              break;
          if (j>=nco1)
          {
            corners_match = FALSE;
            break;
          }
        }
        if (corners_match)
          break;
      }
    }
    if (s1>=ns1)
    {
      PrintDebug("CompareRules: rules not equal: elem%d (son=%d of urule not in erule)\n",id,s0);
      return (NO);
    }
  }
  return (YES);
}

static INT ExtractERule (ELEMENT *elem, ERULE *er)
{
  int nsons = NSONS(elem);
  ELEMENT *sons[MAX_SONS];
  NODE *nodes[MAX_REFINED_CORNERS_DIM];
  int s;

  if (GetNodeContext(elem,nodes))
    REP_ERR_RETURN(1);

  /* get sons with master priority */
  if (GetSons(elem,sons))
    REP_ERR_RETURN(1);

  ER_NSONS(er) = nsons;
  for (s=0; s<nsons; s++)
  {
    ELEMENT *son = sons[s];
    int coe = CORNERS_OF_ELEM(son);
    SHORT corners[MAX_CORNERS_OF_ELEM_DIM];
    int k,j;

    /* corners and dcorners */
    ER_NCO(er,s) = coe;
    for (j=0; j<coe; j++)
    {
      NODE *node = CORNER(son,j);

      for (k=0; k<MAX_REFINED_CORNERS_DIM; k++)
        if (node==nodes[k])
          break;
      ASSERT(k<MAX_REFINED_CORNERS_DIM);
      corners[j] = k;
    }
    ER_DCO(er,s) = Corner2DCorners(CORNERS_OF_ELEM(son),corners);
  }

        #ifdef __DEBUG_ER__
  if (!(REFINE(elem)==COPY && NSONS(elem)>1))
  {
    URULE *ur = MARK2RULEADR(elem,REFINE(elem));

    RuleCompare(ID(elem),ur,er);
  }
        #endif

  return (0);
}

#ifdef ModelP
static int CountIFElements (DDD_OBJ obj)
{
  ELEMENT *elem = (ELEMENT*) obj;

  if (HAS_NO_RULE(elem))
  {
    SETTHEFLAG(elem,FALSE);
    return (0);
  }
  SETTHEFLAG(elem,TRUE);

  ASSERT(global.if_elems<MAX_IFID);

  SETREFINE(elem,global.if_elems);
  global.if_elems++;

  return (0);
}

static int InitMasterRules (DDD_OBJ obj)
{
  ELEMENT *elem   = (ELEMENT*) obj;

  if (THEFLAG(elem))
  {
    ERULE *er               = global.interface_rules+REFINE(elem);
    return (ExtractERule((ELEMENT*)obj,er));
  }
  else
    return (0);
}

static int Gather_ERULE (DDD_OBJ obj, void *data)
{
  ELEMENT *elem   = (ELEMENT*) obj;

  if (THEFLAG(elem))
    return (ExtractERule(elem,(ERULE*)data));
  else
    return (0);
}

static int Scatter_ERULE (DDD_OBJ obj, void *data)
{
  ELEMENT *elem   = (ELEMENT*)obj;

  if (THEFLAG(elem))
  {
    ERULE *er               = global.interface_rules+REFINE(elem);
    memcpy(er,data,sizeof(ERULE));
  }

  return (0);
}

static int Scatter_partial_ERULE (DDD_OBJ obj, void *data)
{
  ELEMENT *elem   = (ELEMENT*)obj;

  if (THEFLAG(elem))
  {
    ERULE *er_ma    = global.interface_rules+REFINE(elem);
    ERULE *er_gh    = (ERULE*)data;
    int s_ma                = ER_NSONS(er_ma);
    int s_gh;

    /* put sons at and of rule */
    for (s_gh=0; s_gh<ER_NSONS(er_gh); s_gh++)
    {
      ASSERT(s_ma<MAX_SONS);

      ER_NCO(er_ma,s_ma) = ER_NCO(er_gh,s_gh);
      ER_DCO(er_ma,s_ma) = ER_DCO(er_gh,s_gh);
      s_ma++;
    }
    ER_NSONS(er_ma) = s_ma;
  }

  return (0);
}

static int ExtractInterfaceERule (DDD_OBJ obj)
{
  ELEMENT *elem   = (ELEMENT*)obj;

  if (THEFLAG(elem))
  {
    ERULE *er               = global.interface_rules+REFINE(elem);
    int id                  = GetRuleID(
                                                                        #ifdef Debug
      elem,
                                                                        #endif
      TAG(elem),er);

    ASSERT(id<MAX_HRID);

    SETREFINE(elem,id);
  }
  else
  {
    global.nelem_not_inspected[TAG(elem)]++;
    global.nelems_not_inspected++;

    IFDEBUG(gm,2)
    int nsons;
    ELEMENT *sons[MAX_SONS];

    if (GetSons(elem,sons))
      REP_ERR_RETURN(1);
    for (nsons=0; nsons<MAX_SONS; nsons++)
      if (sons[nsons]==NULL)
        break;
    PrintDebug(ELEM_INFO("ExtractInterfaceERule",NO,elem,nsons));
    ENDDEBUG
  }

  SETUSED(elem,TRUE);

  return (0);
}

static INT ExtractInterfaceRules (MULTIGRID *mg)
{
  int lev;

  /* TODO (HRR 971211): don't include TOPLEVEL (no elem refined there) */
  for (lev=0; lev<=TOPLEVEL(mg); lev++)
  {
    GRID *grid = GRID_ON_LEVEL(mg,lev);

    /* count interface master and vhghost elements */
    global.if_elems = 0;
    DDD_IFAExecLocal(ElementVHIF, GRID_ATTR(grid), CountIFElements);

    if (global.if_elems>0)
    {
      INT MarkKey;

      PRINTDEBUG(gm,1,("ExtractInterfaceRules: interface elements to consider on level %d: %ld\n",lev,global.if_elems));

      if (MarkTmpMem(global.heap,&MarkKey))
        REP_ERR_RETURN(1);

      /* REMARK (HRR 971207): storage could be Marked and Released in pieces if necessary
                                                      to allow BOTTOM-storage to grow */
      global.interface_rules = GetTmpMem(global.heap,global.if_elems*sizeof(ERULE),MarkKey);
      if (global.interface_rules==NULL)
        REP_ERR_RETURN(1);

      /* init rules of masters */
      DDD_IFAExecLocal(ElementIF, GRID_ATTR(grid), InitMasterRules);

      /* communicate VHghosts --> master */
      DDD_IFAOneway(ElementVHIF, GRID_ATTR(grid), IF_BACKWARD, sizeof(ERULE),
                    Gather_ERULE, Scatter_partial_ERULE);

      /* communicate master --> VHghosts */
      DDD_IFAOneway(ElementVHIF, GRID_ATTR(grid), IF_FORWARD, sizeof(ERULE),
                    Gather_ERULE, Scatter_ERULE);

      /* extract rules from interface elements */
      DDD_IFAExecLocal(ElementVHIF, GRID_ATTR(grid), ExtractInterfaceERule);

      IFDEBUG(gm,1)
      long N_er = 0;
      int tag;

      PrintDebug("ExtractInterfaceRules: level %d: HeapFree is %ld bytes\n",lev,(long)HeapFree(global.heap));

      PrintDebug("ExtractInterfaceRules: number of refrules extracted from interface on level %d:\n",lev);
      for (tag=0; tag<TAGS; tag++)
      {
        long n_er = global.maxrule[tag] - UGMAXRULE(tag);
        N_er += n_er;
        PrintDebug("tag %d: %ld extracted rules\n",tag,n_er);
      }
      PrintDebug("total: %ld extracted rules\n",N_er);
      ENDDEBUG

      if (ReleaseTmpMem(global.heap,MarkKey))
        REP_ERR_RETURN(1);
    }
  }

  return (0);
}
#endif

static INT ExtractRules (MULTIGRID *mg)
{
  ELEMENT *elem;
  ERULE er;
  HRID id;
  INT MarkKey;
  int l,tag,h,maxrules;

  /* for hash table */
  if (MarkTmpMem(global.heap,&MarkKey))
    REP_ERR_RETURN(1);

  if (Hash_Init(MarkKey))
    REP_ERR_RETURN(1);

  global.nelems_inspected = global.nelems_not_inspected = 0;
  for (tag=0; tag<TAGS; tag++)
    global.nelem_inspected[tag] = global.nelem_not_inspected[tag] = 0;

  /* loop elements and extract rules */
        #ifdef ModelP
  /* TODO (HRR 971211): don't include TOPLEVEL (no elem refined there) */
  for (l=0; l<=TOPLEVEL(mg); l++)
    for (elem=FIRSTELEMENT(GRID_ON_LEVEL(mg,l)); elem!=NULL; elem=SUCCE(elem))
      SETUSED(elem,FALSE);
  if (ExtractInterfaceRules(mg))
    REP_ERR_RETURN(1);
        #endif

  /* TODO (HRR 971211): don't include TOPLEVEL (no elem refined there) */
  for (l=0; l<=TOPLEVEL(mg); l++)
    for (elem=FIRSTELEMENT(GRID_ON_LEVEL(mg,l)); elem!=NULL; elem=SUCCE(elem))
                        #ifdef ModelP
      if (!USED(elem))
                        #endif
      if (HAS_NO_RULE(elem))
      {
        if (ExtractERule(elem,&er))
          REP_ERR_RETURN(1);
        id = GetRuleID(
                                                                #ifdef Debug
          elem,
                                                                #endif
          TAG(elem),&er);
        ASSERT(id<MAX_HRID);
        SETREFINE(elem,id);
      }
      else
      {
        global.nelem_not_inspected[TAG(elem)]++;
        global.nelems_not_inspected++;
        IFDEBUG(gm,2)
        int nsons;
        ELEMENT *sons[MAX_SONS];

        if (GetSons(elem,sons))
          REP_ERR_RETURN(1);
        for (nsons=0; nsons<MAX_SONS; nsons++)
          if (sons[nsons]==NULL)
            break;
        PrintDebug(ELEM_INFO("ExtractRules",NO,elem,nsons));
        if (NSONS(elem)!=nsons)
          PrintDebug("------ ERROR: NSONS!=nsons -------\n");
        ENDDEBUG
      }

  global.maxrules = maxrules = 0;
  for (tag=0; tag<TAGS; tag++)
  {
    global.maxrules += global.maxrule[tag];
    maxrules += global.maxrule[tag] - UGMAXRULE(tag);
  }

  if (maxrules>0)
  {
    int n = 0;
    int max_list_len = 0;

    /* make tables of subsequent IDs */
    global.hrule[0] = (HRULE**) GetMem(global.heap,global.maxrules*sizeof(HRULE*),FROM_BOTTOM);
    if (global.hrule[0]==NULL)
      REP_ERR_RETURN(1);
    for (tag=1; tag<TAGS; tag++)
      global.hrule[tag] = global.hrule[tag-1]+global.maxrule[tag-1];

    PRINTDEBUG(gm,1,("ExtractRules: after hrule table allocation: HeapFree is %ld bytes\n",(long)HeapFree(global.heap)));

    for (h=0; h<HASH_SIZE; h++)
    {
      int list_len = 0;
      HRULE *hr;

      for (hr=global.hash_table[h]; hr!=NULL; hr=HR_NEXT(hr))
      {
        global.hrule[HR_TAG(hr)][HR_ID(hr)] = hr;
        list_len++;
        n++;
      }
      max_list_len = MAX(max_list_len,list_len);
    }
    ASSERT(maxrules==n);

    PRINTDEBUG(gm,1,("ExtractRules: max members of hash table entries is %d\n",max_list_len));
  }
  else
    PRINTDEBUG(gm,1,("ExtractRules: realized rules were completely covered by rm rules\n"));

  if (ReleaseTmpMem(global.heap,MarkKey))
    REP_ERR_RETURN(1);

  return (0);
}

static void FindPathForNeighbours (MGIO_RR_RULE *rule, SHORT myID, SHORT Status[MAX_SONS])
{
  SHORT i,nbID;

  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    if (((nbID=rule->sons[myID].nb[i])<FATHER_SIDE_OFFSET) && (Status[nbID]==NB_NOTDONE))
    {
      int *nbPath = &(rule->sons[nbID].path);
      SHORT nbPathDepth;

      /* copy myPath to nbPath */
      *nbPath = rule->sons[myID].path;

      /* complete nbPath */
      nbPathDepth = PATHDEPTH(*nbPath);
      SETNEXTSIDE(*nbPath,nbPathDepth,i);
      SETPATHDEPTH(*nbPath,++nbPathDepth);
      Status[nbID] = NB_TOUCHED;
    }

  /* recursive call for NB_TOUCHED sons */
  for (nbID=1; nbID<rule->nsons; nbID++)
    if (Status[nbID]==NB_TOUCHED)
    {
      Status[nbID] = NB_DONE;
      FindPathForNeighbours(rule,nbID,Status);
    }

  return;
}

static void FillSonPaths (MGIO_RR_RULE *rule)
{
  SHORT Status[MAX_SONS];
  int i;

  /* TODO (HRR 971211): debug recursive path construction
     (taken from GenerateRules but not used by ugio) */

  /* son 0 has trivial path */
  Status[0] = NB_DONE;
  for (i=1; i<rule->nsons; i++)
    Status[i] = NB_NOTDONE;

  /* start recursion with son 0 */
  FindPathForNeighbours(rule,0,Status);

  return;
}

static INT GetFSidesOfCorners (int tag, int n, SHORT corners[MAX_CORNERS_OF_SIDE], SHORT corner_on_side[MAX_CORNERS_OF_SIDE][MAX_SIDES_OF_ELEM])
{
  int coe = CORNERS_OF_TAG(tag);
  int eoe = EDGES_OF_TAG(tag);
  int soe = SIDES_OF_TAG(tag);
  int co,side;

  for (co=0; co<n; co++)
    for (side=0; side<soe; side++)
      corner_on_side[co][side] = FALSE;

  for (co=0; co<n; co++)
    if (corners[co]==coe+CENTER_NODE_INDEX_TAG(tag))
    {
      /* center node: can not be part of a father side */
      return (NO);
    }
    else if (corners[co]<coe)
    {
      /* father corner */
      int fco = corners[co];

      for (side=0; side<soe; side++)
        if (CORNER_OF_SIDE_INV_TAG(tag,side,fco)>=0)
          corner_on_side[co][side] = TRUE;
    }
    else if (corners[co]<(coe+eoe))
    {
      /* edge mid node */
      int ed = corners[co]-coe;
      int i;

                        #ifdef __TWODIM__
      corner_on_side[co][ed] = TRUE;
                        #else
      for (i=0; i<MAX_SIDES_OF_EDGE; i++)
      {
        int sd = SIDE_WITH_EDGE_TAG(tag,ed,i);
        if (sd>=0)
          corner_on_side[co][sd] = TRUE;
      }
                        #endif
    }
                #ifdef __THREEDIM__
  else if (corners[co]<(coe+eoe+soe))
  {
    /* side mid */
    int sd = corners[co]-(coe+eoe);
    corner_on_side[co][sd] = TRUE;
  }
                #endif
  else
    ASSERT(FALSE);                              /* Huh??? */

  return (YES);
}

static INT GetCommonFSide (int nco, int nsi, SHORT corner_on_side[MAX_CORNERS_OF_SIDE][MAX_SIDES_OF_ELEM])
{
  int i,side;

  for (side=0; side<nsi; side++)
  {
    for (i=0; i<nco; i++)
      if (!corner_on_side[i][side])
        break;
    if (i>=nco)
      return (side);
  }
  return (nsi);
}

static INT IsOnFatherSide (int tag, int nsco, SHORT sco[], SHORT *nb)
{
  SHORT sco_on_side[MAX_CORNERS_OF_SIDE][MAX_SIDES_OF_ELEM];

  if (GetFSidesOfCorners(tag,nsco,sco,sco_on_side))
  {
    int fside = GetCommonFSide(nsco,SIDES_OF_TAG(tag),sco_on_side);

    if (fside<SIDES_OF_TAG(tag))
    {
      *nb = FATHER_SIDE_OFFSET+fside;
      return (YES);
    }
  }
  return (NO);
}

static INT SidesMatch (int nsco, SHORT sco0[], SHORT sco1[])
{
  int i;

  /* try each permutation of first with reverse order of second */
  for (i=0; i<nsco; i++)
  {
    int match = TRUE;
    int j;

    for (j=0; j<nsco; j++)
      if (sco0[(i+j)%nsco]!=sco1[nsco-j-1])
      {
        match = FALSE;
        break;
      }
    if (match)
      return (YES);
  }
  return (NO);
}

static void HRule2Mrule (const HRULE *hr, MGIO_RR_RULE *mr)
{
  int coe = CORNERS_OF_TAG(HR_TAG(hr));
  int s,s0;

  /* extracted rules are always irregular */
  /* TODO (HRR 971208): is that really ok (actually not used)? */
  mr->class = GREEN_CLASS;
  mr->nsons = HR_NSONS(hr);

  {int i; for (i=0; i<MGIO_MAX_NEW_CORNERS; i++)
     mr->pattern[i] = 0;}

  /* son tags, corners, sonandnode, pattern */
  for (s=0; s<mr->nsons; s++)
  {
    struct mgio_sondata *sonData = &(mr->sons[s]);
    int nco = HR_NCO(hr,s);
    int j;
    int sd0;

    for (sd0=0; sd0<MAX_SIDES_OF_ELEM; sd0++)
      sonData->nb[sd0] = NOTDONE;

    sonData->tag = REF2TAG(nco);
    DCorners2Corners(nco,HR_DCO(hr,s),sonData->corners);
    for (j=0; j<nco; j++)
    {
      int newco = sonData->corners[j]-coe;

      if (newco<0)
        continue;

      /* TODO (HRR 971208): what about center node pattern for tets? */
      mr->pattern[newco] = 1;
      mr->sonandnode[newco][0] = s;
      mr->sonandnode[newco][1] = j;
    }
  }

  /* son nb */
  for (s0=0; s0<mr->nsons; s0++)
  {
    struct mgio_sondata *sonData0 = mr->sons+s0;
    int sd0;

    for (sd0=0; sd0<MAX_SIDES_OF_ELEM; sd0++)
      if (sonData0->nb[sd0] == NOTDONE)
      {
        int stag0 = sonData0->tag;
        int nsco0 = CORNERS_OF_SIDE_TAG(stag0,sd0);
        SHORT sco0[MAX_CORNERS_OF_SIDE];

        {int j; for (j=0; j<nsco0; j++)
           sco0[j] = sonData0->corners[CORNER_OF_SIDE_TAG(stag0,sd0,j)];}

        if (!IsOnFatherSide(HR_TAG(hr),nsco0,sco0,&(sonData0->nb[sd0])))
        {
          int found = NO;
          int s1;

          /* find matching side of other son */
          for (s1=s0+1; s1<mr->nsons; s1++)
          {
            struct mgio_sondata *sonData1 = mr->sons+s1;
            int sd1;

            for (sd1=0; sd1<MAX_SIDES_OF_ELEM; sd1++)
            {
              int stag1 = sonData1->tag;
              int nsco1 = CORNERS_OF_SIDE_TAG(stag1,sd1);

              if (nsco0==nsco1)
              {
                SHORT sco1[MAX_CORNERS_OF_SIDE];

                {int j; for (j=0; j<nsco1; j++)
                   sco1[j] = sonData1->corners[CORNER_OF_SIDE_TAG(stag1,sd1,j)];}

                if (SidesMatch(nsco0,sco0,sco1))
                {
                  /* neighbour found */
                  sonData0->nb[sd0] = s1;
                  sonData1->nb[sd1] = s0;

                  found = YES;
                  break;
                }
              }
            }
            if (found)
              break;
          }
          ASSERT(found);
        }
      }
  }

  /* son path */
  FillSonPaths(mr);

  return;
}

static void URule2Mrule (const URULE *ur, MGIO_RR_RULE *mr)
{
  int j,k;

  mr->class = ur->class;
  mr->nsons = ur->nsons;
  for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
    mr->pattern[j] = ur->pattern[j];
  for (j=0; j<MGIO_MAX_NEW_CORNERS; j++)
  {
    mr->sonandnode[j][0] = ur->sonandnode[j][0];
    mr->sonandnode[j][1] = ur->sonandnode[j][1];
  }
  for (j=0; j<mr->nsons; j++)
  {
    struct mgio_sondata *sonData = &(mr->sons[j]);

    sonData->tag = ur->sons[j].tag;
    for (k=0; k<MGIO_MAX_CORNERS_OF_ELEM; k++)
      sonData->corners[k] = ur->sons[j].corners[k];
    for (k=0; k<MGIO_MAX_SIDES_OF_ELEM; k++)
      sonData->nb[k] = ur->sons[j].nb[k];
    sonData->path = ur->sons[j].path;
  }
  return;
}

static void WriteDebugInfo (void)
{
  long N_rm=0,N_er=0;
  int tag;

  /* number of rules (rm+er) */
  PrintDebug("------------- Write_RefRules statistics --------------\n");
  PrintDebug("number of refrules:\n");
  for (tag=0; tag<TAGS; tag++)
  {
    long n_rm = UGMAXRULE(tag);
    long n_er = global.maxrule[tag]-UGMAXRULE(tag);
    N_rm += n_rm;
    N_er += n_er;
    PrintDebug("tag %d: %ld rm rules, %ld extracted rules (elems inspected: %ld yes %ld no)\n",
               tag,
               n_rm,
               n_er,
               global.nelem_inspected[tag],
               global.nelem_not_inspected[tag]);
  }
  PrintDebug("total: %ld rm rules, %ld extracted rules (elems inspected: %ld yes %ld no)\n",
             N_rm,
             N_er,
             global.nelem_inspected[tag],
             global.nelem_not_inspected[tag]);

  PrintDebug("------------------------------------------------------\n");

  return;
}

INT NEW_Write_RefRules (MULTIGRID *mg, INT RefRuleOffset[], INT MarkKey, MGIO_RR_RULE **mrule_handle)
{
  MGIO_RR_GENERAL rr_general;
  MGIO_RR_RULE *mrule;
  int BotMarkKey;
  int tag;

  if (mg==NULL)
    REP_ERR_RETURN(1);

  global.heap = MGHEAP(mg);
  PRINTDEBUG(gm,1,("Write_RefRules (er): before any allocation: HeapFree is %ld bytes\n",(long)HeapFree(global.heap)));
  if (Mark(global.heap,FROM_BOTTOM,&BotMarkKey))
    REP_ERR_RETURN(1);

  /* init rule counters (continue with last rule IDs of rm) */
  /* TODO (HRR 971207): this is important for matching existing rules after loading a grid (coarsening)
                                            but CAUTION: refine should never address a rule beyond rm rules! */
  for (tag=0; tag<TAGS; tag++)
    global.maxrule[tag] = UGMAXRULE(tag);

        #if (defined __THREEDIM__) || (defined __DEBUG_ER__)
  /* incomplete refrule set: extract non-existing rules that are realized in mg */
  if (ExtractRules(mg))
    REP_ERR_RETURN(1);
        #endif

  /* write refrules general */
  RefRuleOffset[0] = 0;
  for (tag=0; tag<TAGS; tag++)
  {
    if (tag>0) RefRuleOffset[tag] = RefRuleOffset[tag-1] + global.maxrule[tag-1];
    rr_general.RefRuleOffset[tag] = RefRuleOffset[tag];
  }
  rr_general.nRules = global.maxrules;
  if (Write_RR_General(&rr_general))
    REP_ERR_RETURN(1);

  /* allocate MGIO_RR_RULE table (will stay in scope beyond this module) */
  *mrule_handle = (MGIO_RR_RULE*) GetTmpMem(global.heap,global.maxrules*sizeof(MGIO_RR_RULE),MarkKey);
  if (*mrule_handle==NULL)
    REP_ERR_RETURN(1);

  PRINTDEBUG(gm,1,("Write_RefRules (er): after mrule table allocation: HeapFree is %ld bytes\n",(long)HeapFree(global.heap)));

  /* write refrules */
  mrule = *mrule_handle;
  for (tag=0; tag<TAGS; tag++)
  {
    int i;

    /* rm rules */
    for (i=0; i<UGMAXRULE(tag); i++)
    {
      ASSERT(RefRules[tag]+i!=NULL);
      URule2Mrule(RefRules[tag]+i,mrule);
      mrule++;
    }

    /* er rules */
    for (; i<global.maxrule[tag]; i++)
    {
      ASSERT(global.hrule[tag][i]!=NULL);
      HRule2Mrule(global.hrule[tag][i],mrule);
      mrule++;
    }

  }
  Write_RR_Rules(global.maxrules,*mrule_handle);

  /* free hrules and hrule table */
  if (Release(global.heap,FROM_BOTTOM,BotMarkKey))
    REP_ERR_RETURN(1);

  IFDEBUG(gm,1)
  WriteDebugInfo();
  ENDDEBUG

  PRINTDEBUG(gm,1,("Write_RefRules (er): when finished (storage occupied by MGIO_RR_RULE list): HeapFree is %ld bytes\n",(long)HeapFree(global.heap)));

  return (0);
}

INT ResetRefineTagsBeyondRuleManager (MULTIGRID *mg)
{
  ELEMENT *elem;
  int l,n=0;

  /* TODO (HRR 971211): don't include TOPLEVEL (no elem refined there) */
  for (l=0; l<=TOPLEVEL(mg); l++)
    for (elem=PFIRSTELEMENT(GRID_ON_LEVEL(mg,l)); elem!=NULL; elem=SUCCE(elem))
      if (BEYOND_UG_RULES(elem))
      {
        SETREFINE(elem,COPY);
        n++;
      }
  PRINTDEBUG(gm,1,("ResetRefineTags: done (for %d elements)\n",n));

  return (0);
}
