// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ooppcc.h                                                      */
/*                                                                          */
/* Purpose:   macros and templates for oopp-container classes.              */
/*            (object-oriented pre-processing container classes).           */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: birken@ica3.uni-stuttgart.de                           */
/*            phone: 0049-(0)711-685-7007                                   */
/*            fax  : 0049-(0)711-685-7000                                   */
/*                                                                          */
/* History:   970902 kb  begin                                              */
/*            971015 kb  added ArrayOf, BTreeOf                             */
/*            971020 kb  added SetOf                                        */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */



/****************************************************************************/
#ifdef PtrOf
/****************************************************************************/

#define Ptr(T)  CN(CCAT(T,Ptr))

/* class name for ptr */
#define CPtr    CCAT(PtrOf,Ptr)

/* an object's pointer class is simply a C pointer to the object */
typedef PtrOf * CPtr;


/****************************************************************************/
#undef CPtr
#undef PtrOf
#endif
/****************************************************************************/




/****************************************************************************/
#ifdef ArrayOf
/****************************************************************************/

#define Array(T)   CN(CCAT(T,Array))

/* class name for array */
#define CArray     CCAT(ArrayOf,Array)


#ifndef ArrayAllocate
/* define default allocate function */
#define ArrayAllocate OO_Allocate
#endif

#ifndef ArrayFree
/* define default free function */
#define ArrayFree OO_Free
#endif


/*** array class ***/

#define ClassName CArray
Class_Data_Begin
CN(ArrayOf) *data;
int size;
int used;
Class_Data_End
Method_New_     (int _NEWPARAMS);
void         Method(Free)    (DefThis);
CN(ArrayOf) *Method(GetData) (DefThis);
int          Method(GetSize) (DefThis);


#ifdef ContainerImplementation
Method_New_ (int size _NEWPARAMS)
{
  Construct(This, _CHECKALLOC(This));

  if (size==0)
  {
    This->data = NULL;
  }
  else
  {
    This->data = (CN(ArrayOf) *)ArrayAllocate (sizeof(CN(ArrayOf))*size);
    _CHECKALLOC(This->data);
  }

  This->size = size;
  This->used = 0;

  return(This);
}


void Method(Free) (ParamThis)
{
#ifndef NoArrayFree
  if (This->data!=NULL)
    ArrayFree (This->data);
#endif

  Destruct(This);
}


CN(ArrayOf) *Method(GetData) (ParamThis)
{
  return(This->data);
}


int Method(GetSize) (ParamThis)
{
  return(This->size);
}

#endif

#undef ClassName


/****************************************************************************/
#undef ArrayAllocate
#undef ArrayFree
#ifdef NoArrayFree
        #undef NoArrayFree
#endif
#undef CArray
#undef ArrayOf
#endif
/****************************************************************************/



/****************************************************************************/
#ifdef SegmListOf
/****************************************************************************/

#define SegmList(T)  CN(CCAT(T,SegmList))

/* class names for segmlist and segment */
#define CSegm        CCAT(SegmListOf,Segm)
#define CSegmList    CCAT(SegmListOf,SegmList)


#ifndef SegmSize
/* define default size of each segment */
#define SegmSize 256
#endif





/*** one segment ***/

#define ClassName CSegm
Class_Data_Begin
CN(SegmListOf)  data[SegmSize];
int nItems;

ClassPtr next;
Class_Data_End
Method_New_   (_NEWPARAMS_OR_VOID);
void Method(Free)  (DefThis);
void Method(Print) (DefThis _PRINTPARAMS);


#ifdef ContainerImplementation

Method_New_ (_NEWPARAMS_OR_VOID)
{
  Construct(This, _CHECKALLOC(This));
  This->nItems = 0;
  return(This);
}

void Method(Free) (ParamThis)
{
  Destruct(This);
}

#endif

#undef ClassName





/*** the segm-list ***/

#define ClassName CSegmList
Class_Data_Begin
CSegm    *first;
int nItems;
int nSegms;

int nDiscarded;
Class_Data_End
Method_New_     (_NEWPARAMS_OR_VOID);
CN(SegmListOf) *Method(NewItem) (DefThis);
void    Method(DiscardItem)   (DefThis);
void    Method(Reset)         (DefThis);
int     Method(GetNItems)     (DefThis);
int     Method(GetNDiscarded) (DefThis);
void    Method(GetResources)  (DefThis, int *, int *, size_t *, size_t *);


#ifdef ContainerImplementation

Method_New_ (_NEWPARAMS_OR_VOID)
{
  Construct(This, _CHECKALLOC(This));
  This->first  = NULL;
  This->nItems = 0;
  This->nSegms = 0;
  This->nDiscarded = 0;
  return(This);
}


CN(SegmListOf) *Method(NewItem) (ParamThis)
{
  CN(SegmListOf) *item;
  CSegm          *segm = This->first;

  if (segm==NULL || segm->nItems==SegmSize)
  {
    segm = CALL(New,CSegm) ();
    segm->next = This->first;
    This->first = segm;
    This->nSegms++;
  }

  /* get item address in current segment */
  item = &(segm->data[segm->nItems++]);

  /* count item */
  This->nItems++;

  return(item);
}


int Method(GetNItems) (ParamThis)
{
  return(This->nItems);
}

int Method(GetNDiscarded) (ParamThis)
{
  return(This->nDiscarded);
}

void Method(DiscardItem) (ParamThis)
{
  assert(This!=NULL);
  assert(This->first!=NULL);
  assert(This->first->nItems > 0);

  This->first->nItems--;
  This->nItems--;

  This->nDiscarded++;
}


void Method(Reset) (ParamThis)
{
  CSegm *segm, *next;

  segm = This->first;
  next = NULL;

  while (segm!=NULL)
  {
    next = segm->next;
    CALL(CSegm,Free) (segm);

    segm = next;
  }

  This->first = NULL;
  This->nItems = 0;
  This->nSegms = 0;
  This->nDiscarded = 0;
}


void Method(GetResources) (ParamThis,
                           int *nSegms, int *nItems, size_t *alloc_mem, size_t *used_mem)
{
  size_t allocated=0, used=0;
  CSegm    *segm;

  for (segm=This->first; segm!=NULL; segm=segm->next)
  {
    /* compute memory usage */
    allocated += sizeof(CSegm);
    used += (sizeof(CSegm) -
             (sizeof(CN(SegmListOf)) * (SegmSize-segm->nItems)));
  }

  *nSegms    = This->nSegms;
  *nItems    = This->nItems;
  *alloc_mem = allocated;
  *used_mem  = used;
}

#endif

#undef ClassName



/****************************************************************************/
#undef CSegm
#undef CSegmList
#undef SegmSize
#undef SegmListOf
#endif
/****************************************************************************/



/****************************************************************************/
#ifdef ListOf
/****************************************************************************/

#define List(T)     CN(CCAT(T,List))

/* class names for list and list item */
#define CListItem  CCAT(ListOf,ListItem)
#define CList      CCAT(ListOf,List)


/* method pointers */
#define Find_Method    CCAT(ListOf,FindMethod)
typedef int (*Find_Method)(CN(ListOf) *);



/*** Item of list class ***/

#define ClassName CListItem
Class_Data_Begin
CN(ListOf)  *data;
ClassPtr next;
Class_Data_End
Method_New_   (CN(ListOf) * _NEWPARAMS);
void Method(Print) (DefThis _PRINTPARAMS);


#ifdef ContainerImplementation
Method_New_ (CN(ListOf) *data _NEWPARAMS)
{
  Construct(This, _CHECKALLOC(This));
  This->data = data;
  This->next = NULL;
  return(This);
}

void Method(Print) (ParamThis _PRINTPARAMS)
{
  CALL(ListOf,Print) (This->data _PRINTNEXT);
}
#endif
#undef ClassName




/*** List class ***/

#define ClassName CList
Class_Data_Begin
CN(CListItem) *first;
CN(CListItem) *last;
Class_Data_End
Method_New_   (_NEWPARAMS_OR_VOID);
void Method(Print) (DefThis _PRINTPARAMS);
ClassRef Method(Append) (DefThis, CN(ListOf) *);
CN(ListOf) *Method(Find)   (DefThis, Find_Method);


#ifdef ContainerImplementation
Method_New_ (_NEWPARAMS_OR_VOID)
{
  Construct(This, _CHECKALLOC(This));
  This->first   = NULL;
  This->last    = NULL;
  return(This);
}


void Method(Print) (ParamThis _PRINTPARAMS)
{
  CN(CListItem) *item;

  for(item=This->first; item!=NULL; item=item->next)
    CALL(CListItem,Print) (item _PRINTSAME);
}


ClassRef Method(Append) (ParamThis, CN(ListOf) *data)
{
  CN(CListItem) *item;

  if (data==NULL) return(This);

  /* create item and link it into list */
  item = CALL(New,CListItem) (data);
  if (This->first==NULL)
  {
    This->first = item;
    This->last  = item;
  }
  else
  {
    This->last->next = item;
    This->last = item;
  }

  return(This);
}


CN(ListOf) * Method(Find) (ParamThis, Find_Method find_data)
{
  CN(CListItem) *item;

  for (item=This->first; item!=NULL; item=item->next)
  {
    if ((*find_data)(item->data))
      return(item->data);
  }
  return(NULL);
}

#endif
#undef ClassName


/*** List iterator class ***/

/* iterator is simply a pointer to a CListItem */
#define ListIter(T) CN(CCAT(T,ListIter))
typedef CN (CListItem) * ListIter (ListOf);

#ifndef __ListIterator__
#define __ListIterator__

/* forward iterator */
#define ListFirst(list)   ((list)->first)
#define ListNext(iter)    ((iter)->next)
#define ListIsDone(iter)  ((iter)==NULL)
#define ListCurrent(iter) ((iter)->data)
#endif


/****************************************************************************/
#undef CListItem
#undef CList
#undef Find_Method
#undef ListOf
#endif
/****************************************************************************/




/****************************************************************************/
#ifdef BTreeOf
/****************************************************************************/


/* class names for btree and btree node */
#define CBTreeNode  CCAT(BTreeOf,BTreeNode)
#define CBTree      CCAT(BTreeOf,BTree)

/* this is the corresponding pointer-array class */
#define CPtrArray   CCAT(CN(BTreeOf),PtrArray)


/* method pointers */
#define Iterate_Method         CCAT(BTreeOf,IterateMethod)
typedef void (*Iterate_Method)(CN(BTreeOf) *);
/* not used
   #define Compare_Method    CCAT(BTreeOf,CompareMethod)
   typedef int (*Compare_Method)  (CN(BTreeOf) *, CN(BTreeOf) *);
 */


#ifndef BTreeOrder
/* define default order of BTree */
#define BTreeOrder 8
#endif


/* now the definitions for _all_ BTrees */
#ifndef BTreeGeneral
#define BTreeGeneral

enum BTreeConstant {
  BTREE_ERROR, BTREE_INSERTED, BTREE_FOUND, BTREE_SPLIT_ME
};

#endif





/*** node of BTree class ***/

#define ClassName CBTreeNode
Class_Data_Begin
int nSons;
ClassPtr sons[BTreeOrder+1];

CN(BTreeOf) *data[BTreeOrder];
Class_Data_End
/* no method declaration here, all methods are static */


#ifdef ContainerImplementation

static Method_New_ (CN(BTreeOf) *item,
                    ClassPtr son_l, ClassPtr son_r _NEWPARAMS)
{
  Construct(This, _CHECKALLOC(This));
  This->nSons = 2;
  This->sons[0] = son_l;
  This->sons[1] = son_r;
  This->data[0] = item;
  return(This);
}

static void Method(Print) (ParamThis _PRINTPARAMS)
{
  int i;

  if (This==NULL) return;

  _INDENT; fprintf(fp, "NODE nSons=%d\n", This->nSons);

  for(i=0; i<This->nSons-1; i++)
  {
    if (This->sons[i]!=NULL)
      Method(Print) (This->sons[i]  _PRINTNEXT);

    _INDENT1; CALL(BTreeOf,Print) (This->data[i]  _PRINTNEXT);
  }
  if (This->sons[i]!=NULL)
    Method(Print) (This->sons[i]  _PRINTNEXT);
}


static void Method(Free) (ParamThis)
{
  int i;

  /* free all son nodes */
  for(i=0; i<This->nSons; i++)
  {
    if (This->sons[i]!=NULL)
    {
      Method(Free) (This->sons[i]);
    }
  }

  /* destruct this node */
  Destruct(This);
}


static void Method(Split)  (ParamThis,
                            CN(BTreeOf) **split_item, ClassPtr *new_rnode)
{
  int split, l, r;
  Construct(rnode, _CHECKALLOC(rnode));

  /* compute split index, left part is at most one item smaller */
  split = ((BTreeOrder+1)>>1)-1;
  /*printf("\t\tsplit_index = %d\n", split);*/

  for(l=split+1, r=0; l<This->nSons; l++, r++)
  {
    rnode->sons[r] = This->sons[l];
    rnode->data[r] = This->data[l];
  }
  rnode->sons[r] = rnode->sons[l];

  rnode->nSons = This->nSons-split-1;
  This->nSons  = split+1;
  /*BTreeNodePrint(lnode,9);*/
  /*BTreeNodePrint(rnode,9);*/

  /* set return values */
  *split_item = This->data[split];
  *new_rnode = rnode;
}


static Method(FreeIndex)  (ParamThis, int idx)
{
  int j, n=This->nSons;

  This->sons[n] = This->sons[n-1];
  for(j=n-1; j>idx; j--)
  {
    This->sons[j] = This->sons[j-1];
    This->data[j] = This->data[j-1];
  }

  This->nSons++;
}


static int Method(Insert) (ParamThis,
                           /*Compare_Method cmp_func,*/
                           CN(BTreeOf) *item)
{
  int i, ret, nData = This->nSons-1;

  /* find position in node */
  if (nData < 4)
  {
    /* only a few entries -> use linear search */
    int found=FALSE;

    for(i=0; i<nData && !found; i++)
    {
      /* the func-ptr version is not used here. */
      /* cmp = (*cmp_func) (This->data[i], item);*/

      /* NOTE: the order of arguments for Compare is crucial!
                       first comes the existing key/item, 2nd the new one. */
      int cmp = CALL(CN(BTreeOf), Compare) (This->data[i], item);
      if (cmp==0)
      {
        /* already inserted, return and report FOUND */
        return(BTREE_FOUND);
      }

      if (cmp>0) found = TRUE;
    }
    if (found==TRUE) i--;
  }
  else
  {
    /* many entries -> use binary search */
    int l=0, r=nData-1;

    do {
      int m = (l+r)>>1;

      /* NOTE: the order of arguments for Compare is crucial!
                       first comes the existing key/item, 2nd the new one. */
      int cmp = CALL(CN(BTreeOf), Compare) (This->data[m], item);
      if (cmp==0)
      {
        /* already inserted, return and report FOUND */
        return(BTREE_FOUND);
      }

      if (cmp<0)
        l = m+1;
      else
        r = m-1;

    } while (l<=r);
    /* the correct position is between keys r and l (in that order!) */

    /* i will be set to the index at the right of the correct position */
    i = l;
  }


  if (This->sons[i]!=NULL)
  {
    /* case A: key must be in subtree i, recurse into it */

    ret = Method(Insert) (This->sons[i], /*cmp_func,*/ item);

    if (ret==BTREE_SPLIT_ME)
    {
      ClassPtr new_r;
      CN(BTreeOf) *split_item;
      Method(Split) (This->sons[i], &split_item, &new_r);

      if (i<nData)
      {
        Method(FreeIndex) (This, i);
        This->data[i] = split_item;
        This->sons[i+1] = new_r;
      }
      else
      {
        This->data[i] = split_item;
        This->sons[i+1] = new_r;
        This->nSons++;
      }

      if (This->nSons<BTreeOrder+1)
        ret = BTREE_INSERTED;
    }
  }
  else
  {
    /* case B: there is no subtree for key here, we insert it locally */

    if (i<nData)
    {
      Method(FreeIndex) (This, i);
      This->data[i] = item;
    }
    else
    {
      This->data[i] = item;
      This->sons[i+1] = NULL;
      This->nSons++;
    }
    ret = (This->nSons<BTreeOrder+1) ?
          BTREE_INSERTED : BTREE_SPLIT_ME;
  }

  return(ret);
}



static void Method(Iterate) (ParamThis, Iterate_Method iter_method)
{
  int i;

  /* find position in node */
  for(i=0; i<This->nSons-1; i++)
  {
    if (This->sons[i]!=NULL)
    {
      Method(Iterate) (This->sons[i], iter_method);
    }

    /* call iterator method */
    (*iter_method)(This->data[i]);
  }
  if (This->sons[i]!=NULL)
  {
    Method(Iterate) (This->sons[i], iter_method);
  }
}


static CN(BTreeOf) **Method(Linearize) (ParamThis, CN(BTreeOf) **ptr)
{
  int i;

  for(i=0; i<This->nSons-1; i++)
  {
    if (This->sons[i]!=NULL)
    {
      ptr = Method(Linearize) (This->sons[i], ptr);
    }

    *ptr++ = This->data[i];
  }
  if (This->sons[i]!=NULL)
  {
    ptr = Method(Linearize) (This->sons[i], ptr);
  }

  return(ptr);
}


static void Method(GetResources) (ParamThis,
                                  int *nNodes, size_t *alloc_mem, size_t *used_mem)
{
  size_t allocated=0, used=0;
  int i, nodes=0;

  for(i=0; i<This->nSons; i++)
  {
    if (This->sons[i]!=NULL)
    {
      int i1; size_t s1, s2;
      Method(GetResources) (This->sons[i], &i1, &s1, &s2);
      nodes     += i1;
      allocated += s1;
      used      += s2;
    }
  }

  *nNodes    = nodes     + 1;
  *alloc_mem = allocated + sizeof(CBTreeNode);
  *used_mem  = used      + ((1 + 2*This->nSons) * (sizeof(void *)));
}

#endif

#undef ClassName


/*** BTree class ***/

#define ClassName CBTree
Class_Data_Begin
CN(CBTreeNode) *root;
int nItems;

/*Compare_Method compare_func;*/
Class_Data_End
Method_New_      (/*Compare_Method*/ _NEWPARAMS_OR_VOID);
void         Method(Print)    (DefThis _PRINTPARAMS);
void         Method(Reset)    (DefThis);
int          Method(Insert)   (DefThis, CN(BTreeOf) *);
void         Method(Iterate)  (DefThis, Iterate_Method);
CPtrArray   *Method(GetArray) (DefThis);
void         Method(GetResources)  (DefThis, int *, int *, size_t *, size_t *);


#ifdef ContainerImplementation
Method_New_ (/*Compare_Method compare_func*/ _NEWPARAMS_OR_VOID)
{
  Construct(This, _CHECKALLOC(This));
  This->root   = NULL;
  This->nItems = 0;
  /*This->compare_func = compare_func;*/
  return(This);
}


void Method(Print) (ParamThis _PRINTPARAMS)
{
  CALL(CBTreeNode,Print) (This->root _PRINTSAME);
}


void Method(Reset) (ParamThis)
{
  if (This->root!=NULL)
  {
    CALL(CBTreeNode,Free) (This->root);
  }

  This->root   = NULL;
  This->nItems = 0;
}


int Method(Insert) (ParamThis, CN(BTreeOf) *item)
{
  int ret;

  if (This->root==NULL)
  {
    This->root = CALL(New,CBTreeNode) (item,NULL,NULL);
    This->nItems++;
    return(TRUE);
  }

  ret = CALL(CBTreeNode,Insert) (This->root, /*This->compare_func,*/ item);

  if (ret==BTREE_SPLIT_ME)
  {
    CBTreeNode *new_l, *new_r;
    CN(BTreeOf) *split_item;

    new_l = This->root;
    CALL(CBTreeNode,Split) (new_l, &split_item, &new_r);
    This->root = CALL(New,CBTreeNode) (split_item, new_l, new_r);
  }

  if (ret!=BTREE_FOUND)
    This->nItems++;

  /* return TRUE if item was really inserted */
  return(ret!=BTREE_FOUND);
}


void Method(Iterate) (ParamThis, Iterate_Method iter_method)
{
  if (This->root!=NULL)
    CALL(CBTreeNode,Iterate) (This->root, iter_method);
}


CPtrArray *Method(GetArray) (ParamThis)
{
  CBTreeNode *node;
  int i;
  CPtrArray  *array;
  CN(BTreeOf) **ptr;

  array = CALL(New,CPtrArray) (This->nItems);
  if (This->nItems==0)
    return(array);

  ptr = CALL(CPtrArray,GetData) (array);
  CALL(CBTreeNode,Linearize) (This->root, ptr);

  return(array);
}


void Method(GetResources) (ParamThis,
                           int *nNodes, int *nItems, size_t *alloc_mem, size_t *used_mem)
{
  size_t allocated=0, used=0;
  int nodes=0;

  if (This->root!=NULL)
  {
    CALL(CBTreeNode,GetResources) (This->root, &nodes, &allocated, &used);
  }

  *nNodes    = nodes;
  *nItems    = This->nItems;
  *alloc_mem = allocated + sizeof(CBTree);
  *used_mem  = used      + sizeof(CBTree);
}


#endif

#undef ClassName


/****************************************************************************/
#undef CPtrArray
#undef CBTreeNode
#undef CBTree
/*#undef Compare_Method*/
#undef BTreeOrder
#undef BTreeOf
#endif
/****************************************************************************/



/****************************************************************************/
#ifdef SetOf
#ifndef CSet
/****************************************************************************/

/*
        this class is implemented by using a SegmList for storing the actual
        data and a BTree for fast access, sorting and uniqueness.
 */

/* class name for set */
#define CSet       CCAT(SetOf,Set)


/* set default values for parameters */
#ifndef Set_SegmSize
#define Set_SegmSize 256        /* default size of each segment in SegmList */
#endif
#ifndef Set_BTreeOrder
#define Set_BTreeOrder 32       /* default order of BTree */
#endif



/* instantiate all template types we will need lateron */
#define PtrOf      SetOf
#define ArrayOf    CCAT(SetOf,Ptr)
#define SegmListOf SetOf
#define SegmSize   Set_SegmSize
#define BTreeOf    SetOf
#define BTreeOrder Set_BTreeOrder

#include "ooppcc.h"



/****************************************************************************/

#define Set(T)  CN(CCAT(T,Set))

/* class names for set, segmlist and btree */
#define CSet       CCAT(SetOf,Set)
#define CSegmList  CCAT(SetOf,SegmList)
#define CBTree     CCAT(SetOf,BTree)
#define CPtrArray  CCAT(CN(SetOf),PtrArray)


/*** Set class ***/

#define ClassName CSet
Class_Data_Begin
CSegmList     *list;
CBTree        *tree;

CN(SetOf)     *last_item;             /* temp storage for last call to NewItem() */
Class_Data_End
Method_New_      (_NEWPARAMS_OR_VOID);
void         Method(Print)    (DefThis _PRINTPARAMS);
void         Method(Reset)    (DefThis);
CN(SetOf)   *Method(NewItem)  (DefThis);
int          Method(ItemOK)   (DefThis);
int          Method(GetNItems)     (DefThis);
int          Method(GetNDiscarded) (DefThis);
CPtrArray   *Method(GetArray) (DefThis);
void    Method(GetResources)  (DefThis, int *, int *, int *, size_t *, size_t *);


#ifdef ContainerImplementation

Method_New_ (_NEWPARAMS_OR_VOID)
{
  Construct(This, _CHECKALLOC(This));
  This->list   = CALL(New, CSegmList) ();
  This->tree   = CALL(New, CBTree) ();
  This->last_item = NULL;
  return(This);
}


void Method(Print) (ParamThis _PRINTPARAMS)
{
  /* CALL(CSegmList,Print) (This->list  _PRINTNEXT); */
  CALL(CBTree,Print) (This->tree _PRINTNEXT);
}


void Method(Reset) (ParamThis)
{
  CALL(CSegmList,Reset) (This->list);
  CALL(CBTree,Reset) (This->tree);
}


CN(SetOf) *Method(NewItem) (ParamThis)
{
  This->last_item = CALL(CSegmList,NewItem) (This->list);
  return(This->last_item);
}


int Method(ItemOK) (ParamThis)
{
  if (CALL(CBTree,Insert) (This->tree, This->last_item))
  {
    /* item could be inserted in btree, hence it was a new one. */
    return(TRUE);
  }

  /* item was in tree already, it was an old one. */
  CALL(CSegmList,DiscardItem) (This->list);
  return(FALSE);
}


int Method(GetNItems) (ParamThis)
{
  return(CALL(CSegmList,GetNItems) (This->list));
}

int Method(GetNDiscarded) (ParamThis)
{
  return(CALL(CSegmList,GetNDiscarded) (This->list));
}


CPtrArray *Method(GetArray) (ParamThis)
{
  return(CALL(CBTree,GetArray) (This->tree));
}


void Method(GetResources) (ParamThis,
                           int *nSegms, int *nItems, int *nNodes, size_t *alloc_mem, size_t *used_mem)
{
  size_t allocated=0, used=0;

  CALL(CSegmList,GetResources) (This->list, nSegms, nItems, &allocated, &used);
  *alloc_mem = allocated;
  *used_mem  = used;

  CALL(CBTree,GetResources) (This->tree, nNodes, nItems, &allocated, &used);
  *alloc_mem += allocated;
  *used_mem  += used;

  *alloc_mem += sizeof(CSet);
  *used_mem  += sizeof(CSet);
}


#endif



#undef ClassName


/****************************************************************************/
#undef CSegmList
#undef CBTree
#undef CPtrArray
#undef CSet
#undef SetOf
#endif
#endif
/****************************************************************************/


/****************************************************************************/
