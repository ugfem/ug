// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  cw.c															*/
/*																			*/
/* Purpose:   define global array with predefined control word entries		*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*/
/*																			*/
/* History:   11.01.95 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <stdio.h>

/* define this to exclude extern definition of global arrays */
#define __COMPILE_CW__

#include "devices.h"
#include "debug.h"
#include "general.h"

#include "gm.h"
#include "misc.h"
#include "algebra.h"
#include "ugm.h"
#include "refine.h"
#include "cw.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define CW_EDOBJ                (BITWISE_TYPE(EDOBJ) | BITWISE_TYPE(LIOBJ))
/* take both edges and links in one	        */
#define CW_GROBJ                BITWISE_TYPE(GROBJ)
#define CW_MGOBJ                BITWISE_TYPE(MGOBJ)
#define CW_NDOBJ                BITWISE_TYPE(NDOBJ)
#define CW_VEOBJ                BITWISE_TYPE(VEOBJ)
#define CW_MAOBJ                (BITWISE_TYPE(MAOBJ) | BITWISE_TYPE(COOBJ))
/* take both matrix and connection in one	*/
#define CW_BVOBJ                BITWISE_TYPE(BLOCKVOBJ)

#define CW_VXOBJS               (BITWISE_TYPE(IVOBJ) | BITWISE_TYPE(BVOBJ))
#define CW_ELOBJS               (BITWISE_TYPE(IEOBJ) | BITWISE_TYPE(BEOBJ))
#define CW_GEOMOBJS             (CW_VXOBJS | CW_ELOBJS | CW_NDOBJ | CW_EDOBJ | CW_GROBJ)
/* NOTE: CW_GEOMOBJS and GEOM_OBJECTS differ*/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/* description of a control word predefines */
typedef struct {
  INT used;                                                             /* used this entry					*/
  char *name;                                                           /* name string						*/
  INT control_word_id;                                  /* index in control_words			*/
  unsigned INT offset_in_object ;               /* where in object is it ?			*/
  INT objt_used;                                                /* bitwise object ID				*/
} CONTROL_WORD_PREDEF;

/* description of a control enty predefines */
typedef struct {
  INT used;                                                             /* used this entry					*/
  char *name;                                                           /* name string						*/
  INT control_word;                                             /* index of corresp. controlword	*/
  INT control_entry_id;                                 /* index in control_entries             */
  INT offset_in_word;                                   /* shift in control word			*/
  INT length;                                                   /* number of bits used				*/
  INT objt_used;                                                /* bitwise object ID				*/
} CONTROL_ENTRY_PREDEF;

typedef struct {

  INT read;                                                             /* number of accesses to read		*/
  INT write;                                                            /* number of accesses to write		*/

} CE_USAGE;

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

CONTROL_WORD control_words[MAX_CONTROL_WORDS];
CONTROL_ENTRY control_entries[MAX_CONTROL_ENTRIES];

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static CONTROL_WORD_PREDEF cw_predefines[MAX_CONTROL_WORDS] = {
  CW_INIT(CW_USED,VECTOR_,                        CW_VEOBJ),
  CW_INIT(CW_USED,MATRIX_,                        CW_MAOBJ),
  CW_INIT(CW_USED,BLOCKVECTOR_,           CW_BVOBJ),
  CW_INIT(CW_USED,VERTEX_,                        CW_VXOBJS),
  CW_INIT(CW_USED,NODE_,                          CW_NDOBJ),
  CW_INIT(CW_USED,LINK_,                          CW_EDOBJ),
  CW_INIT(CW_USED,EDGE_,                          CW_EDOBJ),
  CW_INIT(CW_USED,ELEMENT_,                       CW_ELOBJS),
  CW_INIT(CW_USED,FLAG_,                          CW_ELOBJS),
  CW_INIT(CW_USED,PROPERTY_,                      CW_ELOBJS),
  CW_INIT(CW_USED,GRID_,                          CW_GROBJ),
  CW_INIT(CW_USED,GRID_STATUS_,           CW_GROBJ),
  CW_INIT(CW_USED,MULTIGRID_STATUS_,      CW_MGOBJ),
  CW_INIT_UNUSED,
  CW_INIT_UNUSED,
  CW_INIT_UNUSED,
  CW_INIT_UNUSED,
  CW_INIT_UNUSED,
  CW_INIT_UNUSED,
  CW_INIT_UNUSED,
};

static CONTROL_ENTRY_PREDEF ce_predefines[MAX_CONTROL_ENTRIES] = {
  CE_INIT(CE_LOCKED,      VECTOR_,                VOTYPE_,                CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VCOUNT_,                CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VECTORSIDE_,    CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VCLASS_,                CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VDATATYPE_,             CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VNCLASS_,               CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VNEW_,                  CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VCCUT_,                 CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VTYPE_,                 CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VPART_,                 CW_VEOBJ),
  CE_INIT(CE_LOCKED,      VECTOR_,                VCCOARSE_,              CW_VEOBJ),

  CE_INIT(CE_LOCKED,      MATRIX_,                MOFFSET_,               CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                MROOTTYPE_,             CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                MDESTTYPE_,             CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                MDIAG_,                 CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                MSIZE_,                 CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                MNEW_,                  CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                CEXTRA_,                CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                MDOWN_,                 CW_MAOBJ),
  CE_INIT(CE_LOCKED,      MATRIX_,                MUP_,                   CW_MAOBJ),

  CE_INIT(CE_USED,        BLOCKVECTOR_,   BVDOWNTYPE_,    CW_BVOBJ),
  CE_INIT(CE_USED,        BLOCKVECTOR_,   BVLEVEL_,               CW_BVOBJ),
  CE_INIT(CE_USED,        BLOCKVECTOR_,   BVTVTYPE_,              CW_BVOBJ),

  CE_INIT(CE_LOCKED,      GENERAL_,               OBJ_,                   (CW_GEOMOBJS | CW_VEOBJ | CW_MAOBJ)),
  CE_INIT(CE_LOCKED,      GENERAL_,               USED_,                  (CW_GEOMOBJS | CW_VEOBJ | CW_MAOBJ)),
  CE_INIT(CE_LOCKED,      GENERAL_,               THEFLAG_,               (CW_GEOMOBJS | CW_VEOBJ)),
  CE_INIT(CE_LOCKED,      GENERAL_,               LEVEL_,                 CW_GEOMOBJS),

  CE_INIT(CE_LOCKED,      VERTEX_,                MOVE_,                  CW_VXOBJS),
  CE_INIT(CE_LOCKED,      VERTEX_,                MOVED_,                 CW_VXOBJS),
  CE_INIT(CE_LOCKED,      VERTEX_,                ONEDGE_,                CW_VXOBJS),
  CE_INIT(CE_LOCKED,      VERTEX_,                ONSIDE_,                CW_VXOBJS),
  CE_INIT(CE_LOCKED,      VERTEX_,                ONNBSIDE_,              CW_VXOBJS),
  CE_INIT(CE_LOCKED,      VERTEX_,                NOOFNODE_,              CW_VXOBJS),

  CE_INIT(CE_LOCKED,      NODE_,                  NSUBDOM_,               CW_NDOBJ),
  CE_INIT(CE_LOCKED,      NODE_,                  NPROP_,                 CW_NDOBJ),
  CE_INIT(CE_LOCKED,      NODE_,                  MODIFIED_,              (CW_NDOBJ | CW_GROBJ)),
  CE_INIT(CE_LOCKED,      NODE_,                  NTYPE_,                 CW_NDOBJ),

  CE_INIT(CE_USED,        LINK_,                  LOFFSET_,               CW_EDOBJ),

  CE_INIT(CE_USED,        EDGE_,                  AUXEDGE_,               CW_EDOBJ),
  CE_INIT(CE_USED,        EDGE_,                  PATTERN_,               CW_EDOBJ),
  CE_INIT(CE_USED,        EDGE_,                  ADDPATTERN_,    CW_EDOBJ),
  CE_INIT(CE_USED,        EDGE_,                  EDGENEW_,               CW_EDOBJ),
  CE_INIT(CE_USED,        EDGE_,                  EDSUBDOM_,              CW_EDOBJ),
  CE_INIT(CE_USED,        EDGE_,                  NO_OF_ELEM_,    CW_EDOBJ),

  CE_INIT(CE_USED,        ELEMENT_,               REFINE_,                CW_ELOBJS),
  CE_INIT(CE_USED,        ELEMENT_,               ECLASS_,                CW_ELOBJS),
  CE_INIT(CE_USED,        ELEMENT_,               NSONS_,                 CW_ELOBJS),
  CE_INIT(CE_USED,        ELEMENT_,               REFINECLASS_,   CW_ELOBJS),
  CE_INIT(CE_USED,        ELEMENT_,               NEWEL_,                 CW_ELOBJS),
  CE_INIT(CE_USED,        ELEMENT_,               TAG_,                   CW_ELOBJS),

  CE_INIT(CE_USED,        FLAG_,                  MARK_,                  CW_ELOBJS),
  CE_INIT(CE_USED,        FLAG_,                  COARSEN_,               CW_ELOBJS),
  CE_INIT(CE_USED,        FLAG_,                  EBUILDCON_,             CW_ELOBJS),
  CE_INIT(CE_USED,        FLAG_,                  DECOUPLED_,             CW_ELOBJS),
  CE_INIT(CE_USED,        FLAG_,                  UPDATE_GREEN_,  CW_ELOBJS),
  CE_INIT(CE_USED,        FLAG_,                  SIDEPATTERN_,   CW_ELOBJS),
  CE_INIT(CE_USED,        FLAG_,                  MARKCLASS_,             CW_ELOBJS),

  CE_INIT(CE_USED,        PROPERTY_,              SUBDOMAIN_,             CW_ELOBJS),
  CE_INIT(CE_USED,        PROPERTY_,              NODEORD_,               CW_ELOBJS),
  CE_INIT(CE_USED,        PROPERTY_,              PROP_,                  CW_ELOBJS),

        #ifdef ModelP
  CE_INIT(CE_USED,        VECTOR_,                XFERVECTOR_,    CW_VEOBJ),
  CE_INIT(CE_USED,        MATRIX_,                XFERMATX_,              CW_MAOBJ),
        #else /* ModelP */
  CE_INIT_UNUSED,
  CE_INIT_UNUSED,
  CE_INIT_UNUSED,
  CE_INIT_UNUSED,
        #endif /* ModelP */
  CE_INIT_UNUSED,
  CE_INIT_UNUSED,
  CE_INIT_UNUSED,
  CE_INIT_UNUSED,
};

static CE_USAGE ce_usage[MAX_CONTROL_ENTRIES];

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*D
   INT_2_bitpattern	- transform an INT into a bitpattern string

   SYNOPSIS:
   static void INT_2_bitpattern (INT n, char *text)

   PARAMETERS:
   .  n - integer to convert
   .  text - string of size >= 33 for conversion

   DESCRIPTION:
   This function transforms an INT into a bitpattern string consisting of 0s
   and 1s only.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

static void INT_2_bitpattern (INT n, char *text)
{
  INT i;

  memset(text,'0',32*sizeof(char));

  for (i=0; i<32; i++)
    if ((1<<i)&n)
      text[31-i] = '1';
  text[32] = '\0';

  return;
}

/****************************************************************************/
/*D
   ListCWofObject	- print all control entries of an objects control word

   SYNOPSIS:
   void ListCWofObject (const void *obj, INT offset)

   PARAMETERS:
   .  obj - object pointer
   .  offset - controlword offset in (unsigned INT) in object

   DESCRIPTION:
   This function prints the contents of all control entries of an objects control word at a
   given offset.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListCWofObject (const void *obj, INT offset)
{
  INT i,n,ce,last_ce,sub,min,cw_objt,oiw;

  ASSERT(obj!=NULL);
  HEAPFAULT(obj);

  cw_objt = BITWISE_TYPE(OBJT(obj));
  sub = -1;
  last_ce = -1;

  /* print control word entries in ascending order of offsets in word */
  do
  {
    min = MAX_I;
    for (i=0; i<MAX_CONTROL_ENTRIES; i++)
      if (control_entries[i].used)
        if (control_entries[i].objt_used & cw_objt)
          if (control_entries[i].offset_in_object==offset)
          {
            oiw = control_entries[i].offset_in_word;
            if ((oiw<min) && (oiw>=sub))
            {
              if ((oiw==sub) && (i<=last_ce))
                continue;
              ce = i;
              min = oiw;
            }
          }
    if (min==MAX_I)
      break;

    n = CW_READ(obj,ce);
    UserWriteF("  ce %s with offset in cw %3d: %10d\n",control_entries[i].name,min,n);
    sub = min;
    last_ce = ce;
  }
  while (TRUE);

  ASSERT(sub>=0);
}

/****************************************************************************/
/*D
   ListAllCWsOfObject	- print all control entries of all control words of an object

   SYNOPSIS:
   void ListAllCWsOfObject (const void *obj)

   PARAMETERS:
   .  obj - object pointer

   DESCRIPTION:
   This function prints the contents of all control entries of all control words
   of the object. 'ListCWofObject' is called.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ListAllCWsOfObject (const void *obj)
{
  INT i,cw,last_cw,sub,min,cw_objt,offset;

  ASSERT(obj!=NULL);
  HEAPFAULT(obj);

  cw_objt = BITWISE_TYPE(OBJT(obj));
  sub = -1;
  last_cw = -1;

  /* print control word contents in ascending order of offsets */
  do
  {
    min = MAX_I;
    for (i=0; i<MAX_CONTROL_WORDS; i++)
      if (control_words[i].used)
        if (control_words[i].objt_used & cw_objt)
        {
          offset = control_words[i].offset_in_object;
          if ((offset<min) && (offset>=sub))
          {
            if ((offset==sub) && (i<=last_cw))
              continue;
            cw = i;
            min = offset;
          }
        }
    if (min==MAX_I)
      break;

    UserWriteF("cw %s with offset %3d:\n",control_words[cw].name,min);
    ListCWofObject(obj,min);
    sub = min;
    last_cw = cw;
  }
  while (TRUE);

  ASSERT(sub>=0);
}

/****************************************************************************/
/*D
   ListCWofObjectType	- print used pattern of all control entries of an object types control word

   SYNOPSIS:
   static void ListCWofObjectType (INT objt, INT offset)

   PARAMETERS:
   .  obj - object pointer
   .  offset - controlword offset in (unsigned INT) in object

   DESCRIPTION:
   This function prints the used pattern of all control entries of an object types control word at a
   given offset.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

static void ListCWofObjectType (INT objt, INT offset)
{
  INT i,ce,last_ce,sub,min,cw_objt,oiw;
  char bitpat[33];

  cw_objt = BITWISE_TYPE(objt);
  sub = -1;
  last_ce = -1;

  /* print control word entries in ascending order of offsets in word */
  do
  {
    min = MAX_I;
    for (i=0; i<MAX_CONTROL_ENTRIES; i++)
      if (control_entries[i].used)
        if (control_entries[i].objt_used & cw_objt)
          if (control_entries[i].offset_in_object==offset)
          {
            oiw = control_entries[i].offset_in_word;
            if ((oiw<min) && (oiw>=sub))
            {
              if ((oiw==sub) && (i<=last_ce))
                continue;
              ce = i;
              min = oiw;
            }
          }
    if (min==MAX_I)
      break;

    INT_2_bitpattern(control_entries[ce].mask,bitpat);
    printf("  ce %-20s offset in cw %3d, len %3d: %s\n",
           control_entries[ce].name,
           control_entries[ce].offset_in_word,
           control_entries[ce].length,
           bitpat);
    sub = min;
    last_ce = ce;
  }
  while (TRUE);

  if (sub==-1)
    printf(" --- no ce found with objt %d\n",objt);
}

/****************************************************************************/
/*D
   ListAllCWsOfObjectType	- print used pattern of all control entries of all
                                                                control words of an object type

   SYNOPSIS:
   static void ListAllCWsOfObjectType (INT objt)

   PARAMETERS:
   .  obj - object pointer

   DESCRIPTION:
   This function prints the used pattern of all control entries of all control words
   of an object type. 'ListCWofObjectType' is called.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

static void ListAllCWsOfObjectType (INT objt)
{
  INT i,cw,last_cw,sub,min,cw_objt,offset;

  cw_objt = BITWISE_TYPE(objt);
  sub = -1;
  last_cw = -1;

  /* print control word contents in ascending order of offsets */
  do
  {
    min = MAX_I;
    for (i=0; i<MAX_CONTROL_WORDS; i++)
      if (control_words[i].used)
        if (control_words[i].objt_used & cw_objt)
        {
          offset = control_words[i].offset_in_object;
          if ((offset<min) && (offset>=sub))
          {
            if ((offset==sub) && (i<=last_cw))
              continue;
            cw = i;
            min = offset;
          }
        }
    if (min==MAX_I)
      break;

    printf("cw %-20s with offset in object %3d (unsigned INTs):\n",control_words[cw].name,min);
    ListCWofObjectType(objt,min);
    sub = min;
    last_cw = cw;
  }
  while (TRUE);

  if (sub==-1)
    printf(" --- no cw found with objt %d\n",objt);
}

/****************************************************************************/
/*D
   InitPredefinedControlWords	- Initialize control words

   SYNOPSIS:
   INT InitPredefinedControlWords (void)

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes the predefined control words.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if error occured.
   D*/
/****************************************************************************/

static INT InitPredefinedControlWords (void)
{
  INT i,nused;
  CONTROL_WORD *cw;
  CONTROL_WORD_PREDEF *pcw;

  /* clear everything */
  memset(control_words,0,MAX_CONTROL_WORDS*sizeof(CONTROL_WORD));

  nused = 0;
  for (i=0; i<MAX_CONTROL_WORDS; i++)
    if (cw_predefines[i].used)
    {
      pcw = cw_predefines+i;
      ASSERT(pcw->control_word_id<MAX_CONTROL_WORDS);

      nused++;
      cw = control_words+pcw->control_word_id;
      if (cw->used)
      {
        printf("redefinition of control word '%s'\n",pcw->name);
        return(__LINE__);
      }
      cw->used = pcw->used;
      cw->name = pcw->name;
      cw->offset_in_object = pcw->offset_in_object;
      cw->objt_used = pcw->objt_used;
    }

  if (nused!=GM_N_CW)
  {
    printf("InitPredefinedControlWords: nused=%d != GM_N_CW=%d\n",nused,GM_N_CW);
    assert(FALSE);
  }

  return (GM_OK);
}

/****************************************************************************/
/*D
   InitPredefinedControlEntries	- Initialize control word entries

   SYNOPSIS:
   INT InitPredefinedControlEntries (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes the predefined control word entries. Predefined
   entries are not checked for overlap.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if error occured.
   D*/
/****************************************************************************/

static INT InitPredefinedControlEntries (void)
{
  CONTROL_ENTRY *ce,*test_ce;
  CONTROL_WORD *cw;
  CONTROL_ENTRY_PREDEF *pce,*test_pce;
  INT i,j,k,offset,mask,error,nused;

  /* clear everything */
  memset(control_entries,0,MAX_CONTROL_ENTRIES*sizeof(CONTROL_ENTRY));

  error = 0;
  nused = 0;
  for (i=0; i<MAX_CONTROL_ENTRIES; i++)
    if (ce_predefines[i].used)
    {
      pce = ce_predefines+i;
      ASSERT(pce->control_entry_id<MAX_CONTROL_ENTRIES);

      nused++;
      ce = control_entries+pce->control_entry_id;
      if (ce->used)
      {
        printf("redefinition of control entry '%s'\n",pce->name);
        return(__LINE__);
      }
      cw = control_words+pce->control_word;
      ASSERT(cw->used);
      ce->used = pce->used;
      ce->name = pce->name;
      ce->control_word = pce->control_word;
      ce->offset_in_word = pce->offset_in_word;
      ce->length = pce->length;
      ce->objt_used = pce->objt_used;
      ce->offset_in_object = cw->offset_in_object;
      ce->mask = (POW2(ce->length)-1)<<ce->offset_in_word;
      ce->xor_mask = ~ce->mask;

      PRINTDEBUG(gm,1,("ceID %d used %d name %s cw %d\n",
                       i,ce->used,ce->name,ce->control_word));

      ASSERT(ce->objt_used & cw->objt_used);                                    /* ce and cw have! common objects */

      /* set mask in all cws that use some of the ces objects and have the same offset than cw */
      offset = ce->offset_in_object;
      mask   = ce->mask;
      for (k=0; k<MAX_CONTROL_WORDS; k++)
      {
        cw = control_words+k;

        if (!cw->used)
          continue;
        if (!(ce->objt_used & cw->objt_used))
          continue;
        if (cw->offset_in_object!=offset)
          continue;

        /* do other control entries overlap? */
        if (cw->used_mask & mask)
        {
          IFDEBUG(gm,1)
          printf("predef ctrl entry '%s' has overlapping bits with previous ctrl entries:\n",pce->name);
          for (j=0; j<i; j++)
          {
            test_pce = ce_predefines+j;
            test_ce  = control_entries+test_pce->control_entry_id;
            if (test_ce->objt_used & ce->objt_used)
              if (test_ce->offset_in_object==offset)
                if (test_ce->mask & mask)
                  printf(" '%s'",test_pce->name);

          }
          printf("\n");
          ENDDEBUG
            error++;
        }
        cw->used_mask |= mask;
      }
    }

  IFDEBUG(gm,1)
  if (me == master) {
    ListAllCWsOfObjectType(IVOBJ);
    ListAllCWsOfObjectType(IEOBJ);
    ListAllCWsOfObjectType(EDOBJ);
    ListAllCWsOfObjectType(NDOBJ);
    ListAllCWsOfObjectType(VEOBJ);
    ListAllCWsOfObjectType(MAOBJ);
    ListAllCWsOfObjectType(BLOCKVOBJ);
    ListAllCWsOfObjectType(GROBJ);
    ListAllCWsOfObjectType(MGOBJ);
  }
  ENDDEBUG

  /* TODO: enable next lines for error control */
  IFDEBUG(gm,1)
  if (error)
    return (__LINE__);
  ENDDEBUG

  if (nused!=REFINE_N_CE)
  {
    printf("InitPredefinedControlEntries: nused=%d != REFINE_N_CE=%d\n",nused,REFINE_N_CE);
    assert(FALSE);
  }

  return (GM_OK);
}

/****************************************************************************/
/*D
   ResetCEstatistics	- reset counters for read/write control word access

   SYNOPSIS:
   void ResetCEstatistics (void)

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function resets all counters for read/write control word access.
   This is only possible if the code is compiled with #define _DEBUG_CW_
   in gm.h. (Maybe it makes sense to do that only for part of the code.
   Then the warnings should be removed.)

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void ResetCEstatistics (void)
{
        #ifndef _DEBUG_CW_
  PrintErrorMessage('W',"ResetCEstatistics","compile with #ifdef _DEBUG_CW_ in gm.h!");
        #else
  memset(ce_usage,0,MAX_CONTROL_ENTRIES*sizeof(CE_USAGE));
        #endif
}

/****************************************************************************/
/*D
   PrintCEstatistics	- print control word read/write acces statistic to shell

   SYNOPSIS:
   void PrintCEstatistics (void)

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function prints all counters for read/write control word access
   (only if at least one of them is nonzero).
   This is only possible if the code is compiled with #define _DEBUG_CW_
   in gm.h. (Maybe it makes sense to do that only for part of the code.
   Then the warnings should be removed.)

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void PrintCEstatistics (void)
{
  INT i;

        #ifndef _DEBUG_CW_
  PrintErrorMessage('W',"PrintCEstatistics","compile with #ifdef _DEBUG_CW_ in gm.h!");
        #else
  for (i=0; i<MAX_CONTROL_ENTRIES; i++)
    if (control_entries[i].used)
      if (ce_usage[i].read || ce_usage[i].write)
        if (control_entries[i].name!=NULL)
          UserWriteF("ce %-20s: read %10d write %10d\n",control_entries[i].name,ce_usage[i].read,ce_usage[i].write);
        else
          UserWriteF("ce %20d: read %10d write %10d\n",i,ce_usage[i].read,ce_usage[i].write);
        #endif
}

/****************************************************************************/
/*D
   ReadCW	- function to replace CW_READ macro and does extended error checks

   SYNOPSIS:
   unsigned INT ReadCW (const void *obj, INT ceID)

   PARAMETERS:
   .  obj - object pointer
   .  ceID - control entry ID

   DESCRIPTION:
   This function is to replace the CW_READ and CW_READ_STATIC macros of gm.h and does extended
   error checks:~
   .n   obj != NULL
   .n   HEAPFAULT
   .n   ceID in valid range
   .n   control entry used
   .n   if ce uses a geom object the object type is checked
   Additionally the read accesses to a control entry are counted.

   CAUTION:
   set #define _DEBUG_CW_ to replace CW_READ by ReadCW but be aware of the
   terribble slowing down of the program in this case (no large problems!).

   RETURN VALUE:
   unsigned INT
   .n   number read from the control entry of the object
   D*/
/****************************************************************************/

unsigned INT ReadCW (const void *obj, INT ceID)
{
  CONTROL_ENTRY *ce;
  unsigned INT off_in_obj,mask,i,off_in_wrd,cw,cw_objt;

  ASSERT(obj!=NULL);
  HEAPFAULT(obj);
  if ((ceID<0) || (ceID>=MAX_CONTROL_ENTRIES))
  {
    printf("ReadCW: ceID=%d out of range\n",ceID);
    assert(FALSE);
  }

  ce_usage[ceID].read++;

  ce = control_entries+ceID;

  if (!ce->used)
  {
    printf("ReadCW: ceID=%d unused\n",ceID);
    assert(FALSE);
  }

  cw_objt = BITWISE_TYPE(OBJT(obj));
  if (!(cw_objt & ce->objt_used))
  {
    printf("ReadCW: invalid objt %d for ce %s\n",OBJT(obj),ce->name);
    assert(FALSE);
  }

  off_in_wrd = ce->offset_in_word;
  off_in_obj = ce->offset_in_object;
  mask = ce->mask;
  cw = ((unsigned INT *)(obj))[off_in_obj];
  i = cw & mask;
  i = i>>off_in_wrd;

  return (i);
}

/****************************************************************************/
/*D
   WriteCW	- function to replace CW_WRITE macro and does extended error checks

   SYNOPSIS:
   void WriteCW (void *obj, INT ceID, INT n)

   PARAMETERS:
   .  obj - object pointer
   .  ceID - control entry ID
   .  n - number to write to the objects control entry

   DESCRIPTION:
   This function is to replace the CW_WRITE and CW_WRITE_STATIC macros of gm.h and does extended
   error checks:~
   .n   obj != NULL
   .n   HEAPFAULT
   .n   ceID in valid range
   .n   control entry used
   .n   if ce uses a geom object the object type is checked
   .n   n small enough for length of control entry
   Additionally the write accesses to a control entry are counted.

   CAUTION:
   set #define _DEBUG_CW_ to replace CW_WRITE by WriteCW but be aware of the
   terribble slowing down of the program in this case (no large problems!).

   RETURN VALUE:
   unsigned INT
   .n   number read from the control entry of the object
   D*/
/****************************************************************************/

void WriteCW (void *obj, INT ceID, INT n)
{
  CONTROL_ENTRY *ce;
  unsigned INT off_in_obj,mask,i,j,off_in_wrd,cw_objt,xmsk;
  unsigned INT *pcw;

  ASSERT(obj!=NULL);
  HEAPFAULT(obj);
  if ((ceID<0) || (ceID>=MAX_CONTROL_ENTRIES))
  {
    printf("WriteCW: ceID=%d out of range\n",ceID);
    assert(FALSE);
  }

  ce_usage[ceID].write++;

  ce = control_entries+ceID;

  if (!ce->used)
  {
    printf("WriteCW: ceID=%d unused\n",ceID);
    assert(FALSE);
  }

  cw_objt = BITWISE_TYPE(OBJT(obj));

  /* special: SETOBJT cannot be checked */
  if (cw_objt==BITWISE_TYPE(0)) {
    if (ceID!=OBJ_CE)
      if (cw_objt != ce->objt_used) {
        printf("WriteCW: objt 0 but no SETOBJT access\n");
        assert(FALSE);
      }
  }
  else if (!(cw_objt & ce->objt_used))
  {
    printf("WriteCW: invalid objt %d for ce %s\n",OBJT(obj),ce->name);
    assert(FALSE);
  }

  off_in_wrd = ce->offset_in_word;
  off_in_obj = ce->offset_in_object;
  mask = ce->mask;
  xmsk = ce->xor_mask;
  pcw = ((unsigned INT *)(obj)) + off_in_obj;
  i = (*pcw) & xmsk;
  j = n<<off_in_wrd;

  if (j>mask)
  {
    printf("WriteCW: value=%d exceeds max=%d for %d\n",
           n,POW2(ce->length)-1,ce->name);
    assert(FALSE);
  }

  j = j & mask;

  *pcw = i | j;
}

/****************************************************************************/
/*D
   AllocateControlEntry	-  Allocates space in object control words

   SYNOPSIS:
   INT AllocateControlEntry (INT cw_id, INT length, INT *ce_id);

   PARAMETERS:
   .  cw_id - id of a control word
   .  length - number of bits to allocate
   .  ce_id -  returns identifier of control entry descriptor

   DESCRIPTION:
   This function allocates 'length' consecutive bits in the control word of an
   object identified through the `control word id` 'cw_id'.
   It returns '0' and a valid id in 'ce_id' if space was available.
   The 'ce_id' can be used to read and write the requested bits with the
   'CW_READ' and 'CW_WRITE' macros as in the following example.

   The control word ids of all UG objects are defined in 'gm.h' and usually
   have the name `<object>`'_CW', except e.g. the 'ELEMENT' which has three
   words that are used bitwise.

   EXAMPLE:

   The following code fragment allocates 'NORDER_LEN' bits in the 'flag' word
   of the 'ELEMENT' structure.

   .vb
   INT ce_NORDER;

   if (AllocateControlEntry(FLAG_CW,NORDER_LEN,&ce_NORDER) != GM_OK) return (1);
   .ve

   The following macros then read and write the requested bits

   .vb
   #define NORDER(p)      CW_READ(p,ce_NORDER)
   #define SETNORDER(p,n) CW_WRITE(p,ce_NORDER,n)
   .ve


   RETURN VALUE:
   INT
   .n    GM_OK if ok
   .n    GM_ERROR if error occured.
   D*/
/****************************************************************************/

INT AllocateControlEntry (INT cw_id, INT length, INT *ce_id)
{
  INT free, i, offset;
  CONTROL_ENTRY *ce;
  CONTROL_WORD *cw;
  unsigned INT mask;

  /* check input */
  if ((length<0)||(length>=32)) return(GM_ERROR);
  if ((cw_id<0)||(cw_id>=MAX_CONTROL_WORDS)) return(GM_ERROR);

  /* it is sufficient to check only the control entries control word
     multiple object types are only allowed for predefines */
  cw = control_words+cw_id;

  /* find unused entry */
  for (i=0; i<MAX_CONTROL_ENTRIES; i++)
    if (!control_entries[i].used) break;
  if (i==MAX_CONTROL_ENTRIES) return(GM_ERROR);
  free = i;
  ce = control_entries+free;

  /* lets see if enough consecutive bits are available */
  mask = POW2(length)-1;
  for (i=0; i<=32-length; i++)
  {
    if ((mask&cw->used_mask)==0) break;
    mask <<= 1;
  }
  if (i>32-length) return(GM_ERROR);
  offset = i;

  /* fill new entry */
  *ce_id = free;
  ce->used = 1;
  ce->control_word = cw_id;
  ce->offset_in_object = cw->offset_in_object;
  ce->offset_in_word = offset;
  ce->length = length;
  ce->mask = mask;
  ce->xor_mask = ~mask;
  ce->name = NULL;
  ce->objt_used = cw->objt_used;

  /* remember used bits */
  cw->used_mask |= mask;

  /* ok, exit */
  return(GM_OK);
}

/****************************************************************************/
/*D
   FreeControlEntry - Frees space in object control words

   SYNOPSIS:
   INT FreeControlEntry (INT ce_id);

   PARAMETERS:
   .  ce_id - control entry descriptor to free

   DESCRIPTION:
   This function frees space in object control words that has been allocated
   with 'AllocateControlEntry'.

   RETURN VALUE:
   INT
   .n     GM_OK if ok
   .n     GM_ERROR if error occured.
   D*/
/****************************************************************************/

INT FreeControlEntry (INT ce_id)
{
  CONTROL_ENTRY *ce;
  CONTROL_WORD *cw;

  /* check parameter */
  if ((ce_id<0)||(ce_id>=MAX_CONTROL_ENTRIES)) return(GM_ERROR);
  ce = control_entries+ce_id;
  cw = control_words+ce->control_word;

  /* check if locked */
  if (ce->used == 2)
    return(GM_ERROR);

  /* free used bits */
  cw->used_mask &= ce->xor_mask;

  /* free control entry */
  ce->used = 0;

  /* ok, exit */
  return(GM_OK);
}

INT PrintCW (void)
{
  INT i;
  CONTROL_ENTRY *ce;

  for (i=0; i<MAX_CONTROL_ENTRIES; i++)
    if (control_entries[i].used)
    {
      ce = control_entries+i;
      UserWriteF("ce %-20s in cw %-20s offset %3d length %3d\n",
                 ce->name,control_words[ce->control_word].name,ce->offset_in_object,ce->length);
    }
  return(GM_OK);
}


/****************************************************************************/
/*D
   InitCW - init cw.c file

   SYNOPSIS:
   INT InitCW (void)

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes the control word manager.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   > 0 line in which error occured.
   D*/
/****************************************************************************/

INT InitCW (void)
{
  if (InitPredefinedControlWords())
    return (__LINE__);
  if (InitPredefinedControlEntries())
    return (__LINE__);

        #ifdef _DEBUG_CW_
  ResetCEstatistics();
        #endif

  return (GM_OK);
}
