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

/* define this to exclude extern definition of global arrays */
#define __COMPILE_CW__

#include "devices.h"
#include "debug.h"
#include "general.h"

#include "switch.h"
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

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {
  INT used;                                                             /* used this entry					*/
  INT control_entry_id;                                 /* index in control_entries             */
  INT control_word ;                                            /* pointer to corresponding controlw*/
  INT offset_in_word;                                   /* shift in control word			*/
  INT length;                                                   /* number of bits used				*/
} predefined_control_entry;

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

CONTROL_WORD control_words[MAX_CONTROL_WORDS] = {
  {VECTOR_OFFSET, 0},
  {MATRIX_OFFSET, 0},
  {GENERAL_OFFSET, 0},
  {VERTEX_OFFSET, 0},
  {NODE_OFFSET, 0},
  {LINK_OFFSET, 0},
  {EDGE_OFFSET, 0},
  {ELEMENT_OFFSET, 0},
  {FLAG_OFFSET, 0},
  {GRID_CW_OFFSET, 0},
  {GRID_STATUS_OFFSET, 0},
  {MULTIGRID_STATUS_OFFSET, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}
} ;

CONTROL_ENTRY control_entries[MAX_CONTROL_ENTRIES];

predefined_control_entry predefines[MAX_CONTROL_ENTRIES] = {
  {1,VTYPE_CE,            VECTOR_CW,      VTYPE_SHIFT,            VTYPE_LEN               },
  {1,VBUILDCON_CE,        VECTOR_CW,      VBUILDCON_SHIFT,        VBUILDCON_LEN   },
  {1,VCUSED_CE,           VECTOR_CW,      VCUSED_SHIFT,           VCUSED_LEN              },
  {1,VCOUNT_CE,           VECTOR_CW,      VCOUNT_SHIFT,           VCOUNT_LEN              },
  {1,VECTORSIDE_CE,       VECTOR_CW,      VECTORSIDE_SHIFT,       VECTORSIDE_LEN  },
  {1,VCLASS_CE,           VECTOR_CW,      VCLASS_SHIFT,           VCLASS_LEN      },
  {1,VDATATYPE_CE,        VECTOR_CW,      VDATATYPE_SHIFT,        VDATATYPE_LEN   },
  {1,VNCLASS_CE,          VECTOR_CW,      VNCLASS_SHIFT,          VNCLASS_LEN     },
  {1,VNEW_CE,             VECTOR_CW,      VNEW_SHIFT,             VNEW_LEN        },
  {1,VCNEW_CE,            VECTOR_CW,      VCNEW_SHIFT,            VCNEW_LEN       },
  {1,VCNB_CE,                     VECTOR_CW,      VCNB_SHIFT,                     VCNB_LEN                },
  {1,VCCUT_CE,            VECTOR_CW,      VCCUT_SHIFT,            VCCUT_LEN               },

  {1,MOFFSET_CE,          MATRIX_CW,      MOFFSET_SHIFT,          MOFFSET_LEN             },
  {1,MROOTTYPE_CE,        MATRIX_CW,      MROOTTYPE_SHIFT,        MROOTTYPE_LEN   },
  {1,MDESTTYPE_CE,        MATRIX_CW,      MDESTTYPE_SHIFT,        MDESTTYPE_LEN   },
  {1,MDIAG_CE,            MATRIX_CW,      MDIAG_SHIFT,            MDIAG_LEN               },
  {1,MTYPE_CE,            MATRIX_CW,      MTYPE_SHIFT,            MTYPE_LEN               },
  {1,MUSED_CE,            MATRIX_CW,      MUSED_SHIFT,            MUSED_LEN               },
  {1,MSIZE_CE,            MATRIX_CW,      MSIZE_SHIFT,            MSIZE_LEN               },
  {1,MNEW_CE,             MATRIX_CW,      MNEW_SHIFT,             MNEW_LEN                },
  {1,CEXTRA_CE,           MATRIX_CW,      CEXTRA_SHIFT,           CEXTRA_LEN              },
  {1,MDOWN_CE,            MATRIX_CW,      MDOWN_SHIFT,            MDOWN_LEN               },
  {1,MUP_CE,              MATRIX_CW,      MUP_SHIFT,              MUP_LEN                 },

  {1,BVDOWNTYPE_CE,       BLOCKVECTOR_CW, BVDOWNTYPE_SHIFT,       BVDOWNTYPE_LEN  },

  {1,OBJ_CE,              GENERAL_CW,     OBJ_SHIFT,              OBJ_LEN                 },
  {1,USED_CE,             GENERAL_CW,     USED_SHIFT,             USED_LEN                },
  {1,TAG_CE,              GENERAL_CW,     TAG_SHIFT,              TAG_LEN                 },
  {1,LEVEL_CE,            GENERAL_CW,     LEVEL_SHIFT,            LEVEL_LEN               },
  {1,THEFLAG_CE,          GENERAL_CW,     THEFLAG_SHIFT,          THEFLAG_LEN             },

  {1,VERTEX_GEN,          VERTEX_CW,      GENERAL_SHIFT,          GENERAL_LEN             },
  {1,MOVE_CE,             VERTEX_CW,      MOVE_SHIFT,             MOVE_LEN                },
  {1,MOVED_CE,            VERTEX_CW,      MOVED_SHIFT,            MOVED_LEN               },
  {1,ONEDGE_CE,           VERTEX_CW,      ONEDGE_SHIFT,           ONEDGE_LEN              },

  {1,NODE_GEN,            NODE_CW,        GENERAL_SHIFT,          GENERAL_LEN             },
  {1,CLASS_CE,            NODE_CW,        CLASS_SHIFT,            CLASS_LEN               },
  {1,NPROP_CE,            NODE_CW,        NPROP_SHIFT,            NPROP_LEN               },
  {1,MODIFIED_CE,         NODE_CW,        MODIFIED_SHIFT,         MODIFIED_LEN    },
  {1,NTYPE_CE,            NODE_CW,        NTYPE_SHIFT,            NTYPE_LEN               },

  {1,LINK_GEN,            LINK_CW,        GENERAL_SHIFT,          GENERAL_LEN             },
  {1,LOFFSET_CE,          LINK_CW,        LOFFSET_SHIFT,          LOFFSET_LEN             },

  {1,EDGE_GEN,            EDGE_CW,        GENERAL_SHIFT,          GENERAL_LEN             },
  {1,EOFFSET_CE,          EDGE_CW,        LOFFSET_SHIFT,          LOFFSET_LEN             },
  {0,0,0,0,0},
  {1,NOOFELEM_CE,         EDGE_CW,        NOOFELEM_SHIFT,         NOOFELEM_LEN    },
  {1,AUXEDGE_CE,          EDGE_CW,        AUXEDGE_SHIFT,          AUXEDGE_LEN             },
  {1,PATTERN_CE,          EDGE_CW,        PATTERN_SHIFT,          PATTERN_LEN             },
  {1,ADDPATTERN_CE,       EDGE_CW,        ADDPATTERN_SHIFT,       ADDPATTERN_LEN  },
  {1,EDGENEW_CE,          EDGE_CW,        EDGENEW_SHIFT,          EDGENEW_LEN             },

  {1,ELEMENT_GEN,         ELEMENT_CW,     GENERAL_SHIFT,          GENERAL_LEN             },
  {1,REFINE_CE,           ELEMENT_CW,     REFINE_SHIFT,           REFINE_LEN              },
  {1,ECLASS_CE,           ELEMENT_CW,     ECLASS_SHIFT,           ECLASS_LEN              },
  {1,NSONS_CE,            ELEMENT_CW,     NSONS_SHIFT,            NSONS_LEN               },
  {1,NEWEL_CE,            ELEMENT_CW,     NEWEL_SHIFT,        NEWEL_LEN           },
  {1,REFINECLASS_CE,      ELEMENT_CW,     REFINECLASS_SHIFT,      REFINECLASS_LEN },

  {1,MARK_CE,             FLAG_CW,        MARK_SHIFT,             MARK_LEN                },
  {1,COARSEN_CE,          FLAG_CW,        COARSEN_SHIFT,          COARSEN_LEN             },
  {1,EBUILDCON_CE,        FLAG_CW,        EBUILDCON_SHIFT,        EBUILDCON_LEN   },
  {1,DECOUPLED_CE,        FLAG_CW,        DECOUPLED_SHIFT,        DECOUPLED_LEN   },
  /* TODO: delete next line */
  /*	{1,EDGEPATTERN_CE,	FLAG_CW,	EDGEPATTERN_SHIFT,	EDGEPATTERN_LEN	}, */
  {1,SIDEPATTERN_CE,      FLAG_CW,        SIDEPATTERN_SHIFT,      SIDEPATTERN_LEN },
  {1,MARKCLASS_CE,        FLAG_CW,        MARKCLASS_SHIFT,        MARKCLASS_LEN   },


        #ifdef ModelP
  {1,XFERLINK_CE,         LINK_CW,        XFERLINK_SHIFT,         XFERLINK_LEN    },
  {1,XFERVECTOR_CE,       VECTOR_CW,      XFERVECTOR_SHIFT,       XFERVECTOR_LEN  },
  {1,XFERNODE_CE,         NODE_CW,        XFERNODE_SHIFT,         XFERNODE_LEN    },
  {1,XFERMATX_CE,         MATRIX_CW,      XFERMATX_SHIFT,         XFERMATX_LEN    },
        #else /* ModelP */
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
        #endif /* ModelP */
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
  {0,0,0,0,0},
}; /* last entry used: 72 for XFERLINK_CE */

/* free entry: 36 (EXTRA for edge does not exist any more) */

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

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

INT InitPredefinedControlEntries (void)
{
  INT i;
  CONTROL_ENTRY *ce;
  CONTROL_WORD *cw;
  predefined_control_entry *pce;

  for (i=0; i<MAX_CONTROL_ENTRIES; i++)
  {
    ce = control_entries+i;
    ce->used = 0;
    ce->control_word = 0;
    ce->offset_in_word = 0;
    ce->length = 0;
    ce->offset_in_object = 0;
    ce->mask = 0;
    ce->xor_mask = 0;
  }
  for (i=0; i<MAX_CONTROL_ENTRIES; i++)
    if (predefines[i].used)
    {
      pce = predefines+i;
      ce = control_entries+pce->control_entry_id;
      if (ce->used)
      {
        PrintErrorMessage('E',"InitPredefinedControlEntries","redefinition of control entry");
        return(__LINE__);
      }
      cw = control_words+pce->control_word;
      ce->used = 1;
      ce->control_word = pce->control_word;
      ce->offset_in_word = pce->offset_in_word;
      ce->length = pce->length;
      ce->offset_in_object = cw->offset_in_object;
      ce->mask = (POW2(ce->length)-1)<<ce->offset_in_word;
      ce->xor_mask = ~ce->mask;
                        #ifdef Debug
      /* check for overlapping control entries */
      if (cw->used_mask & ce->mask)
      {
        int j;
        CONTROL_ENTRY *test_ce;
        predefined_control_entry *test_pce;
        PRINTDEBUG(gm,1,("control_entry[%d] has overlapping bits with previous control_entries:\n",i));
        for (j=0; j<i; j++)
        {
          test_pce = predefines+j;
          test_ce = control_entries+test_pce->control_entry_id;
          if (test_ce->mask & ce->mask)
            PRINTDEBUG(gm,1,(" %d",j));

        }
        PRINTDEBUG(gm,1,("\n"));
      }
                        #endif
      cw->used_mask |= ce->mask;
    }
  return (GM_OK);
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

  /* free used bits */
  cw->used_mask &= ce->xor_mask;

  /* free control entry */
  ce->used = 0;

  /* ok, exit */
  return(GM_OK);
}
