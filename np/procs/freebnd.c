// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  freebnd.c														*/
/*																			*/
/* Purpose:   moving free boundaries		                                                                */
/*																			*/
/* Author:	  Henrik Rentz-Reichert			                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   01.09.97 begin, ug version 3.7								*/
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

#include "general.h"
#include "compiler.h"

/* dev */
#include "devices.h"

/* gm */
#include "gm.h"
#include "evm.h"

/* np */
#include "np.h"
#include "udm.h"

/* own header */
#include "freebnd.h"

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


/****************************************************************************/
/*D
   MoveFreeBoundary - move free boundary according to positions given in vd

   SYNOPSIS:
   INT MoveFreeBoundary (MULTIGRID *mg, const VECDATA_DESC *vd);

   PARAMETERS:
   .  mg - multigrid
   .  vd - data descriptor describing nodes only and having DIM components describing
                the new global coordinates of the node

   DESCRIPTION:
   All vertices of nodes are moved to the new global position given in vd.
   Vertices with move < DIM will be skipped automatically.
   CAUTION: you must not change the topology of the boundary triangulation! This
   will not be checked by the function.

   RETURN VALUE:
   INT
   .n		0:		ok
   .n		n>0:	else

   SEE ALSO:
   BNDP_Move
   D*/
/****************************************************************************/

INT MoveFreeBoundary (MULTIGRID *mg, INT level, const VECDATA_DESC *vd)
{
  VECTOR *vec;
  VERTEX *vert;
  INT lev;

        #ifdef ModelP
  /* TODO: parallel version */
  PrintErrorMessage('E',"MoveFreeBoundary","parallel not implemented");
  return (1);
        #endif

  if (VD_ncmps_in_otype_mod(vd,NODEVEC,NON_STRICT)<DIM)
    REP_ERR_RETURN(1);

  if (!VD_SUCC_COMP(vd))
    REP_ERR_RETURN(1);

  for (lev=0; lev<=level; lev++)
    for (vec=FIRSTVECTOR(GRID_ON_LEVEL(mg,lev)); vec!=NULL; vec=SUCCVC(vec))
      if ((lev==level) || FINE_GRID_DOF(vec))
      {
        vert = MYVERTEX((NODE*)VOBJECT(vec));
        if (!VD_ISDEF_IN_TYPE(vd,VTYPE(vec)))
          continue;
        if (OBJT(vert)!=BVOBJ)
          continue;
        if (MOVE(vert)!=DIM)
          continue;
        if (MoveFreeBoundaryVertex(mg,vert,VVALUEPTR(vec,VD_CMP_OF_TYPE(vd,VTYPE(vec),0))))
          REP_ERR_RETURN (1);
      }

  if (FinishMovingFreeBoundaryVertices(mg))
    REP_ERR_RETURN(1);

  return (0);
}

/****************************************************************************/
/*D
   StoreMGgeom - store all vertex local and global coordinates

   SYNOPSIS:
   INT StoreMGgeom (const MULTIGRID *mg, const VECDATA_DESC *vd)

   PARAMETERS:
   .  mg - multigrid
   .  vd - data descriptor describing ALL nodes having at least 2*DIM comps

   DESCRIPTION:
   This function stores all vertex local and global coordinates. For later
   restore by 'RestoreMGgeom'.

   RETURN VALUE:
   INT
   .n		0:		ok
   .n		n>0:	else

   SEE ALSO:
   RestoreMGgeom
   D*/
/****************************************************************************/

INT StoreMGgeom (const MULTIGRID *mg, const VECDATA_DESC *vd)
{
  NODE *nd;
  VECTOR *vec;
  DOUBLE *pos;
  INT l,vt;

  if (VD_ncmps_in_otype(vd,NODEVEC)<2*DIM)
    REP_ERR_RETURN (1);

  if (!VD_SUCC_COMP(vd))
    REP_ERR_RETURN(1);

  for (l=0; l<=TOPLEVEL(mg); l++)
    for (nd=FIRSTNODE(GRID_ON_LEVEL(mg,l)); nd!=NULL; nd=SUCCN(nd))
    {
      vec = NVECTOR(nd);
      ASSERT(vec!=NULL);
      vt = VTYPE(vec);
      ASSERT(VD_ISDEF_IN_TYPE(vd,vt));

      pos = CVECT(MYVERTEX(nd));
      V_DIM_COPY(pos,VVALUEPTR(vec,VD_CMP_OF_TYPE(vd,vt,0)));
      pos = LCVECT(MYVERTEX(nd));
      V_DIM_COPY(pos,VVALUEPTR(vec,VD_CMP_OF_TYPE(vd,vt,DIM)));
    }
  return (0);
}

/****************************************************************************/
/*D
   RestoreMGgeom - restore all vertex local and global coordinates

   SYNOPSIS:
   INT RestoreMGgeom (MULTIGRID *mg, const VECDATA_DESC *vd)

   PARAMETERS:
   .  mg - multigrid
   .  vd - data descriptor describing ALL nodes having at least 2*DIM comps

   DESCRIPTION:
   This function restores all vertex local and global coordinates previously
   stored by 'StoreMGgeom'.

   RETURN VALUE:
   INT
   .n		0:		ok
   .n		n>0:	else

   SEE ALSO:
   StoreMGgeom
   D*/
/****************************************************************************/

INT RestoreMGgeom (MULTIGRID *mg, const VECDATA_DESC *vd)
{
  VERTEX *vtx;
  NODE *nd;
  VECTOR *vec;
  INT l,vt;

  if (VD_ncmps_in_otype(vd,NODEVEC)<2*DIM)
    REP_ERR_RETURN (1);

  if (!VD_SUCC_COMP(vd))
    REP_ERR_RETURN(1);

  for (l=0; l<=TOPLEVEL(mg); l++)
    for (nd=FIRSTNODE(GRID_ON_LEVEL(mg,l)); nd!=NULL; nd=SUCCN(nd))
    {
      vec = NVECTOR(nd);
      ASSERT(vec!=NULL);
      vt = VTYPE(vec);
      ASSERT(VD_ISDEF_IN_TYPE(vd,vt));

      vtx = MYVERTEX(nd);

      if (MOVE(vtx)==DIM)
        if (SetVertexGlobalAndLocal(vtx,
                                    VVALUEPTR(vec,VD_CMP_OF_TYPE(vd,vt,0)),
                                    VVALUEPTR(vec,VD_CMP_OF_TYPE(vd,vt,DIM))))
          REP_ERR_RETURN (1);
    }
  return (0);
}
