// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef ModelP


#include <stdio.h>


#include "gm.h"
#include "parallel.h"


void ddd_test (int mode, MULTIGRID *theMG)
{
  ELEMENT   *elem;
  VECTOR    *vec;
  VSEGMENT  *vseg;
  MATRIX *mat;
  int level = 0;
  int nmat;

  InitMemMgr(theMG->theHeap);
  InitCurrMG(theMG);

  DDD_XferBegin();

  if (me==0)
  {
    for(elem=FIRSTELEMENT(GRID_ON_LEVEL(theMG,0)); elem!=NULL; elem=SUCCE(elem))
    {
      DDD_XferCopyObjX(DDD_OBJ(elem), 1, 7,(OBJT(elem)==BEOBJ) ? BND_SIZE(TAG(elem)) : INNER_SIZE(TAG(elem)));
    }
  }

  DDD_XferEnd();
}

/*
   void ddd_test (int mode, MULTIGRID *theMG)
   {
        ELEMENT   *elem;
        VECTOR    *vec;
        VSEGMENT  *vseg;
        MATRIX *mat;
        int       level = 0;
        int       nmat;

        InitMemMgr(theMG->theHeap);
        InitCurrMG(theMG);
        DDD_XferBegin();
                                                                                                                                                                        if (me==0)
        {
                for(vec=FIRSTVECTOR(GRID_ON_LEVEL(theMG,0)); vec!=NULL; vec=SUCCVC(vec))
                {
   #ifdef Debug
                        printf("%2d:  test(): v=%x Xfer size=%d\n",me,vec,theMG->theFormat->VectorSizes[VTYPE(vec)]);
   #endif
                        DDD_XferCopyObjX(DDD_OBJ(vec), 1, 7,theMG->theFormat->VectorSizes[VTYPE(vec)]);
                }
        }
        DDD_XferEnd();
   }
 */
/*
   void ddd_test (int mode, MULTIGRID *theMG)
   {
        ELEMENT   *elem;
        VECTOR    *vec;
        NODE      *node;
        LINK      *link;
        VSEGMENT  *vseg;
        MATRIX *mat;
        int       level = 0;
        int       nmat,nlink;

        InitMemMgr(theMG->theHeap);
        InitCurrMG(theMG);
        DDD_XferBegin();
                                                                                                                                                                        if (me==0)
        {
                for(node=FIRSTNODE(GRID_ON_LEVEL(theMG,0)); node!=NULL; node=SUCCN(node))
                {
                        DDD_XferCopyObj(DDD_OBJ(node), 1, 7);
                }
        }
        DDD_XferEnd();
   }

 */

#endif /* ModelP */
