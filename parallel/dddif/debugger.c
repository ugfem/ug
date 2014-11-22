// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/****************************************************************************/

#ifdef ModelP

#include <config.h>
#include <stdio.h>

#include "parallel.h"
#include "memmgr.h"
#include "general.h"
#include "ugm.h"
#include "ugdevices.h"
#include "namespace.h"

#include "debugger.h"

USING_UG_NAMESPACES
USING_PPIF_NAMESPACE

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/


/****************************************************************************/

void NS_DIM_PREFIX ddd_DisplayContext (void)
{
  int i, last=-1;
  char sep[2] = "";
  char buf[20];

  /* only master should display context */
  if (me!=master)
    return;

  UserWrite("current context: (");
  for(i=0; i<procs+1; i++)
  {
    if (i>=procs || !CONTEXT(i))
    {
      if (last+1==i-1)
      {
        sprintf(buf, "%s%d", sep, last+1); UserWrite(buf);
        sep[0] = ',';
      }
      if (last+1<i-1)
      {
        sprintf(buf, "%s%d-%d", sep, last+1, i-1); UserWrite(buf);
        sep[0] = ',';
      }
      last = i;
    }
  }
  UserWrite(")\n");
}


/****************************************************************************/
/*
   dddif_DisplayMemoryUsage -

   SYNOPSIS:
   static void dddif_DisplayMemoryUsage (void);

   PARAMETERS:
   .  void

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void dddif_DisplayMemoryUsage (void)
{
  memmgr_Report();

  UserWriteF("mem for interfaces:  %8ld bytes\n",
             (unsigned long) DDD_IFInfoMemoryAll()
             );

  UserWriteF("mem for couplings:   %8ld bytes\n",
             (unsigned long) DDD_InfoCplMemory()
             );
}


/****************************************************************************/
/*
   ddd_pstat -

   SYNOPSIS:
   void ddd_pstat (char *arg);

   PARAMETERS:
   .  arg

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/


void NS_DIM_PREFIX ddd_pstat (char *arg)
{
  int cmd;

  if (arg==NULL)
    return;

  cmd = arg[0];

  switch (cmd)
  {
  case 'X' :
    dddif_PrintGridRelations(dddctrl.currMG);
    break;

  case 'm' :
    SYNC_CONTEXT;
    dddif_DisplayMemoryUsage();
    SYNC_END;
    break;

  case 'c' :
    DDD_ConsCheck();
    UserWrite("\n");
    break;

  case 's' :
    SYNC_CONTEXT;
    DDD_Status();
    UserWrite("\n");
    SYNC_END;
    break;

  case 't' :
    if (me==master)
    {
      /* display ddd types */
      DDD_TypeDisplay(TypeVector);
      DDD_TypeDisplay(TypeIVertex);
      DDD_TypeDisplay(TypeBVertex);
      DDD_TypeDisplay(TypeNode);
                                #ifdef __THREEDIM__
      DDD_TypeDisplay(TypeEdge);
                                #endif

                                #ifdef __TWODIM__
      DDD_TypeDisplay(TypeTrElem);
      DDD_TypeDisplay(TypeTrBElem);
      DDD_TypeDisplay(TypeQuElem);
      DDD_TypeDisplay(TypeQuBElem);
                                #endif

                                #ifdef __THREEDIM__
      DDD_TypeDisplay(TypeTeElem);
      DDD_TypeDisplay(TypeTeBElem);
      DDD_TypeDisplay(TypePyElem);
      DDD_TypeDisplay(TypePyBElem);
      DDD_TypeDisplay(TypePrElem);
      DDD_TypeDisplay(TypePrBElem);
      DDD_TypeDisplay(TypeHeElem);
      DDD_TypeDisplay(TypeHeBElem);
                                #endif

      /* display dependent types */
      DDD_TypeDisplay(TypeMatrix);
                                #ifdef __TWODIM__
      DDD_TypeDisplay(TypeEdge);
                                #endif
    }
    break;


  case 'i' :
  {
    DDD_IF ifId = atoi(arg+1);

    SYNC_CONTEXT;
    if (ifId<=0)
      DDD_IFDisplayAll();
    else
      DDD_IFDisplay(ifId);

    UserWrite("\n");
    SYNC_END;
  }
  break;

  case 'l' :
    SYNC_CONTEXT;
    DDD_ListLocalObjects();
    UserWrite("\n");
    SYNC_END;
    break;

  case 'b' :
    buggy(dddctrl.currMG);
    UserWrite("BUGGY: returning control to caller\n");
    break;
  }
}


/****************************************************************************/

/*
        buggy - interactive debugging tool for distributed grids ug/ddd

        960603 kb  copied from pceq, modified, integrated
 */


/****************************************************************************/


/****************************************************************************/
/*
   buggy_ShowCopies -

   SYNOPSIS:
   static void buggy_ShowCopies (DDD_HDR hdr);

   PARAMETERS:
   .  hdr -

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void buggy_ShowCopies (DDD_HDR hdr)
{
  int   *p, i;

  p = DDD_InfoProcList(hdr);
  for(i=0; p[i]!=-1; i+=2)
  {
    printf("%4d:    copy on %3d with prio %d\n",
           me, p[i], p[i+1]);
  }
}



/****************************************************************************/


/****************************************************************************/
/*
   buggy_ElemShow -

   SYNOPSIS:
   static void buggy_ElemShow (ELEMENT *e);

   PARAMETERS:
   .  e -

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void buggy_ElemShow (ELEMENT *e)
{
  ELEMENT *SonList[MAX_SONS];
  int i;

  printf("%4d:    ID=%06d LEVEL=%02d corners=%03d\n", me,
         ID(e), LEVEL(e), CORNERS_OF_ELEM(e));

  if (EFATHER(e))
    printf("%4d:    father=%08x\n", me,
           DDD_InfoGlobalId(PARHDRE(EFATHER(e))));

  if (PREDE(e))
    printf("%4d:    pred=%08x\n", me,
           DDD_InfoGlobalId(PARHDRE(PREDE(e))));

  if (SUCCE(e))
    printf("%4d:    succ=%08x\n", me,
           DDD_InfoGlobalId(PARHDRE(SUCCE(e))));

  for(i=0; i<SIDES_OF_ELEM(e); i++)
  {
    if (NBELEM(e,i)!=NULL)
    {
      printf("%4d:    nb[%d]=%08x\n", me,
             i, DDD_InfoGlobalId(PARHDRE(NBELEM(e,i))));
    }
  }

  if (GetAllSons(e,SonList)==0)
  {
    for(i=0; SonList[i]!=NULL; i++)
    {
      printf("%4d:    son[%d]=%08x prio=%d\n", me,
             i,
             DDD_InfoGlobalId(PARHDRE(SonList[i])),
             DDD_InfoPriority(PARHDRE(SonList[i]))
             );
    }
  }
}


/****************************************************************************/
/*
   buggy_NodeShow -

   SYNOPSIS:
   static void buggy_NodeShow (NODE *n);

   PARAMETERS:
   .  n -

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void buggy_NodeShow (NODE *n)
{
  int i;

  printf("%4d:    ID=%06d LEVEL=%02d\n", me,
         ID(n), LEVEL(n));

  /* print coordinates of that node */
  printf("%4d:    VERTEXID=%06d LEVEL=%02d", me,
         ID(MYVERTEX(n)), LEVEL(MYVERTEX(n)));
  for(i=0; i<DIM; i++)
  {
    printf(" x%1d=%11.4E",i, (float)(CVECT(MYVERTEX(n))[i]) );
  }
  printf("\n");


  if (NFATHER(n))
    printf("%4d:    father=%08x\n", me,
           DDD_InfoGlobalId(PARHDR((NODE *)NFATHER(n))));

  if (PREDN(n))
    printf("%4d:    pred=%08x\n", me,
           DDD_InfoGlobalId(PARHDR(PREDN(n))));

  if (SUCCN(n))
    printf("%4d:    succ=%08x\n", me,
           DDD_InfoGlobalId(PARHDR(SUCCN(n))));
}


/****************************************************************************/
/*
   buggy_Search -

   SYNOPSIS:
   static void buggy_Search (MULTIGRID *theMG, DDD_GID gid);

   PARAMETERS:
   .  the MG -
   .  gid -

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void buggy_Search (MULTIGRID *theMG, DDD_GID gid)
{
  int level, found;

  found = FALSE;
  for(level=0; level<=TOPLEVEL(theMG); level++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,level);
    ELEMENT *e;
    NODE    *n;

    /* search for element */
    for(e=PFIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
    {
      if (DDD_InfoGlobalId(PARHDRE(e))==gid)
      {
        printf("%4d: ELEMENT gid=%08x, adr=%p, level=%d\n",
               me, gid, e, level);
        buggy_ShowCopies(PARHDRE(e));
        buggy_ElemShow(e);
        found = TRUE;
      }
    }


    /* search for node */
    for(n=PFIRSTNODE(theGrid); n!=NULL; n=SUCCN(n))
    {
      if (DDD_InfoGlobalId(PARHDR(n))==gid)
      {
        printf("%4d: NODE gid=%08x, adr=%p, level=%d\n",
               me, gid, n, level);
        buggy_ShowCopies(PARHDR(n));
        buggy_NodeShow(n);
        found = TRUE;
      }
    }
  }

  if (!found)
  {
    DDD_HDR hdr = DDD_SearchHdr(gid);

    if (hdr!=NULL)
    {
      printf("%4d: DDDOBJ gid=%08x, typ=%d, level=%d\n",
             me, gid, DDD_InfoType(hdr), DDD_InfoAttr(hdr));
      buggy_ShowCopies(hdr);
    }
    else
    {
      printf("%4d: unknown gid=%08x\n", me, gid);
    }
  }
}



/****************************************************************************/

/*
   static void bug_eg_item (Elem *item)
   {
        char  fname[30];
        FILE  *f;

        sprintf(fname,"buggy%08x.%03d.gnu",item->InfoGlobalId(),me);
        f = fopen(fname,"w");
        if (f==0)
        {
                Cout<<me<<": can't open "<<fname<<cr;
                return;
        }

        {
                ConnIterator(Side) siter(item->conn_dn());
                Side *side;
                for(; siter; ++siter)
                {
                        side = siter.rel()->obj();
                        Cout << me<<" side "<< hex << side->InfoGlobalId() <<dec<< cr;
                        {
                                ConnIterator(Node) niter(side->conn_dn());
                                Node *node;
                                for(; niter; ++niter)
                                {
                                        node = niter.rel()->obj();
                                        Cout << me <<"    node " << hex
                                                 << node->InfoGlobalId()
                                                 << dec << cr;
                                        fprintf(f,"%f %f\n",node->x[0],node->x[1]);
                                }
                                fprintf(f,"\n");
                        }
                }
        }

        fclose(f);
   }


   static void bug_elem_graphics (GridSegment &target, unsigned int gid)
   {
        Elem  *item=0;

        {
                Iterator(Elem) iter(target.elem_list());

                for (; iter && item==0; ++iter)
                {
                        if (iter.obj()->InfoGlobalId()==gid)
                        {
                                item = iter.obj();
                        }
                }
        }

        if (item==0)
        {
                Iterator(Elem) iter(target.ghostlist);

                for (; iter && item==0; ++iter)
                {
                        if (iter.obj()->InfoGlobalId()==gid)
                        {
                                item = iter.obj();
                        }
                }

        }

        if (item==0)
        {
                Cout<<me<<": no current elem/ghost"<<cr;
                return;
        }

        bug_eg_item(item);
   }
 */


/****************************************************************************/
/*
   buggy_help -

   SYNOPSIS:
   static void buggy_help (void);

   PARAMETERS:
   .  void

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

/****************************************************************************/

static void buggy_help (void)
{
  printf(
    " *\n"
    " * BUGGY ug debugger\n"
    " *\n"
    " *   x or q   quit\n"
    " *   p<no>    change current processor\n"
    " *   l        list DDD objects on current proc\n"
    " *   <gid>    change to object with gid\n"
    " *   ? or h   this help message\n"
    " *\n");
}



/****************************************************************************/


/****************************************************************************/
/*
   buggy -

   SYNOPSIS:
   void buggy (MULTIGRID *theMG);

   PARAMETERS:
   .  the MG

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void NS_DIM_PREFIX buggy (MULTIGRID *theMG)
{
  char buff[100];
  DDD_GID gid;
  int proc, cmd;

  Synchronize();

  if (me==0)
  {
    printf("%04d: started buggy.\n", me);
    fflush(stdout);
  }

  proc = 0;
  gid = 0x0;
  do
  {
    if (me==0)
    {
      do {
        printf("%04d: buggy> ", proc);
        fflush(stdout);
        scanf("%s", buff);
      } while (buff[0] == 0);

      switch (buff[0])
      {
      case 'x' :
      case 'q' :
        proc = -1;
        cmd = 0;
        break;

      case 'p' :
        proc = (int) strtol(buff+1, 0, 0);
        cmd = 1;
        break;

      case 'l' :
        cmd = 2;
        break;

      case '?' :
      case 'h' :
        cmd=99;
        break;

      default :
        cmd = 3;
        gid = /*(DDD_GID)*/ strtol(buff, 0, 0);
      }
    }

    Broadcast(&cmd, sizeof(int));
    Broadcast(&proc, sizeof(int));
    Broadcast(&gid, sizeof(unsigned int));

    if (me==proc)
    {
      switch (cmd)
      {
      case 99 :
        buggy_help();
        break;

      case 2 :
        DDD_ListLocalObjects();
        break;

      default :
        buggy_Search(theMG, gid);
        /*
                                                bug_ghost(target, gid);
                                                bug_side(target, gid);
         */
        break;
      }
    }

    fflush(stdout);
    Synchronize();
  }
  while (proc>=0);
}


/****************************************************************************/
/*                                                                          */
/* Function:  PrintGridRelations                                            */
/*                                                                          */
/* Purpose:   print certain information about a grid in order to test       */
/*            the formal approach to parallelisation.                       */
/*                                                                          */
/****************************************************************************/

#define PREFIX "__"

void NS_DIM_PREFIX dddif_PrintGridRelations (MULTIGRID *theMG)
{
  ELEMENT *e, *enb;
  GRID    *theGrid = GRID_ON_LEVEL(theMG,TOPLEVEL(theMG));
  INT j;

  SYNC_CONTEXT;
  for(e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
  {
    printf(PREFIX "master(e%07x, p%02d).\n", EGID(e), me);

    for(j=0; j<SIDES_OF_ELEM(e); j++)
    {
      enb = NBELEM(e,j);
      if (enb!=NULL)
      {
        printf(PREFIX "nb(e%07x, e%07x).\n", EGID(e), EGID(enb));
      }
    }
  }
  SYNC_END;

}

#undef PREFIX


/****************************************************************************/



#endif /* ModelP */
