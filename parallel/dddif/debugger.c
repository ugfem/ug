// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/****************************************************************************/

#ifdef ModelP

#include <stdio.h>

#include "parallel.h"
#include "general.h"

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)


/****************************************************************************/

void buggy (MULTIGRID *);


/****************************************************************************/

void ddd_DisplayContext (void)
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



void ddd_pstat (int cmd)
{
  switch (cmd)
  {
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

  case 'i' :
    SYNC_CONTEXT;
    DDD_IFDisplay();
    UserWrite("\n");
    SYNC_END;
    break;

  case 'l' :
    SYNC_CONTEXT;
    DDD_ListLocalObjects();
    UserWrite("\n");
    SYNC_END;
    break;

  case 'b' :
    buggy(dddctrl.currMG);
    UserWrite("BUGGY: aborted\n");
    break;
  }
}


/****************************************************************************/

/*
        buggy - interactive debugging tool for distributed grids ug/ddd

        960603 kb  copied from pceq, modified, integrated
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

static void buggy_ElemShow (ELEMENT *e)
{
  int i;

  if (EFATHER(e))
    printf("%4d:    father=%08x\n", me,
           DDD_InfoGlobalId(PARHDRE(EFATHER(e))));

  for(i=0; i<SIDES_OF_ELEM(e); i++)
  {
    if (NBELEM(e,i)!=NULL)
    {
      printf("%4d:    nb[%d]=%08x\n", me,
             i, DDD_InfoGlobalId(PARHDRE(NBELEM(e,i))));
    }
  }

  for(i=0; i<SONS_OF_ELEM(e); i++)
  {
    if (SON(e,i)!=NULL)
    {
      printf("%4d:    son[%d]=%08x\n", me,
             i, DDD_InfoGlobalId(PARHDRE(SON(e,i))));
    }
  }
}


static void buggy_ElemSearch (MULTIGRID *theMG, DDD_GID gid)
{
  int level;

  for(level=0; level<=TOPLEVEL(theMG); level++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,level);
    ELEMENT *e;
    for(e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
    {
      if (DDD_InfoGlobalId(PARHDRE(e))==gid)
      {
        printf("%4d: ELEMENT gid=%08x, adr=%08x, level=%d\n",
               me, gid, e, level);
        buggy_ShowCopies(PARHDRE(e));
        buggy_ElemShow(e);
      }
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
                                        fprintf(f,"%lf %lf\n",node->x[0],node->x[1]);
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

void buggy (MULTIGRID *theMG)
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
        gid = (DDD_GID) strtol(buff, 0, 0);
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
        buggy_ElemSearch(theMG, gid);
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



#endif /* ModelP */
