// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdio.h>
#include "../main/params.h"
#include "../main/defs.h"
#include "../main/structs.h"

void merge_goals(double *merged_goal, double *goal, int ndims, int part_type)
{
  char buffer[100];

  if (ndims==1)
  {
    if (part_type==BI_1 || part_type==BI_2)
    {
      merged_goal[0] = goal[0];
      merged_goal[1] = goal[1];
    }
  }
  else if (ndims==2)
  {
    if (part_type==QUAD_1)
    {
      merged_goal[0] = goal[0]+goal[2];
      merged_goal[1] = goal[1]+goal[3];
    }
    else if (part_type==QUAD_2 || part_type==QUAD_3 || part_type==QUAD_4)
    {
      merged_goal[0] = goal[0]+goal[1];
      merged_goal[1] = goal[2]+goal[3];
    }
  }
  else
  {
    sprintf(buffer,"FATAL: merge_goals() failed, ndims=%d\n", ndims);
  }
}

int partition_type(int xlen, int ylen, int ndims, int *array2proc,
                   int *xlens, int *ylens)
{
  int part_type;
  char buf[100];

  if (ndims == 1)
  {
    if (xlen>=ylen)
    {
      part_type     = BI_1;
      array2proc[1] = MED(xlen);
      xlens[0]      = MED(xlen);
      xlens[1]      = xlen-MED(xlen);
      ylens[0]      = ylen;
      ylens[1]      = ylen;
    }
    else if (xlen<ylen)
    {
      part_type = BI_2;
      array2proc[1] = MED(ylen)*xlen;
      xlens[0]      = xlen;
      xlens[1]      = xlen;
      ylens[0]      = MED(ylen);
      ylens[1]      = ylen-MED(ylen);
    }
    else
    {
      sprintf(buf,"FATAL: partition_type failed, ndims=%d\n", ndims);
      UserWrite(buf);
    }
  }
  else if (ndims == 2)
  {
    if ((xlen>=ylen)&&(xlen<=3*ylen))
    {
      part_type = QUAD_1;
      array2proc[1] = MED(xlen);
      array2proc[2] = MED(ylen)*xlen;
      array2proc[3] = MED(ylen)*xlen+MED(xlen);
      xlens[0]      = MED(xlen);
      xlens[1]      = xlen-MED(xlen);
      xlens[2]      = MED(xlen);
      xlens[3]      = xlen-MED(xlen);
      ylens[0]      = MED(ylen);
      ylens[1]      = MED(ylen);
      ylens[2]      = ylen-MED(ylen);
      ylens[3]      = ylen-MED(ylen);
    }
    else if ((xlen<ylen)&&(3*xlen>=ylen))
    {
      part_type = QUAD_2;
      array2proc[1] = MED(xlen);
      array2proc[2] = MED(ylen)*xlen;
      array2proc[3] = MED(ylen)*xlen+MED(xlen);
      xlens[0]      = MED(xlen);
      xlens[1]      = xlen-MED(xlen);
      xlens[2]      = MED(xlen);
      xlens[3]      = xlen-MED(xlen);
      ylens[0]      = MED(ylen);
      ylens[1]      = MED(ylen);
      ylens[2]      = ylen-MED(ylen);
      ylens[3]      = ylen-MED(ylen);
    }
    else if (xlen>3*ylen)
    {
      part_type = QUAD_3;
      array2proc[1] = MED(MED(xlen));
      array2proc[2] = MED(xlen);
      array2proc[3] = MED(xlen-MED(xlen));
      xlens[0]      = MED(MED(xlen));
      xlens[1]      = MED(xlen)-MED(MED(xlen));
      xlens[2]      = MED(xlen-MED(xlen));
      xlens[3]      = xlen-MED(xlen)-MED(xlen-MED(xlen));
      ylens[0]      = ylen;
      ylens[1]      = ylen;
      ylens[2]      = ylen;
      ylens[3]      = ylen;
    }
    else if (3*xlen<ylen)
    {
      part_type = QUAD_4;
      array2proc[1] = MED(MED(ylen));
      array2proc[2] = MED(ylen);
      array2proc[3] = MED(ylen-MED(ylen));
      xlens[0]      = xlen;
      xlens[1]      = xlen;
      xlens[2]      = xlen;
      xlens[3]      = xlen;
      ylens[0]      = MED(MED(ylen));
      ylens[1]      = MED(ylen)-MED(MED(ylen));
      ylens[2]      = MED(ylen-MED(ylen));
      ylens[3]      = ylen-MED(ylen)-MED(ylen-MED(ylen));
    }
    else
    {
      sprintf(buf,"FATAL: partition_type failed, ndims=%d\n", ndims);
      UserWrite(buf);
    }
  }
  else if (ndims == 3)
  {
    if (xlen>ylen && xlen<7*ylen)
    {
      part_type = 1;
      array2proc[1] = MED(MED(ylen));
      array2proc[2] = MED(ylen);
      array2proc[3] = MED(ylen-MED(ylen));
      array2proc[4] = MED(MED(ylen));
      array2proc[5] = MED(ylen);
      array2proc[6] = MED(ylen-MED(ylen));
      array2proc[7] = MED(ylen-MED(ylen));
      xlens[0]      = xlen;
      xlens[1]      = xlen;
      xlens[2]      = xlen;
      xlens[3]      = xlen;
      xlens[4]      = xlen;
      xlens[5]      = xlen;
      xlens[6]      = xlen;
      xlens[7]      = xlen;
      ylens[0]      = MED(MED(ylen));
      ylens[1]      = MED(ylen)-MED(MED(ylen));
      ylens[2]      = MED(ylen-MED(ylen));
      ylens[3]      = ylen-MED(ylen)-MED(ylen-MED(ylen));
      ylens[4]      = MED(MED(ylen));
      ylens[5]      = MED(ylen)-MED(MED(ylen));
      ylens[6]      = MED(ylen-MED(ylen));
      ylens[7]      = ylen-MED(ylen)-MED(ylen-MED(ylen));
    }
    else
    {
      sprintf(buf,"FATAL: partition_type failed, ndims=%d\n", ndims);
      UserWrite(buf);
    }
  }
  else return(-1);


  return(part_type);
}


/* map() maps array positions to single numbers */

int map (int x, int y, int ndims, int part_type,
         int *xlens, int *ylens)
{
  int i;
  int len;
  int xpos,ypos;

  if (ndims == 1)
  {
    if (part_type==BI_1)
      if (x<xlens[0]) return(0);
      else return(1);
    if (part_type==BI_2)
      if (y<ylens[0]) return(0);
      else return(1);
  }
  else if (ndims == 2)
  {
    if (part_type==QUAD_1 || part_type==QUAD_2)
    {
      if (x<xlens[0]) xpos=0;
      else xpos=1;
      if (y<ylens[0]) ypos=0;
      else ypos=1;

      return(2*ypos+xpos);
    }
    else if (part_type==QUAD_3)
    {
      len = xlens[0];
      for (i=0; i<3; i++)
      {
        if (x<len) break;
        len += xlens[i+1];
      }

      return(i);
    }
    else if (part_type==QUAD_4)
    {
      len = ylens[0];
      for (i=0; i<3; i++)
      {
        if (y<len) break;
        len += ylens[i+1];
      }

      return(i);
    }
  }
}

/* map the subarray position back to array position */

int suba2a(short *assignment, int *loc2glob, int i, short *subassign, int ndims,
           int xlen, int ylen, int set, int *xlens, int *ylens, int part_type)
{
  extern int DEBUG_ARRAY;
  short oldass,newass;
  int s;
  char buffer[100];

  if (loc2glob!=NULL)
  {
    oldass=assignment[loc2glob[i]];
  }
  else
    oldass=0;
  s=*subassign;

  if (ndims==1)
  {
    if (part_type==BI_1)
    {
      if (set==0)
        newass = s/xlens[set]*xlen+s%xlens[set];
      else if (set==1)
        newass = s/xlens[set]*xlen+xlens[0]+s%xlens[set];
    }
    else if (part_type==BI_2)
    {
      if (set==0)
        newass = s;
      else if (set==1)
        newass = (s/xlens[set]+ylens[0])*xlen+s%xlens[set];
    }
    else
    {
      sprintf(buffer,"FATAL: suba2a() failed ndims=%d part_type=%d\n",
              ndims, part_type);
      UserWrite(buffer);
    }
  }
  else if (ndims==2)
  {
    if (part_type==QUAD_1 || part_type==QUAD_2)
    {
      if (set==0)
        newass = s/xlens[set]*xlen+s%xlens[set];
      else if (set==1)
        newass = s/xlens[set]*xlen+xlens[0]+s%xlens[set];
      else if (set==2)
        newass = (s/xlens[set]+ylens[0])*xlen+s%xlens[set];
      else if (set==3)
        newass = (s/xlens[set]+ylens[0])*xlen+xlens[2]+s%xlens[set];
    }
    else if (part_type==QUAD_3)
    {
      if (set==0)
        newass = s/xlens[set]*xlen+s%xlens[set];
      else if (set==1)
        newass = s/xlens[set]*xlen+xlens[0]+s%xlens[set];
      else if (set==2)
        newass = s/xlens[set]*xlen+xlens[0]+xlens[1]+s%xlens[set];
      else if (set==3)
        newass = s/xlens[set]*xlen+xlens[0]+xlens[1]+xlens[2]+s%xlens[set];
    }
    else if (part_type==QUAD_4)
    {
      if (set==0)
        newass = s/xlens[set]*xlen+s%xlens[set];
      else if (set==1)
        newass= (s/xlens[set]+ylens[0])*xlen+s%xlens[set];
      else if (set==2)
        newass = (s/xlens[set]+ylens[0]+ylens[1])*xlen+s%xlens[set];
      else if (set==3)
        newass = (s/xlens[set]+ylens[0]+ylens[1]+ylens[2])*xlen
                 +s%xlens[set];
    }
    else
    {
      sprintf(buffer,"FATAL: suba2a() failed ndims=%d part_type=%d\n",
              ndims, part_type);
      UserWrite(buffer);
    }
  }
  else
  {
    sprintf(buffer,"FATAL: suba2a() failed ndims=%d part_type=%d\n",
            ndims, part_type);
    UserWrite(buffer);
  }

# ifdef __DEBUG__
  if (DEBUG_ARRAY>1)
  {
    sprintf(buffer,"suba2a: newass=%d oldass=%d subass=%d xlen=%d subylen=%d xlen=%d ylen=%d\n",
            newass,oldass,s,xlens[set],ylens[set],xlen,ylen);
    UserWrite(buffer);
  }
# endif
  return(newass);
}



/* compute partitioning dimension for this set in */
/* the next recursion step                        */
int compute_sub_ndims (int xlen, int ylen, int ndims)
{
  int procs;

  procs = xlen*ylen;

  if (ndims==1)
  {
    if (procs>=2) return(1);
    else return(0);
  }
  else if (ndims==2)
  {
    if (procs>=4) return(2);
    else if (procs>=2) return(1);
    else return(0);
  }
  else if (ndims==3)
  {
    if (procs>=8) return(3);
    else if (procs>=4) return(2);
    else if (procs>=2) return(1);
    else return(0);
  }
}
