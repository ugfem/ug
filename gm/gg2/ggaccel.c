// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ggaccel.c                                                     */
/*                                                                          */
/* Purpose:   accelerating grid generator					                        */
/*                                                                          */
/* Author:    Dirk Feuchter			                                        */
/*                                                                          */
/* History:   10.05.95					                                        */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/




/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include  <stdlib.h>
#include  "ggaccel.h"
#include  "misc.h"
#include  "gm.h"
#include  "ugm.h"
#include  "ggm.h"


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/
static QUADTREETYP *startpointer;

static MULTIGRID *MG;

static COORD startwidth;

static SOURCETYP *source;

static BALTREETYP **q, *btree_rootpointer;

static int del_edg_fnd;

static MG_GGDATA *myMGdata;

static GG_PARAM *myPars;

/* data for CVS */
static char rcsid[] = "$Header$";

static INT QuObj;
static INT ScObj;
static INT QfclObj;
static INT EttObj;


/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* Function:  GetGGObjIDs		                                                */
/*                                                                          */
/* Purpose:   necessary for dynamic object data in ug 3.0					*/
/*            the different data structures for the grid generator get their*/
/*			  id«s															*/
/*                                                                          */
/*                                                                          */
/****************************************************************************/

static INT InitAccelObjs (MULTIGRID *theMG)
{

  Mark(MGHEAP(theMG),FROM_TOP);

  QuObj    = GetFreeOBJT();
  ScObj    = GetFreeOBJT();
  QfclObj  = GetFreeOBJT();
  EttObj   = GetFreeOBJT();

  return (0);
}


/****************************************************************************/
/*                                                                          */
/* Function:  TerminateAccel		                                        */
/*                                                                          */
/* Purpose:                                                                                                                             */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

INT TerminateAccel (MULTIGRID *theMG, INT flag)
{

  theMG->freeObjects[QuObj] = NULL;
  theMG->freeObjects[ScObj]  = NULL;
  theMG->freeObjects[QfclObj]  = NULL;
  theMG->freeObjects[EttObj]  = NULL;

  ReleaseOBJT(QuObj);
  ReleaseOBJT(ScObj);
  ReleaseOBJT(QfclObj);
  ReleaseOBJT(EttObj);

  Release(MGHEAP(theMG),FROM_TOP);

  return (0);
}


/****************************************************************************/
/*                                                                          */
/* Function:  length_of_edge                                                    */
/*                                                                          */
/* Purpose:   evaluates length of an edge									*/
/*                                                                          */
/* Input:     FRONTCOMP* edgefc, edgefc_aim: frontcomponents of the edge    */
/*                                                                          */
/* Output:    returns length of the edge between edgefc and edgefc_aim      */
/*                                                                          */
/****************************************************************************/

static float length_of_edge(FRONTCOMP* edgefc, FRONTCOMP* edgefc_aim)
{
  VERTEX *theVertex;
  float xc,yc,sxc,syc,dist2;


  theVertex = MYVERTEX(FRONTN(edgefc));
  sxc = XC(theVertex);
  syc = YC(theVertex);

  theVertex = MYVERTEX(FRONTN(edgefc_aim));
  xc = XC(theVertex);
  yc = YC(theVertex);

  dist2 = (xc-sxc)*(xc-sxc)+(yc-syc)*(yc-syc);
  return(dist2);
}

/****************** end of function length_of_edge **************************/





/****************************************************************************/
/*                                                                          */
/* Function:  length_of_angle                                                   */
/*                                                                          */
/* Purpose:   evaluates size of an Advancing Front angle					*/
/*                                                                          */
/* Input:     FRONTCOMP* anglefc_ori, anglefc, anglefc_aim: frontcomponents */
/*            which belong to an angle of the Advancing Front               */
/*                                                                          */
/* Output:    returns size of the angle between anglefc_ori, anglefc and    */
/*            anglefc_aim. The angles {0¡,...,90¡,...180¡,...,270¡,...360¡} */
/*			  are given by floats     {-1.0,..,0.0,...,1.0,...,2.0,...,3.0} */
/*                                                                          */
/****************************************************************************/

float length_of_angle
  (FRONTCOMP* anglefc_ori, FRONTCOMP* anglefc, FRONTCOMP* anglefc_aim)
{
  VERTEX *theVertex;
  COORD xc,yc,px,py,sx,sy,angle;

  theVertex = MYVERTEX(FRONTN(anglefc));
  xc = XC(theVertex);
  yc = YC(theVertex);
  theVertex = MYVERTEX(FRONTN(anglefc_ori));
  px = xc - XC(theVertex);
  py = yc - YC(theVertex);
  theVertex = MYVERTEX(FRONTN(anglefc_aim));
  sx = XC(theVertex) - xc;
  sy = YC(theVertex) - yc;

  /* angle > 180 degrees? */
  if ((py*sx-px*sy) > SMALLCOORD)                       /* check (py*sx-px*sy)>0 */
  {
    angle = (px*sx+py*sy)/sqrt((sx*sx+sy*sy)*(px*px+py*py));
    angle = 2.0 - angle;
    return(angle);
  }

  angle = (px*sx+py*sy)/sqrt((sx*sx+sy*sy)*(px*px+py*py));
  return(angle);

}

/********************* end of function length_of_edge ***********************/








/****************************************************************************/
/*                                                                          */
/* Function:  PointInTriangle                                                   */
/*                                                                          */
/* Purpose:   determines whether a node lies within in a triangle plus an   */
/*            epsilon environment											*/
/*                                                                          */
/* Input:     COORD x[3], COORD y[3]: coordinates of the triangle                       */
/*			  COORD pt[DIM] : coordinates of the node                                               */
/*                                                                          */
/* Output:    YES or NO, it depends whether the node is within or out of    */
/*			  the triangle													*/
/*                                                                          */
/****************************************************************************/

static INT PointInTriangle (COORD pt[DIM], COORD x[3], COORD y[3])
{
  int i,index;
  COORD sx,sy,hlp;
  COORD eps;
  COORD x_eps[3], y_eps[3];


  /* Increasing of the search-triangle with eps to find nodes, */
  /* which lie on its edges */

  eps = myPars->epsi;


  for (i=0; i<3; i++)
  {
    x_eps[i] = x[i];
    y_eps[i] = y[i];
  }


  /* find x_min */
  if (x_eps[0] <  x_eps[1])
    index = 0;
  else
    index = 1;

  if  (x_eps[index] < x_eps[2] )
    x_eps[index] -= eps;
  else
    x_eps[2] -= eps;

  /* find y_min */
  if (y_eps[0] <  y_eps[1])
    index = 0;
  else
    index = 1;

  if  (y_eps[index] < y_eps[2] )
    y_eps[index] -= eps;
  else
    y_eps[2] -= eps;

  /* find x_max */
  if (x_eps[0] >  x_eps[1])
    index = 0;
  else
    index = 1;

  if  (x_eps[index] > x_eps[2] )
    x_eps[index] += eps;
  else
    x_eps[2] += eps;

  /* find y_max */
  if (y_eps[0] >  y_eps[1])
    index = 0;
  else
    index = 1;

  if  (y_eps[index] > y_eps[2] )
    y_eps[index] += eps;
  else
    y_eps[2] += eps;

  for (i=0; i<3; i++)
  {
    sx = x_eps[(i+1)%3] - x_eps[i];
    sy = y_eps[(i+1)%3] - y_eps[i];
    hlp = ((sy*(pt[0]-x_eps[i])-sx*(pt[1]-y_eps[i]))/(sx*sx+sy*sy));
    if (hlp > SMALLCOORD)
      return (NO);
  }
  return (YES);

}

/****************** end of function PointInTriangle *************************/







/****************************************************************************/
/*                                                                          */
/* Function:  PointInCircle                                                             */
/*                                                                          */
/* Purpose:   determines whether a node lies within in a circle                         */
/*                                                                          */
/* Input:     COORD x, COORD y: coordinates of the middle of the circle         */
/*			  COORD pt[DIM] : coordinates of the node                                               */
/*                                                                          */
/* Output:    YES or NO, it depends whether the node is within or out of    */
/*			  the triangle													*/
/*                                                                          */
/****************************************************************************/

static INT PointInCircle (COORD pt[DIM], COORD x, COORD y, COORD searchrad2)
{
  /* searchrad2 is the squared radius */

  COORD dx,dy;

  dx = pt[0] - x;
  dy = pt[1] - y;

  if (searchrad2-(dx*dx+dy*dy)>SMALLCOORD )
    /* check of searchrad2>(dx*dx+dy*dy) */
    return (YES);
  else
    return (NO);
}

/******************* end of function PointInCircle **************************/







/****************************************************************************/
/*                                                                          */
/* Function:  distance                                                                  */
/*                                                                          */
/* Purpose:   evaluates the distance between two nodes                                          */
/*                                                                          */
/* Input:     FRONTCOMP *p1, FRONTCOMP *p2: the two nodes                                       */
/*                                                                          */
/* Output:    the squared distance in floating point                            */
/*                        Using the square root is not necessary						*/
/*                                                                          */
/****************************************************************************/

static COORD distance(FRONTCOMP *p1, FRONTCOMP *p2)
{
  COORD dist;

  dist = pow( (XC(MYVERTEX(FRONTN(p2))) - XC(MYVERTEX(FRONTN(p1)))), 2 ) +
         pow( (YC(MYVERTEX(FRONTN(p2))) - YC(MYVERTEX(FRONTN(p1)))), 2 ) ;
  return(dist);
}

/******************* end of function PointInCircle **************************/




/****************************************************************************/
/*                                                                          */
/* Function:  search                                                                    */
/*                                                                          */
/* Purpose:   runs through the quadtree and provides the essential pointer  */
/*			  concerning the frontcomponent new_p							*/
/*			  the function calls itself (recursively used)					*/
/*                                                                          */
/* Input:     QUADTREETYP *q_pointer: pointer at a quadtree structure       */
/*            SOURCETYP *so: pointer at a source respectivley a left        */
/*                           inferior corner								*/
/*			  COORD *wi: the width of one of the four subquadrangles                */
/*            FRONTCOMP *new_p: pointer at the new frontcomponent			*/
/*                                                                          */
/* Output:    pointer on a quadtree structure                                           */
/*                                                                          */
/****************************************************************************/

static QUADTREETYP* search(QUADTREETYP *q_pointer, SOURCETYP *so, COORD *wi,
                           FRONTCOMP *new_p)
{
  unsigned char helpc;

  if ( YC(MYVERTEX(FRONTN(new_p))) <  ( so->y + *wi) )
  {
    if ( XC(MYVERTEX(FRONTN(new_p))) < ( so->x + *wi) )
    {
      /* third quadrant  */
      helpc = q_pointer->q_flag;
      helpc <<= 7;                   /* Bit 0 */
      if ( helpc == 0 )                   /* Quadreepointer  !!! */
      {
        *wi = (*wi)/2;
        q_pointer = search( q_pointer->q_array[0], so, wi, new_p);
        return(q_pointer);
      }
      else                   /* Nodepointer !!! */
      {
        q_pointer->q_flag <<= 4 ;
        q_pointer->q_flag >>= 4 ;                         /* Highpart := 0 !!! */
        helpc = 0;                         /* 0 for third quadrant  */
        helpc <<= 4;
        q_pointer->q_flag |= helpc;
        return(q_pointer);
      }
    }             /* if ( third quadrant) */

    else
    {
      /* second quadrant  */
      helpc = q_pointer->q_flag;
      helpc >>= 1;  helpc <<= 7;                   /* Bit 1 */
      if ( helpc == 0 )                   /* Quadtreepointer !!! */
      {
        so->x = so->x + *wi;
        *wi = (*wi)/2;
        q_pointer = search( q_pointer->q_array[1], so, wi, new_p );
        return(q_pointer);
      }
      else
      {
        q_pointer->q_flag <<= 4 ; q_pointer->q_flag >>= 4;
        helpc = 1;                         /* 1 for second quadrant  */
        helpc <<=4 ;
        q_pointer->q_flag |= helpc;
        return(q_pointer);
      }
    }              /* else (second quadrant) */
  }        /* if (third or second quadrant) */
  else
  {
    if ( XC(MYVERTEX(FRONTN(new_p))) >= ( so->x + *wi) )
    {

      /* first quadrant  */

      helpc = q_pointer->q_flag;
      helpc >>= 2;  helpc <<= 7;                   /* Bit 2 */
      if ( helpc == 0 )                   /* Quadtreepointer exists !!! */
      {
        so->x = so->x + *wi;
        so->y = so->y + *wi;
        *wi = (*wi)/2;
        q_pointer = search( q_pointer->q_array[2], so, wi, new_p );
        return(q_pointer);
      }
      else
      {
        q_pointer->q_flag <<= 4 ; q_pointer->q_flag >>= 4;
        helpc = 2;                         /* 2 for first quadrant  */
        helpc <<=4 ;
        q_pointer->q_flag |= helpc;
        return(q_pointer);
      }

    }              /* if first quadrant */

    else
    {
      /* fourth quadrant  */

      helpc = q_pointer->q_flag;
      helpc >>= 3;  helpc <<= 7;                   /* Bit 3 */
      if ( helpc == 0 )                   /* Quadtreepointer !!! */
      {
        so->y = so->y + *wi;
        *wi = (*wi)/2;
        q_pointer = search( q_pointer->q_array[3], so, wi, new_p );
        return(q_pointer);
      }
      else
      {
        q_pointer->q_flag <<= 4 ; q_pointer->q_flag >>= 4;
        helpc = 3;                         /* 3 for fourth quadrant  */
        helpc <<=4 ;
        q_pointer->q_flag |= helpc;
        return(q_pointer);
      }
    }             /* else (fourth quadrant ) */

  }        /* else (first or fourth quadrant */
}

/************************ end of function search ****************************/




/****************************************************************************/
/*                                                                          */
/* Function:  environment_search                                                */
/*                                                                          */
/*                                                                          */
/* Purpose:   searchs for frontcomponents wich are close or                             */
/*            within the suggested triangle									*/
/*                                                                          */
/* Input:     INDEPFRONTLIST *theIFL: the Independent Front List, the new   */
/*									  suggested triangle belongs to         */
/*            QUADTREETYP *q_pointer: pointer at a quadtree structure       */
/*            SOURCETYP *so: pointer at a source respectivley a left        */
/*                                    inferior corner					    */
/*			  FRONTCOMP* thefoundPoints[MAXNPOINTS] : Array for the FC«s,   */
/*			                          found within the small searching      */
/*									  rectangle								*/
/*			  FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS] :  Array for   */
/*									  the FC«s, found within the big        */
/*									  searching rectangle - necessary to    */
/*									  detect long cutting edges             */
/*            SOURCETYP *search_sq_ld,*search_sq_ru,*big_search_sq_ld,      */
/*                      *big_search_sq_ru : coordinates of the different    */
/*									  sources								*/
/*            COORD xt[3], COORD yt[3] : coordinates of the sugg. triangle	*/
/*			  COORD searchradis : radius of the circle around the new       */
/*									  suggested node to detect also close   */
/*									  cutting edges							*/
/*            int *foundpoints, int *ii: reference parameters for number of */
/*									  found points                                                  */
/*																			*/
/* Output:	  FRONTCOMP* thefoundPoints[MAXNPOINTS] : Array with the FC«s,  */
/*			                          found within the small searching      */
/*									  rectangle								*/
/*			  FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS] :  Array with  */
/*									  the FC«s, found within the big        */
/*									  searching rectangle - necessary to    */
/*									  detect long cutting edges             */
/*			                          found within the small searching      */
/*									  rectangle								*/
/*            int *foundpoints : Number of  FC«s, found within the small    */
/*			                          searching rectangle					*/
/*            int *ii          : Number of  FC«s, found within the big      */
/*			                          searching rectangle					*/
/*																			*/
/****************************************************************************/

static void environment_search(INDEPFRONTLIST *theIFL,
                               QUADTREETYP *q_pointer, SOURCETYP* so,
                               FRONTCOMP *thefoundPoints[MAXNPOINTS],
                               FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS],
                               COORD wi,
                               SOURCETYP *search_sq_ld, SOURCETYP *search_sq_ru,
                               SOURCETYP *big_search_sq_ld,
                               SOURCETYP *big_search_sq_ru,
                               COORD xt[3], COORD yt[3],
                               COORD searchradis, int *foundpoints, int *ii)
{
  QFCLISTTYP *hn_pointer;
  int i,touch_possible;
  unsigned char helpc;
  COORD pt[DIM];


  for(i=0; i<=3; i++)
  {
    switch (i)
    {
    case 0 : break;
    case 1 : so->x += wi;break;
    case 2 : so->y += wi;break;
    case 3 : so->x -= wi;break;
    }

    touch_possible = 1;
    if (  ( so->x + wi < big_search_sq_ld->x ) ||
          ( so->x > big_search_sq_ru->x )      ||
          ( so->y + wi < big_search_sq_ld->y ) ||
          ( so->y > big_search_sq_ru->y )
          )
      touch_possible = 0;


    if (touch_possible)
    {
      /*  touch, intersection or covering                */
      helpc = q_pointer->q_flag;
      /* Is there any node or do we have to go deaper ?   */
      helpc >>= i; helpc <<= 7;
      if (helpc == 0)                    /* if "quadtreepointer" */
      {
        wi = wi/2;
        environment_search( theIFL, q_pointer->q_array[i], so,
                            thefoundPoints, theIntersectfoundPoints,
                            wi, search_sq_ld, search_sq_ru,
                            big_search_sq_ld, big_search_sq_ru,
                            xt, yt, searchradis, foundpoints, ii );
        so->y -= wi;
        wi = wi*2;
      }

      else
      {
        hn_pointer = q_pointer->q_array[i];
        /* this is the way, how the nodepointer of the quadtree becomes */
        /*  a void-pointer !!!                                    */

        /* Does this square have a node  ?     */
        /* Is it necessary to take it in consideration ???  */
        if (hn_pointer != NULL)
        {

          /* Is the node within the  big searching-square    */
          /* maximal 4 comparisons !!!                     */

          if (   ( XC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                   >= big_search_sq_ld->x  )                     &&

                 ( XC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                   <= big_search_sq_ru->x  )                     &&

                 ( YC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                   >= big_search_sq_ld->y  )                     &&

                 ( YC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                   <= big_search_sq_ru->y  )
                 )
          {
            if (   ( XC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                     < search_sq_ld->x  )                      ||

                   ( XC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                     > search_sq_ru->x  )                      ||

                   ( YC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                     < search_sq_ld->y  )                      ||

                   ( YC(MYVERTEX(FRONTN(FROC(hn_pointer))))
                     > search_sq_ru->y  )
                   )
            /* == "but the node is not within the small */
            /* searching-square" !!!   */
            {
              /* what«s about its successor ? */
              /* so we«ve found a candidate for possible */
              /* "FrontLineIntersction" */
              while ( hn_pointer != NULL)
              {
                if ( MYIFL(MYFL(FROC(hn_pointer))) == theIFL )
                {
                  theIntersectfoundPoints[*ii] = FROC(hn_pointer);
                  (*ii)++;
                  if
                  (

                    (XC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                     <  big_search_sq_ld->x  )  ||

                    (XC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                     >  big_search_sq_ru->x  )  ||

                    (YC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                     <  big_search_sq_ld->y  )  ||

                    (YC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                     >  big_search_sq_ru->y  )

                  )

                  {
                    theIntersectfoundPoints[*ii] =
                      PREDFC(FROC(hn_pointer));

                    (*ii)++;
                  }
                }
                hn_pointer = NXT(hn_pointer);
              }

            }                                     /* of "if the node is not within the samll s..." */

            else                                     /* the node is within the small searching-square */
            {
              pt[0] = XC(MYVERTEX(FRONTN(FROC(hn_pointer))));
              pt[1] = YC(MYVERTEX(FRONTN(FROC(hn_pointer))));
              if ( ( (pt[0] == xt[0] ) &&
                     (pt[1] == yt[0] ) )  ||

                   ( (pt[0] == xt[1] ) &&
                     (pt[1] == yt[1] ) )  ||

                   ( (pt[0] == xt[2] ) &&
                     (pt[1] == yt[2] ) )
                   )

              {
                continue;                                                 /*it is a node of the triangle itself*/
              }




              else
              {
                /* this node must be considered !!!*/
                if ( PointInTriangle    (pt,xt,yt) ||
                     PointInCircle (pt,xt[2],yt[2],
                                    searchradis*searchradis )
                     )
                {
                  while (hn_pointer != NULL)
                  {
                    if ( MYIFL(MYFL(FROC(hn_pointer)))
                         == theIFL )
                    {
                      thefoundPoints[*foundpoints] =
                        FROC(hn_pointer);

                      (*foundpoints)++;
                    }
                    hn_pointer = NXT(hn_pointer);
                  }
                }                                                 /* of in_triangle */

                else
                {
                  /* what«s about its successor ?   */
                  /* Well It is in any case within the big    */
                  /* searching-square therefore no comparisons */
                  /* necessary!!!                              */
                  /* so we«ve found a candidate for possible  */
                  /* "FrontLineIntersection" */

                  while ( hn_pointer != NULL)
                  {
                    if ( MYIFL(MYFL(FROC(hn_pointer)))
                         == theIFL )
                    {
                      theIntersectfoundPoints[*ii] =
                        FROC(hn_pointer);

                      (*ii)++;

                      if
                      (

                        (XC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                         <  big_search_sq_ld->x  )  ||

                        (XC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                         >  big_search_sq_ru->x  )  ||

                        (YC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                         <  big_search_sq_ld->y  )  ||

                        (YC(MYVERTEX(FRONTN(PREDFC(FROC(hn_pointer)))))
                         >  big_search_sq_ru->y  )

                      )

                      {

                        theIntersectfoundPoints[*ii] =
                          PREDFC(FROC(hn_pointer));

                        (*ii)++;
                      }
                    }
                    hn_pointer = NXT(hn_pointer);
                  }

                }
              }                                           /* else it is not a node of the triangle itself */
            }                                     /* von else: the node is within the small s... !!!*/
          }                                     /* von maximal 4 comparisions */
        }                         /* if ... != NULL */
      }                   /* else */
    }             /* if touch,... */
  }       /* for-loop */

  theIntersectfoundPoints[*ii] = NULL;

}

/****************** end of function environment_search **********************/







/****************************************************************************/
/*                                                                          */
/* Function:  insert                                                                    */
/*                                                                          */
/* Purpose:   inserts a node into the quadtree                                                          */
/*                                                                          */
/* Input:     QFCLISTTYP *p_new : this FC must be inserted in the quadtree      */
/*			  QUADTREETYP *q_place : the inserting place, determined by the */
/*									 function search()                      */
/*			  SOURCETYP *src : concerning source							*/
/*			  COORD wi : width of the subquadrangles of the concerning      */
/*				         depth												*/
/*																			*/
/****************************************************************************/

static void insert(QFCLISTTYP *p_new, QUADTREETYP *q_place,
                   SOURCETYP *src, COORD wi)
{
  int lauf;
  unsigned char helpc1, helpc2, helpc3, helpc_qu_flg;
  QUADTREETYP *q_pointer;
  QFCLISTTYP *n_zeiger;


  helpc1= q_place->q_flag;
  helpc1 >>= 4;
  n_zeiger = q_place->q_array[helpc1];        /* in wich square ? --> helpc1*/
  if ( n_zeiger == NULL )
  {
    /* in this quadrant there are no nodes !!!*/
    q_place->q_array[helpc1] = p_new;
  }       /* n_zeiger == NULL */

  else if
  ((XC(MYVERTEX(FRONTN(FROC(n_zeiger))))==XC(MYVERTEX(FRONTN(FROC(p_new)))))&&
   (YC(MYVERTEX(FRONTN(FROC(n_zeiger))))==YC(MYVERTEX(FRONTN(FROC(p_new)))))  )
  {
    /* insert at top of list: */
    NXT(p_new) = n_zeiger;
    q_place->q_array[helpc1] = p_new;
  }

  else
  {
    /* there is already (!)another node(!) */
    q_pointer = GetFreeObject( MG, QuObj);
    if ( q_pointer == NULL )
    {
      q_pointer = GetMem(MG->theHeap,sizeof(QUADTREETYP),FROM_TOP);
      if ( q_pointer == NULL )
      {
        PrintErrorMessage('E',"bnodes"," ERROR: No memory !!! error in quadtreefunction <insert>");
      }
    }
    SETOBJT(q_pointer,QuObj);
    for (lauf=0; lauf<=3; lauf++)
    {
      q_pointer->q_array[lauf] = NULL;
    }
    q_pointer->q_flag = 15;             /* d.i. "00001111" */

    q_place->q_array[helpc1] = q_pointer;
    /* Attention! Now the flag, which characterisizes the new quadtreepoiner */
    /* has to be set from 1 to 0 !!!    */
    switch(helpc1)
    {
    case 0 : helpc_qu_flg = 254; break;                    /* "1111_1110" */
    case 1 : helpc_qu_flg = 253; break;                    /* "1111_1101" */
    case 2 : helpc_qu_flg = 251; break;                    /* "1111_1011" */
    case 3 : helpc_qu_flg = 247; break;                    /* "1111_0111" */
    }             /* von switch(helpc1) */
                  /* The following and-operation sets the  concerning bit from "1" */
                  /* (= nodepointer) to "0" (= quadtreepointer) !!! */

    q_place->q_flag &= helpc_qu_flg;

    switch (helpc1)
    {
    case 0 : break;
    case 1 : src->x += wi; break;
    case 2 : src->x += wi; src->y += wi; break;
    case 3 : src->y += wi;
    }

    wi = wi/2;              /* width-adaptation !!! */
    helpc2 = 0;

    if ( (YC(MYVERTEX(FRONTN(FROC(n_zeiger)))) < (src->y + wi)) &&
         (XC(MYVERTEX(FRONTN(FROC(n_zeiger)))) >= (src->x + wi))   )
      helpc2 = 1;

    else if (YC(MYVERTEX(FRONTN(FROC(n_zeiger)))) >= (src->y + wi))
    {
      if (XC(MYVERTEX(FRONTN(FROC(n_zeiger)))) >= (src->x + wi))
        helpc2 = 2;
      else
        helpc2 = 3;
    }

    /* for all quadrants together:    */
    helpc2 <<= 4;
    q_pointer->q_flag |= helpc2;
    helpc2 >>= 4;

    insert(n_zeiger, q_pointer, src, wi);

    /* Besides we must consider, whether we can insert the point p_new or have
       /* to refine again */
    helpc3 = 0;
    if ( (YC(MYVERTEX(FRONTN(FROC(p_new)))) < (src->y + wi)) &&
         (XC(MYVERTEX(FRONTN(FROC(p_new)))) >= (src->x + wi))   )
      helpc3 = 1;

    else if (YC(MYVERTEX(FRONTN(FROC(p_new)))) >= (src->y + wi))
    {
      if (XC(MYVERTEX(FRONTN(FROC(p_new)))) >= (src->x + wi))
        helpc3 = 2;
      else
        helpc3 = 3;
    }

    /* for all quadrants together     */
    helpc3 <<= 4;
    q_pointer->q_flag <<= 4;
    q_pointer->q_flag >>= 4;
    q_pointer->q_flag |= helpc3;

    insert(p_new,q_pointer,src,wi);

  }       /* else */
}

/************************* end of function insert ***************************/






/****************************************************************************/
/*                                                                          */
/* Function:  delete_node                                                               */
/*                                                                          */
/* Purpose:   deletes a frontcomponent within the quadtree                                      */
/*                                                                          */
/* Input:     QUADTREETYP *q_pointer : pointer at the concerning                        */
/*									   structure							*/
/*			  FRONTCOMP *p_del : this frontcomponent shall be deleted               */
/*			  COORD wi : width of the subquadrangles of the concerning      */
/*				         depth												*/
/*			  SOURCETYP *so : concerning source								*/
/*			  int *stop, QFCLISTTYP **nd_mem : necessary for deleting               */
/*													process                 */
/*																			*/
/****************************************************************************/

static void delete_node(QUADTREETYP *q_pointer, FRONTCOMP *p_del, COORD width,
                        SOURCETYP *so, int *stop, QFCLISTTYP **nd_mem)
{
  unsigned char helpc, hc;
  QFCLISTTYP *nodepointer_qfcl, *onebefore;
  int place, zeroes, i;

  if  ( YC(MYVERTEX(FRONTN(p_del))) < ( so->y + width) )
  {       /* the node must be in 3. or 2.quadrant */
    if ( XC(MYVERTEX(FRONTN(p_del))) < (so->x + width) )
    {
      /* 3. quadrant */
      place = 0;
    }
    else
    {
      /* 2. quadrant */
      place = 1;
      so->x += width;
    }
  }       /* of  3. or 2.quadrant */

  else if ( XC(MYVERTEX(FRONTN(p_del))) >= (so->x + width) )
  {
    /* 1. quadrant */
    place = 2;
    so->x+=width;
    so->y+=width;
  }

  else
  {
    place = 3;
    so->y+=width;
  }       /* 4.quadrant*/

  helpc = q_pointer->q_flag; ; helpc >>= place; helpc <<= 7;

  if (helpc == 0)  /*  here is a quadtreepointer !!! */
  {
    delete_node( q_pointer->q_array[place] , p_del, width/2, so, stop, nd_mem);
  }

  else  /* the node should be here */
  {
    nodepointer_qfcl = q_pointer->q_array[place];
    if ( nodepointer_qfcl == NULL )             /* there is no node at all !!! */
    {
      PrintErrorMessage('E',"bnodes","Error: I cannot delete a node, which  doesn«t exist!!!");
    }
    else
    {               /* here we are !!! */
      if ( FROC(nodepointer_qfcl) == p_del )
      {
        q_pointer->q_array[place] = NXT(nodepointer_qfcl);
        PutFreeObject(MG,nodepointer_qfcl);
      }
      else
      {
        do
        {
          if ( NXT(nodepointer_qfcl) != NULL )
          {
            onebefore = nodepointer_qfcl;
            nodepointer_qfcl = NXT(nodepointer_qfcl);
          }
          else
            PrintErrorMessage('E',"bnodes","ERR: in delete_node QFCL: node doesn«t exist !");
        }
        while ( FROC(nodepointer_qfcl) != p_del );

        NXT(onebefore) = NXT(nodepointer_qfcl);
        PutFreeObject(MG, nodepointer_qfcl);
      }

    }

  }




  if (*nd_mem != NULL)
  {
    /* flag and pointer of the tree are adapted!!! */
    switch (place)
    {
    case 0 : q_pointer->q_flag |= 1; break;
    case 1 : q_pointer->q_flag |= 2; break;
    case 2 : q_pointer->q_flag |= 4; break;
    case 3 : q_pointer->q_flag |= 8; break;
    }
    q_pointer->q_array[place] = *nd_mem;
    *nd_mem = NULL;
  }

  if ( *stop ==  0 )
  {
    /* Attention! Do we now have a quadtreepointer with only one QFCL !? */
    zeroes = 0;
    for (i=0; i<=3; i++)
    {
      if ( q_pointer->q_array[i] == NULL )
        zeroes++;
      else
        place = i;
    }

    hc = q_pointer->q_flag;
    hc >>= place; hc <<= 7;


    if ( ( zeroes == 3 ) && ( hc > 0 ) && ( width != startwidth/2 )  )
    /* 1) The critical case happened, that there are exactly three        */
    /*   nullpointer belonging to a quadtreepointer,                      */
    /* 2) but the fourth pointer is no quadtreepoiner                     */
    /* 3) and we are not at the startlevel !!!                            */
    {

      nodepointer_qfcl = q_pointer->q_array[place];
      *nd_mem = nodepointer_qfcl;

      PutFreeObject( MG, q_pointer);

    }
    else *stop = 1;
  }               /* von stop */
}

/********************** end of function delete_node *************************/









/****************************************************************************/
/*                                                                          */
/* Function:  delete_quadtree                                                   */
/*                                                                          */
/* Purpose:   deletes the  complete quadtree                                                            */
/*                                                                          */
/* Input:     QUADTREETYP *q_pointer : pointer at the concerning quadtree       */
/*									   structure, beginning with the root	*/
/*                                                                          */
/****************************************************************************/

static void delete_quadtree(QUADTREETYP *q_pointer)
{
  unsigned char helpc;
  int i;


  for(i=0; i<4; i++)
  {
    helpc = q_pointer->q_flag;
    helpc >>= i;
    helpc <<= 7;

    if (helpc == 0)
    {
      delete_quadtree(q_pointer->q_array[i]);
    }
  }

  PutFreeObject( MG, q_pointer);
}

/******************** end of function delete_quadtree ***********************/









/****************************************************************************/
/*                                                                          */
/* Function:  FCTreeUpdate                                                              */
/*                                                                          */
/* Purpose:   distinguishs the calls of DELETE_ND and InsertQuadtree            */
/*                                                                          */
/* Input:     FRONTCOMP *fc : concerning frontcomponent                                         */
/*			  int n : flag for distinction									*/
/*                                                                          */
/****************************************************************************/

static void FCTreeUpdate(FRONTCOMP *fc, int n)
{
  if ( n < 1 )
    DELETE_ND(fc);
  else
    InsertQuadtree(fc,n);
}

/******************** end of function FCTreeUpdate **************************/

















/****************************************************************************/
/*                                                                          */
/* Function:  InsertQuadtree                                                    */
/*                                                                          */
/* Purpose:   prepairs quadtree inserting process and distinguishs the          */
/*			  number of frontcomponents to be inserted						*/
/*                                                                          */
/* Input:     FRONTCOMP *pon : pointer at the array of frontcomponents      */
/*            int ncomp: number of frontcomponents which must be inserted   */
/*                                                                          */
/****************************************************************************/

static void InsertQuadtree(FRONTCOMP *pon, int ncomp)
{
  SOURCETYP *srce;
  QFCLISTTYP *p_new;
  int i;
  QUADTREETYP *qz_s;
  COORD actual_width;


  srce = GetFreeObject( MG, ScObj);
  if ( srce == NULL )
  {
    srce = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_TOP);
    if ( srce == NULL )
    {PrintErrorMessage('E',"bnodes","ERROR: No memory !!! in InsertQuadtree");}
  }
  else
    i=1;

  SETOBJT(srce,ScObj);

  for ( i=0; i<ncomp; i++)
  {
    srce->x = source->x; srce->y = source->y;
    actual_width = startwidth/2;
    qz_s = search( startpointer, srce, &actual_width, &pon[i]);

    p_new = GetFreeObject( MG, QfclObj);

    if ( p_new == NULL )
    {
      p_new = GetMem(MG->theHeap,sizeof(QFCLISTTYP),FROM_TOP);
      if ( p_new == NULL )
      {
        PrintErrorMessage('E',"bnodes","ERR:No memory! -> quadtreefunction <InsertQuadtree>");
      }
    }

    SETOBJT(p_new,QfclObj);
    FROC(p_new) = &pon[i];
    NXT(p_new) = NULL;

    insert(p_new, qz_s, srce, actual_width);
  }
  PutFreeObject(MG, srce);
}

/******************* end of function InsertQuadtree *************************/





/****************************************************************************/
/*                                                                          */
/* Function:  DELETE_ND                                                                 */
/*                                                                          */
/* Purpose:   prepairs quadtree deleting process                                                */
/*                                                                          */
/* Input:     FRONTCOMP *delete_p : pointer at the frontcomponent which     */
/*									will be deleted							*/
/*                                                                          */
/****************************************************************************/

static void DELETE_ND ( FRONTCOMP *delete_p )
{
  int aufhoeren;
  SOURCETYP *srce;
  QFCLISTTYP *nd_mem;


  srce = GetFreeObject( MG, ScObj);
  if ( srce == NULL )
  {
    srce = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_TOP);
    if ( srce == NULL )
    {PrintErrorMessage('E',"bnodes","ERROR: No memory !!! in InsertQuadtree");}
  }
  SETOBJT(srce, ScObj);

  aufhoeren = 0;
  srce->x = source->x; srce->y = source->y;
  nd_mem = NULL;
  delete_node( startpointer, delete_p, startwidth/2, srce, &aufhoeren, &nd_mem);
  PutFreeObject(MG, srce);

}

/********************** end of function DELETE_ND ***************************/





/****************************************************************************/
/*                                                                          */
/* Function:  btree_ins                                                                 */
/*                                                                          */
/* Purpose:   inserts a frontcomponent in a balanced tree dependent on the  */
/*			  the concerning angle size respectively edge length            */
/*                                                                          */
/* Input:     FRONTCOMP *basefc : frontcomponent, which will be inserted        */
/*			  float x : appertaining angle size respectively edge length    */
/*			  BALTREETYP **p, int *h : necessary for balancing              */
/*                                                                          */
/****************************************************************************/

static void btree_ins(FRONTCOMP *basefc, float x, BALTREETYP **p, int *h)
{
  BALTREETYP *p1,*p2;

  if ( *p == NULL )
  {
    *p = GetFreeObject( MG, EttObj);
    if ( *p == NULL )
    {
      *p = GetMem(MG->theHeap,sizeof(BALTREETYP),FROM_TOP);
      if ( *p == NULL )
      {PrintErrorMessage('E',"bnodes"," ERROR: No memory !!! in btree_ins");}
    }
    SETOBJT(*p, EttObj);


    *h = 1;
    (*p)->eFC = basefc;
    (*p)->key = x;
    (*p)->left = NULL;
    (*p)->right = NULL;
    (*p)->bal = 0;
  }

  else if ( x <= (*p)->key)
  {
    btree_ins(basefc, x, &((*p)->left), h);
    if (*h == 1 )             /* left branch became longer */
      switch ( (*p)->bal )
      {
      case  1 : (*p)->bal = 0;*h = 0; break;
      case  0 : (*p)->bal = -1; break;
      case -1 :                  /* renewed balancing is necessary */

        p1 = (*p)->left;

        if ( p1->bal == -1 )
        {
          /* simple LL-rotation */
          (*p)->left = p1->right;
          p1->right = *p;
          (*p)->bal = 0;  *p = p1;
        }

        else
        {
          /* double LR-rotation */
          p2 = p1->right;
          p1->right = p2->left; p2->left = p1;
          (*p)->left = p2->right; p2->right = *p;
          if ( p2->bal == -1 ) (*p)->bal =   1;
          else (*p)->bal =   0;
          if ( p2->bal ==  1 ) p1->bal = -1;
          else p1->bal =  0;
          *p = p2;
        }
        (*p)->bal = 0; *h = 0;

        break;
      }           /* of switch */

  }       /* of else if ( x <= (*p)->key ) */

  else if ( x > (*p)->key )
  {
    btree_ins(basefc, x, &((*p)->right), h);
    if (*h == 1 )             /* left branch became longer */
      switch ( (*p)->bal )
      {
      case -1 : (*p)->bal = 0;*h = 0; break;
      case  0 : (*p)->bal = +1; break;
      case  1 :                  /* renewed balancing is necessary */

        p1 = (*p)->right;

        if ( p1->bal == +1 )
        {
          /* simple RR-rotation */
          (*p)->right = p1->left;
          p1->left = *p;
          (*p)->bal = 0;  *p = p1;
        }
        else
        {
          /* double RL-rotation */
          p2 = p1->left;
          p1->left = p2->right; p2->right = p1;
          (*p)->right = p2->left; p2->left = *p;
          if ( p2->bal == +1 ) (*p)->bal =  -1;
          else (*p)->bal =   0;
          if ( p2->bal == -1 ) p1->bal = +1;
          else p1->bal =  0;
          *p = p2;
        }

        (*p)->bal = 0; *h = 0;
        break;
      }           /* of switch */

  }       /* of else if ( x > *p->key ) */

}

/********************** end of function btree_ins ***************************/




/****************************************************************************/
/*                                                                          */
/* Function:  balance1                                                                  */
/*                                                                          */
/* Purpose:   necessary for balancing basetree								*/
/*            in particular RR and RL-rotation                                                  */
/*																			*/
/****************************************************************************/

static void balance1( BALTREETYP **p, int *h )
{
  BALTREETYP *p1, *p2;
  int b1,b2;

  switch ( (*p)->bal )
  {
  case -1 :        (*p)->bal = 0; break;
  case  0 :        (*p)->bal = +1; *h = 0; break;
  case  1 :                     /* renewed balancing is necessary */

    p1 = (*p)->right;
    b1 = p1->bal;

    if ( b1 >= 0 )
    {
      /* simple RR-rotation */
      (*p)->right = p1->left;
      p1->left = *p;
      if ( b1 == 0 )
      {
        (*p)->bal = 1;
        p1->bal = -1;
        *h = 0;
      }
      else
      {
        (*p)->bal = 0;
        p1->bal = 0;                                                 /*!!! 0 instead of -1 in Wirths book */
      }

      *p = p1;
    }

    else
    {
      /* double RL-rotation */
      p2 = p1->left;
      b2 = p2->bal;
      p1->left = p2->right;
      p2->right = p1;
      (*p)->right = p2->left;
      p2->left = *p;
      if ( b2 == +1 )
        (*p)->bal =  -1;
      else
        (*p)->bal =   0;
      if ( b2 == -1 )
        p1->bal = +1;
      else
        p1->bal =  0;
      *p = p2;
      p2->bal = 0;
    }

    break;
  }       /* of switch */
}

/********************** end of function balance1 ****************************/




/****************************************************************************/
/*                                                                          */
/* Function:  balance2                                                                  */
/*                                                                          */
/* Purpose:   necessary for balancing basetree								*/
/*            in particular LR and LL-rotation                                                  */
/*																			*/
/****************************************************************************/

static void balance2( BALTREETYP **p, int *h )
{
  BALTREETYP *p1, *p2;
  int b1, b2;

  switch ( (*p)->bal )
  {
  case  1 :        (*p)->bal = 0; break;
  case  0 :        (*p)->bal = -1; *h = 0; break;
  case -1 :                     /* renewed balancing is necessary */

    p1 = (*p)->left;
    b1 = p1->bal;

    if ( b1 <= 0 )
    {
      /* simple LL-rotation */
      (*p)->left = p1->right;
      p1->right = *p;
      if ( b1 == 0 )
      {
        (*p)->bal = -1;
        p1->bal = +1;
        *h = 0;
      }
      else
      {
        (*p)->bal = 0;
        p1->bal = 0;
      }
      *p = p1;
    }
    else
    {
      /* double LR-rotation */
      p2 = p1->right;
      b2 = p2->bal;
      p1->right = p2->left;
      p2->left = p1;
      (*p)->left = p2->right;
      p2->right = *p;
      if ( b2 == -1 )
        (*p)->bal =  +1;
      else
        (*p)->bal =   0;

      if ( b2 == +1 )
        p1->bal = -1;
      else
        p1->bal =  0;

      *p = p2;
      p2->bal = 0;
    }

    break;
  }       /* of switch */
}

/********************** end of function balance2 ****************************/



/****************************************************************************/
/*                                                                          */
/* Function:  del                                                                               */
/*                                                                          */
/* Purpose:   necessary for deleting process in basetree					*/
/*																			*/
/****************************************************************************/

static void del( BALTREETYP **r, int *h )
{
  if ( (*r)->right != NULL )
  {
    del( &((*r)->right), h );
    if ( *h == 1 ) balance2( r, h);
  }
  else
  {
    (*q)->key = (*r)->key;
    (*q)->eFC = (*r)->eFC;
    *r = (*r)->left; *h = 1;
  }
}

/************************* end of function del ******************************/










/****************************************************************************/
/*                                                                          */
/* Function:  delete                                                                    */
/*                                                                          */
/* Purpose:   deletes a frontcomponent in a balanced tree dependent on the  */
/*			  the concerning angle size respectively edge length            */
/*                                                                          */
/* Input:     FRONTCOMP *basefc : frontcomponent, which will be deleted         */
/*			  float x : appertaining angle size respectively edge length    */
/*			  BALTREETYP **p, int *h : necessary for balancing              */
/*                                                                          */
/****************************************************************************/

static void delete( FRONTCOMP* basefc, float x, BALTREETYP **p, int *h )
{
  if ( ( x < (*p)->key ) && ( (*p)->left != NULL ) )
  {
    delete( basefc, x, &((*p)->left), h );
    if ( del_edg_fnd == 1 )
      if ( *h == 1 )
        balance1( p, h );
  }

  else if ( ( x > (*p)->key ) && ( (*p)->right != NULL ) )
  {
    delete( basefc, x, &((*p)->right), h );
    if ( del_edg_fnd == 1 )
      if ( *h == 1 )
        balance2( p, h );
  }

  else if ( ( (*p)->key == x ) && ( (*p)->eFC != basefc ) )
  {
    if ( ( del_edg_fnd == 0  ) && ( (*p)->left != NULL ) )
    {
      delete( basefc, x, &((*p)->left), h );
      if ( del_edg_fnd == 1 )
        if ( *h == 1 )
          balance1( p, h );
    }

    if ( ( del_edg_fnd == 0  ) && ( (*p)->right != NULL ) )
    {
      delete( basefc, x, &((*p)->right), h );
      if ( del_edg_fnd == 1 )
        if ( *h == 1 )
          balance2( p, h );
    }
  }


  else if ( (*p)->key == x )
  {
    del_edg_fnd = 1;

    q = p;

    if ( (*q)->right == NULL )
    {
      if      ( (*q)->left == NULL )
      {
        *p = NULL; *h = 1;
      }
      else
      {
        *p = (*q)->left;
        *h = 1;
      }
    }

    else if ( (*q)->left == NULL )
    {
      *p = (*q)->right; *h = 1;
    }

    else
    {
      del( &((*q)->left), h );
      if ( *h == 1 )
        balance1( p, h );
    }

  }


}

/************************* end of function delete ******************************/




/****************************************************************************/
/*                                                                          */
/* Function:  BaseTreeUpdate                                                    */
/*                                                                          */
/* Purpose:   prepairs updating of the basetree and distinguishs the        */
/*			  base criterion and besides the kind of tree manipulation      */
/*                                                                          */
/* Input:	  FRONTCOMP* P : predecessor of the concerning frontcomponent Q */
/*            FRONTCOMP* Q : The frontcomponent which will be inserted or   */
/*							 deleted within the balanced tree               */
/*			  FRONTCOMP* S : successor of the concerning frontcomponent Q   */
/*			  int ch       : flag to distinguish deleting and inserting     */
/*                                                                          */
/****************************************************************************/

static void BaseTreeUpdate( FRONTCOMP* P, FRONTCOMP* Q, FRONTCOMP* S, int ch,
                            int anglecrit, int edgecrit)
{
  int h;
  float edge_lth;
  float angle_wth;

  if (anglecrit)
  {
    angle_wth = length_of_angle(P, Q, S);
    h = 0;

    if ( ch == 0 )
    {
      del_edg_fnd = 0;
      delete( Q, angle_wth, &btree_rootpointer, &h );
      if ( del_edg_fnd == 0 )
      {PrintErrorMessage('E',"bnodes","ERROR: node not found in Edgetree");}
    }
    else
      btree_ins(Q, angle_wth, &btree_rootpointer, &h);
  }

  else if (edgecrit)
  {
    edge_lth = length_of_edge(Q, S);
    h = 0;

    if ( ch == 0 )
    {
      del_edg_fnd = 0;
      delete( Q, edge_lth, &btree_rootpointer, &h );
      if ( del_edg_fnd == 0 )
      {PrintErrorMessage('E',"bnodes","ERROR: node not found in Edgetree");}
    }
    else
      btree_ins(Q, edge_lth, &btree_rootpointer, &h);
  }
}

/*********************** end of function BaseTreeUpdate ************************/







/****************************************************************************/
/*                                                                          */
/* Function:  AccelUpdate                                                               */
/*                                                                          */
/* Purpose:   prepairs updating of the quadtree and the basetree            */
/*			  distinguishs the different cases of new elemenmts             */
/*                                                                          */
/* Input:	  FRONTCOMP* theFC,thenewFC,the_old_succ : the concerning FC«s  */
/*            int cas : flag for distinction of the different cases             */
/*			  int anglecrit,edgecrit : flag for distinction of edge and     */
/*									   angle criterion                      */
/*                                                                          */
/****************************************************************************/

void AccelUpdate( FRONTCOMP* theFC,  FRONTCOMP* thenewFC, FRONTCOMP* the_old_succ, int cas,  int anglecrit,  int edgecrit )
{
  switch (cas)
  {

  case NORMALCASE :

    BaseTreeUpdate(theFC, SUCCFC(theFC), SUCCFC(thenewFC),
                   1, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(theFC), theFC, SUCCFC(thenewFC),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(theFC), theFC, SUCCFC(theFC),
                   1, anglecrit, edgecrit);

    BaseTreeUpdate(theFC, SUCCFC(thenewFC), SUCCFC(SUCCFC(thenewFC)),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(thenewFC, SUCCFC(thenewFC), SUCCFC(SUCCFC(thenewFC)),
                   1, anglecrit, edgecrit);

    FCTreeUpdate(thenewFC, 1);

    break;


  case LEFTNEIGHBOURCASE :

    BaseTreeUpdate(thenewFC, theFC, the_old_succ,
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(thenewFC), thenewFC, theFC,
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(thenewFC), thenewFC, the_old_succ,
                   1, anglecrit, edgecrit);

    BaseTreeUpdate(theFC, the_old_succ, SUCCFC(SUCCFC(thenewFC)),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(thenewFC, the_old_succ, SUCCFC(SUCCFC(thenewFC)),
                   1, anglecrit, edgecrit);

    FCTreeUpdate(theFC,0);

    break;



  case RIGHTNEIGHBOURCASE :
    BaseTreeUpdate(theFC, the_old_succ, thenewFC,
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(the_old_succ, thenewFC, SUCCFC(thenewFC),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(theFC, thenewFC, SUCCFC(thenewFC),
                   1, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(theFC), theFC, the_old_succ,
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(theFC), theFC, thenewFC,
                   1, anglecrit, edgecrit);

    FCTreeUpdate(the_old_succ,0);

    break;


  case ININTERCASE :
    BaseTreeUpdate(PREDFC(thenewFC), thenewFC, SUCCFC(SUCCFC(theFC)),
                   0, anglecrit, edgecrit);

    FCTreeUpdate(SUCCFC(theFC), 1);

    BaseTreeUpdate(theFC, SUCCFC(theFC), SUCCFC(SUCCFC(theFC)),
                   1, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(thenewFC), thenewFC, SUCCFC(thenewFC),
                   1, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(theFC), theFC, SUCCFC(thenewFC),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(PREDFC(theFC), theFC, SUCCFC(theFC),
                   1, anglecrit, edgecrit);

    BaseTreeUpdate(theFC, SUCCFC(thenewFC), SUCCFC(SUCCFC(thenewFC)),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(thenewFC, SUCCFC(thenewFC), SUCCFC(SUCCFC(thenewFC)),
                   1, anglecrit, edgecrit);

    break;


  case FINALCASE :
    BaseTreeUpdate(PREDFC(theFC), theFC, SUCCFC(theFC),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(theFC, SUCCFC(theFC), SUCCFC(SUCCFC(theFC)),
                   0, anglecrit, edgecrit);

    BaseTreeUpdate(SUCCFC(theFC), SUCCFC(SUCCFC(theFC)), theFC,
                   0, anglecrit, edgecrit);

    FCTreeUpdate(theFC,0);

    FCTreeUpdate(SUCCFC(theFC),0);

    FCTreeUpdate(SUCCFC(SUCCFC(theFC)),0);

    break;


  default : PrintErrorMessage('E',"bnodes"," ERROR: This case is not allowed! <AccelUpdate>");
  }
}

/************************ end of function AccelUpdate **************************/







/****************************************************************************/
/*                                                                          */
/* Function:  AccelInit                                                                 */
/*                                                                          */
/* Purpose:   Init function of the accelerator                              */
/*                                                                          */
/* Input:	  GRID *the_Grid : pointer at the Grid                          */
/*            int anglecrit, edgecrit : flag for distinction of edge crit.  */
/*										and angle criterion                 */
/*                                                                          */
/****************************************************************************/

int AccelInit(GRID *the_Grid, int anglecrit, int edgecrit, GG_PARAM *params)
{
  int l;
  DOMAIN *ilDomeno;
  INDEPFRONTLIST *edge_theIFL;
  FRONTLIST *edge_theFL;
  FRONTCOMP *edge_theFC;



  myPars = params;

  /* Adaption to the Multigrid: */
  MG = GetCurrentMultigrid();
  if (MG == NULL) PrintErrorMessage('E',"bnodes","no multigrid received");

  InitAccelObjs(MG);


  ilDomeno = MG->theDomain;


  del_edg_fnd = 0;



  startpointer = GetFreeObject( MG, QuObj);
  if ( startpointer == NULL )
  {
    startpointer = GetMem(MG->theHeap,sizeof(QUADTREETYP),FROM_TOP);
    if ( startpointer == NULL )
    {PrintErrorMessage('E',"bnodes","ERROR: No memory !!!");}
  }

  SETOBJT(startpointer,QuObj);
  startpointer->q_flag = 15;       /* ==  "00001111" !!!! */
  for ( l=0; l<=3; l++ )
  {
    startpointer->q_array[l] = NULL;
  }

  source = GetFreeObject( MG, ScObj);
  if ( source == NULL )
  {
    source = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_TOP);
    if ( source == NULL )
    {PrintErrorMessage('E',"bnodes","ERROR: No memory !!!");return(1);}
  }
  SETOBJT(source,ScObj);
  source->x = ilDomeno->MidPoint[0] - ilDomeno->radius;
  source->y = ilDomeno->MidPoint[1] - ilDomeno->radius;

  startwidth = ilDomeno->radius + ilDomeno->radius;

  /* for the edge_tree_initialzation */

  btree_rootpointer = NULL;

  myMGdata = GetMGdataPointer(MYMG(the_Grid));
  for (edge_theIFL=STARTIFL(myMGdata); edge_theIFL!=NULL;
       edge_theIFL=SUCCIFL(edge_theIFL))
  {
    for (edge_theFL=STARTFL(edge_theIFL); edge_theFL!=NULL;
         edge_theFL=SUCCFL(edge_theFL))
    {
      for (edge_theFC=STARTFC(edge_theFL); edge_theFC!=NULL;
           edge_theFC=SUCCFC(edge_theFC))
      {
        BaseTreeUpdate( PREDFC(edge_theFC), edge_theFC, SUCCFC(edge_theFC),
                        1, anglecrit, edgecrit);
        FCTreeUpdate(edge_theFC,1);
        if ( edge_theFC == LASTFC(edge_theFL) )
          break;
      }
    }
  }

  return(0);
}

/************************* end of function AccelInit ***************************/













/****************************************************************************/
/*                                                                          */
/* Function:  AccelFCTreeSearch                                                 */
/*                                                                          */
/* Purpose:   prepairs the search for frontcomponents wich ar close or      */
/*            within the suggested triangle									*/
/*                                                                          */
/* Input:     INDEPFRONTLIST *theIFL: the Independent Front List, the new   */
/*									  suggested triangle belongs to         */
/*			  FRONTCOMP* thefoundPoints[MAXNPOINTS] : Array for the FC«s,   */
/*			                          found within the small searching      */
/*									  rectangle								*/
/*			  FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS] :  Array for   */
/*									  the FC«s, found within the big        */
/*									  searching rectangle - necessary to    */
/*									  detect long cutting edges             */
/*			  COORD xt[3], COORD yt[3] : coordinates of the sugg. triangle	*/
/*			  COORD searchradis : radius of the circle around the new       */
/*									  suggested node to detect also close   */
/*									  cutting edges							*/
/*																			*/
/* Output:    returns the number of the interesting frontcomponents in the  */
/*			  small searching rectangle										*/
/*                                                                          */
/****************************************************************************/

int AccelFCTreeSearch(INDEPFRONTLIST *theIFL, FRONTCOMP* thefoundPoints[MAXNPOINTS],
                      FRONTCOMP *theIntersectfoundPoints[MAXNPOINTS], COORD xt[3],
                      COORD yt[3], COORD searchradis)
{
  SOURCETYP *srce;
  SOURCETYP *search_sq_ld, *search_sq_ru;
  SOURCETYP *big_search_sq_ld, *big_search_sq_ru;
  float maxsidelength;
  int foundpoints, ii;


  srce = GetFreeObject( MG, ScObj);
  if ( srce == NULL )
  {
    srce = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_TOP);
    if ( srce == NULL )
    {PrintErrorMessage('E',"bnodes","ERROR: No memory !!! in InsertQuadtree");}
  }
  SETOBJT(srce, ScObj);
  srce->x = source->x; srce->y = source->y;

  search_sq_ru = GetFreeObject( MG, ScObj);
  if ( search_sq_ru == NULL )
  {
    search_sq_ru = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_TOP);
    if ( search_sq_ru == NULL )
    {PrintErrorMessage('E',"bnodes","No memory !!! in InsertQuadtree");}
  }
  SETOBJT(search_sq_ru, ScObj);

  search_sq_ld = GetFreeObject( MG, ScObj);
  if ( search_sq_ld == NULL )
  {
    search_sq_ld = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_TOP);
    if ( search_sq_ld == NULL )
    {PrintErrorMessage('E',"bnodes","ERROR: No memory !!! in InsertQuadtree");}
  }
  SETOBJT(search_sq_ld, ScObj);


  big_search_sq_ru = GetFreeObject( MG, ScObj);
  if ( big_search_sq_ru == NULL )
  {
    big_search_sq_ru = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_BOTTOM);
    if ( big_search_sq_ru == NULL )
    {PrintErrorMessage('E',"bnodes","ERROR: No memory !!! in InsertQuadtree");}
  }
  SETOBJT(big_search_sq_ru, ScObj);

  big_search_sq_ld = GetFreeObject( MG, ScObj);
  if ( big_search_sq_ld == NULL )
  {
    big_search_sq_ld = GetMem(MG->theHeap,sizeof(SOURCETYP),FROM_BOTTOM);
    if ( big_search_sq_ld == NULL )
    {PrintErrorMessage('E',"bnodes"," ERROR: No memory !!! in InsertQuadtree");}
  }
  SETOBJT(big_search_sq_ld, ScObj);



  /* Now we build the searching rectangle !!! */
  search_sq_ld->x = min ( xt[2] - searchradis, xt[0]);
  search_sq_ld->x = min ( search_sq_ld->x, xt[1]);
  search_sq_ld->y = min ( yt[2] - searchradis, yt[0]);
  search_sq_ld->y = min ( search_sq_ld->y, yt[1]);
  search_sq_ru->x = max ( xt[2] + searchradis, xt[0]);
  search_sq_ru->x = max ( search_sq_ru->x, xt[1]);
  search_sq_ru->y = max ( yt[2] + searchradis, yt[0]);
  search_sq_ru->y = max ( search_sq_ru->y, yt[1]);

  maxsidelength = myPars->h_global;

  big_search_sq_ld->x = search_sq_ld->x - maxsidelength;
  big_search_sq_ld->y = search_sq_ld->y - maxsidelength;
  big_search_sq_ru->x = search_sq_ru->x + maxsidelength;
  big_search_sq_ru->y = search_sq_ru->y + maxsidelength;


  /* Thus the key-hole is still considered */



  foundpoints = 0; ii = 0;
  environment_search( theIFL, startpointer, srce, thefoundPoints,
                      theIntersectfoundPoints, startwidth/2, search_sq_ld,
                      search_sq_ru, big_search_sq_ld, big_search_sq_ru, xt, yt,
                      searchradis, &foundpoints, &ii );
  PutFreeObject(MG, srce);
  PutFreeObject(MG, search_sq_ru);
  PutFreeObject(MG, search_sq_ld);

  return (foundpoints);

}

/********************* end of function AccelFCTreeSearch ***********************/






/****************************************************************************/
/*                                                                          */
/* Function:  AccelBaseTreeSearch                                               */
/*                                                                          */
/* Purpose:   delivers the FC with the smallest edge respectively the FC    */
/*            with the smallest angle as the base node for the next tiangle */
/*                                                                          */
/* Output:    FRONTLIST *myList: reference parameter for the pointer at the */
/*                               Frontlist of the found FC                      */
/*                                                                          */
/* Return:    the "ideal" FC for the next basis during the triangulation    */
/*                                                                          */
/****************************************************************************/

FRONTCOMP* AccelBaseTreeSearch(FRONTLIST** myList)
{
  BALTREETYP* p;
  p = btree_rootpointer;
  if ( p == NULL )
  {
    return(NULL);
  }
  else
  {
    while ( p->left != NULL )
      p = p->left;
    *myList = p->eFC->myFL;
    return (p->eFC);
  }
} /* of AccelBaseTreeSearch() */

/******************* end of function AccelBaseTreeSearch **********************/
