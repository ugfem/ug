// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdio.h>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>

#include <template.hh>
#include <array.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>
#include <meshing/ruler2.hh>

static double maxdouble = 1.0e37;
static double SMALLDOUBLE = 1.0e-12;

static int CalcNewPoint(double xt[3],double yt[3],double h)
{
  double xp1p2,yp1p2,lp1p2;

  xp1p2 = xt[0] -xt[1];
  yp1p2 = yt[0] -yt[1];
  lp1p2 = sqrt(xp1p2*xp1p2+yp1p2*yp1p2);

  xt[2] = 0.5 * (xt[0] + xt[1]) + yp1p2/lp1p2/1.0;
  yt[2] = 0.5 * (yt[0] + yt[1]) - xp1p2/lp1p2/1.0;

  return(0);
}

static int CheckNewElement(ARRAY<Point2d> &lpoints, ARRAY<ILINE> &llines,double xt[3],double yt[3])
{
  int i;
  double xM,yM,xR,yR,xQ,yQ,denominator,lambda1,lambda2;

  xM = 0.5 * (xt[0] + xt[1]);
  yM = 0.5 * (yt[0] + yt[1]);

  for(i=2; i<=llines.Size(); i++)
  {
    xR = lpoints[llines[i].I1()].X();
    yR = lpoints[llines[i].I1()].Y();
    xQ = lpoints[llines[i].I2()].X();
    yQ = lpoints[llines[i].I2()].Y();

    denominator = (yR-yQ)*(xt[2]-xM) - (xR-xQ)*(yt[2]-yM);

    if(fabs(denominator)<SMALLDOUBLE)
      continue;
    if(fabs(xt[2]-xM)<SMALLDOUBLE)
    {
      lambda2 = (xR-xM)/(xR-xQ);
      lambda1 = (yR - yM - lambda2 * (yR - yQ)) / (yt[2] - yM);
    }
    else
    {
      lambda2 = ((yR - yM) * (xt[2] - xM) - (xR - xM) * (yt[2] - yM)) / denominator;
      lambda1 = (xR - xM - lambda2 * (xR - xQ)) / (xt[2] - xM);
    }
    if( (lambda1<1.0-2*SMALLDOUBLE) && (lambda1>2*SMALLDOUBLE) && (lambda2<1.0-2*SMALLDOUBLE) && (lambda2>2*SMALLDOUBLE) )
      return(0);
  }
  return(1);
}

static int Cross_Point(double p0[2], double n0[2], double p1[2], double n1[2], double *cp)
{
  double lambda0, a;

  a = - n0[0] * n1[1] + n0[1] * n1[0];

  lambda0 = ( n1[1] * (p0[0] - p1[0]) - n1[0] * (p0[1] - p1[1]) ) / a;
  /*	lambda1 = ( n0[1] * (p0[0] - p1[0]) - n0[0] * (p0[1] - p1[1]) ) / a;

          cp0[0] = p0[0] + lambda0 * n0[0];
          cp0[1] = p0[1] + lambda0 * n0[1];

          cp1[0] = p1[0] + lambda1 * n1[0];
          cp1[1] = p1[1] + lambda1 * n1[1];

          if(sqrt( (cp0[0]-cp1[0])*(cp0[0]-cp1[0]) + (cp0[1]-cp1[1])*(cp0[1]-cp1[1]) )>SMALLDOUBLE)
                  printf("%s\n", "ERROR");*/

  cp[0] = p0[0] + lambda0 * n0[0];
  cp[1] = p0[1] + lambda0 * n0[1];
  return(0);
}


double Calc_Circumcircusmidpoint(double xt[3],double yt[3],double *midp)
{
  double p0[2],p1[2],p2[2],n0[2],n1[2],n2[2],cp01[2],cp12[2],cp20[2],l;

  /* edgemidpoints */
  p0[0] = 0.5 * ( xt[0] + xt[1] );
  p0[1] = 0.5 * ( yt[0] + yt[1] );

  p1[0] = 0.5 * ( xt[1] + xt[2] );
  p1[1] = 0.5 * ( yt[1] + yt[2] );

  p2[0] = 0.5 * ( xt[2] + xt[0] );
  p2[1] = 0.5 * ( yt[2] + yt[0] );

  /* normals */
  n0[0] = - (yt[1] - yt[0]);
  n0[1] =   (xt[1] - xt[0]);
  l = sqrt(n0[0]*n0[0]+n0[1]*n0[1]);
  n0[0] = n0[0] / l;
  n0[1] = n0[1] / l;

  n1[0] = - (yt[2] - yt[1]);
  n1[1] =   (xt[2] - xt[1]);
  l = sqrt(n1[0]*n1[0]+n1[1]*n1[1]);
  n1[0] = n1[0] / l;
  n1[1] = n1[1] / l;

  n2[0] = - (yt[0] - yt[2]);
  n2[1] =   (xt[0] - xt[2]);
  l = sqrt(n2[0]*n2[0]+n2[1]*n2[1]);
  n2[0] = n2[0] / l;
  n2[1] = n2[1] / l;

  Cross_Point(p0, n0, p1, n1, cp01);
  Cross_Point(p1, n1, p2, n2, cp12);
  Cross_Point(p2, n2, p0, n0, cp20);

  /*	if( (sqrt( (cp01[0]-cp12[0])*(cp01[0]-cp12[0]) + (cp01[1]-cp12[1])*(cp01[1]-cp12[1]) )>SMALLDOUBLE) ||
              (sqrt( (cp12[0]-cp20[0])*(cp12[0]-cp20[0]) + (cp12[1]-cp20[1])*(cp12[1]-cp20[1]) )>SMALLDOUBLE) ||
              (sqrt( (cp20[0]-cp01[0])*(cp20[0]-cp01[0]) + (cp20[1]-cp01[1])*(cp20[1]-cp01[1]) )>SMALLDOUBLE))
          {
                  printf("%s\n", "ERROR");
                  midp[0] = maxdouble;
                  midp[1] = maxdouble;
          }
          else
          {
                  midp[0] = cp01[0];
                  midp[1] = cp01[1];
          }*/

  midp[0] = cp01[0];
  midp[1] = cp01[1];
  return(0.0);
}

int GenerateTriangle (ARRAY<Point2d> & lpoints, ARRAY<ILINE> & llines,
                      ARRAY<Element> & elements,
                      ARRAY<INDEX> & dellines,
                      double h)
{
  int i, j, point_i, old, id[4];
  double xt[3],yt[3],np_midp[2],np[2],min_circ,np_circ,circ,midp[2];
  Element element;
  Point2d npoint;

  elements.SetSize (0);
  dellines.SetSize (0);

  min_circ = maxdouble;
  point_i = -1;

  xt[0] = lpoints[llines[1].I1()].X();
  yt[0] = lpoints[llines[1].I1()].Y();
  xt[1] = lpoints[llines[1].I2()].X();
  yt[1] = lpoints[llines[1].I2()].Y();

  CalcNewPoint(xt,yt,h);

  /* first: test the new point */
  if(CheckNewElement(lpoints,llines,xt,yt))
  {
    Calc_Circumcircusmidpoint(xt,yt,np_midp);
    np_circ = sqrt( (np_midp[0]-xt[0])*(np_midp[0]-xt[0]) + (np_midp[1]-yt[0])*(np_midp[1]-yt[0]) );
    np[0] = xt[2];
    np[1] = yt[2];
  }
  else
    np_circ = maxdouble;

  for(i=2; i<=lpoints.Size(); i++)
  {
    xt[2] = lpoints[i].X();
    yt[2] = lpoints[i].Y();

    if( ( ( (xt[1]-xt[0])*(yt[2]-yt[0])-(yt[1]-yt[0])*(xt[2]-xt[0]) )>SMALLDOUBLE ) )
      if(CheckNewElement(lpoints,llines,xt, yt))
      {
        Calc_Circumcircusmidpoint(xt, yt, midp);
        /*				circ0 = sqrt( (midp[0]-xt[0])*(midp[0]-xt[0]) + (midp[1]-yt[0])*(midp[1]-yt[0]) );
                                        circ1 = sqrt( (midp[0]-xt[1])*(midp[0]-xt[1]) + (midp[1]-yt[1])*(midp[1]-yt[1]) );
                                        circ2 = sqrt( (midp[0]-xt[2])*(midp[0]-xt[2]) + (midp[1]-yt[2])*(midp[1]-yt[2]) );
                                        if( (sqrt((circ0-circ1)*(circ0-circ1))<SMALLDOUBLE) &&
                                            (sqrt((circ1-circ2)*(circ1-circ2))<SMALLDOUBLE) &&
                                                (sqrt((circ2-circ0)*(circ2-circ0))<SMALLDOUBLE) )
                                                circ = circ1;
                                        else
                                                printf("%s\n", "ERROR");*/

        circ = sqrt( (midp[0]-xt[0])*(midp[0]-xt[0]) + (midp[1]-yt[0])*(midp[1]-yt[0]) );
        if(min_circ>circ)
        {
          min_circ = circ;
          point_i = i;
        }
      }
  }

  if(min_circ>1.5*np_circ)
  {
    npoint.X() = np[0];
    npoint.Y() = np[1];
    point_i = lpoints.Append (npoint);
  }

  element.PNum(1) = 1;
  element.PNum(2) = 2;
  element.PNum(3) = point_i;
  element.SetNP(3);
  elements.Append (element);

  if(point_i==-1)
    return(-1);

  id[1] = id[2] = id[3] = -1;
  old = llines.Size();
  for(i=1; i<=old; i++)
  {
    for(j=0; j<3; j++)
    {
      if( (element.PNum(j%3+1)==llines[i].I1()) && (element.PNum((j+1)%3+1)==llines[i].I2()))
      {
        dellines.Append(i);
        id[j+1] = 1;
      }
    }
  }
  for(j=0; j<3; j++)
  {
    if(id[j+1]==-1)
      llines.Append( ILINE(element.PNum((j+1)%3+1), element.PNum(j%3+1)) );
  }
  return(1);
}
