// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   nginterface.cc                                               */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <strstream.h>
#include <math.h>
#include <new.h>
#include <ctype.h>
#include <string.h>
#include <time.h>



#include <myadt.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>
#include <linalg/linalg.hh>

// #include <meshing/meshtool.hh>
// #include <meshing/meshing2.hh>
// #include <meshing/meshsurf.hh>

#include <meshing/global.hh>
#include <meshing/meshing3.hh>

static int LGM_DEBUG = 0;

extern "C"
{
#include "../../gginterface.h"
}

int testmode;

extern "C"
void UserWriteF (char * ch, ...);


void MyError (char * ch)
{
  UserWriteF (ch);
}



static ARRAY<Point3d> * points;
static ARRAY<Element> * volelements;
static int nbp,disp;
static double vol0;

class my_meshing3 : public Meshing3
{
public:
  my_meshing3 (char * rulefilename);
  virtual int SavePoint (const Point3d & p);
  virtual void SaveElement (const Element & elem);
  virtual void Get_Local_h_3d(double *in,double *out);

  friend int StartNetgen (double h, int smooth, int display,int prism);
};



my_meshing3 :: my_meshing3 (char * rulefilename)
  : Meshing3 (rulefilename)
{
  ;
}

int my_meshing3 :: SavePoint (const Point3d & p)
{
  return points -> Append (p);
}

#define SMALL 1e-10
#define MAXDOUBLE 1e200
#define PI 3.141592654

#define Lenght(vec)             sqrt(vec[0]*vec[0]      \
                                     +vec[1]*vec[1]  \
                                     +vec[2]*vec[2])

#define Cross(vec,vec1,vec2)    vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]; \
  vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]; \
  vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

#define Minus(sol,vec1,vec2)    sol[0] = vec1[0] - vec2[0];     \
  sol[1] = vec1[1] - vec2[1];     \
  sol[2] = vec1[2] - vec2[2];


#define V3_ADD(A,B,C)                              {(C)[0] = (A)[0] + (B)[0];\
                                                    (C)[1] = (A)[1] + (B)[1];\
                                                    (C)[2] = (A)[2] + (B)[2];}

#define V3_SCALE(c,C)                              {(C)[0] = (c)*(C)[0];\
                                                    (C)[1] = (c)*(C)[1];\
                                                    (C)[2] = (c)*(C)[2];}

#define V3_EUKLIDNORM_OF_DIFF(A,B,b)    (b) = (sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])+((A)[2]-(B)[2])*((A)[2]-(B)[2]))));

#define V3_CLEAR(A)                                {(A)[0] = 0.0; (A)[1]= 0.0; (A)[2] = 0.0;}

static double IN_CIRC(Element elem,float vol)
{
  double p1[3],p2[3],p3[3],p4[3],p5[3],n1[3],n2[3],n3[3],s1,s2,s3,s4,in_circ,rho;
  int i,n;

  p1[0] = (double)points->Get(elem.PNum(1)).X();
  p1[1] = (double)points->Get(elem.PNum(1)).Y();
  p1[2] = (double)points->Get(elem.PNum(1)).Z();

  p2[0] = (double)points->Get(elem.PNum(2)).X();
  p2[1] = (double)points->Get(elem.PNum(2)).Y();
  p2[2] = (double)points->Get(elem.PNum(2)).Z();

  p3[0] = (double)points->Get(elem.PNum(3)).X();
  p3[1] = (double)points->Get(elem.PNum(3)).Y();
  p3[2] = (double)points->Get(elem.PNum(3)).Z();

  p4[0] = (double)points->Get(elem.PNum(4)).X();
  p4[1] = (double)points->Get(elem.PNum(4)).Y();
  p4[2] = (double)points->Get(elem.PNum(4)).Z();

  Minus(n1,p2,p1);
  Minus(n2,p3,p1);
  Cross(n3,n1,n2);
  s1 = Lenght(n3)/2;

  Minus(n1,p4,p1);
  Minus(n2,p2,p1);
  Cross(n3,n1,n2);
  s2 = Lenght(n3)/2;

  Minus(n1,p4,p1);
  Minus(n2,p3,p1);
  Cross(n3,n1,n2);
  s3 = Lenght(n3)/2;

  Minus(n1,p4,p2);
  Minus(n2,p4,p3);
  Cross(n3,n1,n2);
  s4 = Lenght(n3)/2;

  in_circ = - 3 * vol / (s1 + s2 + s3 + s4);

  return(in_circ);
}

static double OUT_CIRC(Element elem,float vol)
{
  double p0[3],p1[3],p2[3],p3[3],n1[3],out_circ;
  double A,B,C,D,E,F,a,b,c,nn;
  int i;

  p0[0] = (double)points->Get(elem.PNum(1)).X();
  p0[1] = (double)points->Get(elem.PNum(1)).Y();
  p0[2] = (double)points->Get(elem.PNum(1)).Z();

  p1[0] = (double)points->Get(elem.PNum(2)).X();
  p1[1] = (double)points->Get(elem.PNum(2)).Y();
  p1[2] = (double)points->Get(elem.PNum(2)).Z();

  p2[0] = (double)points->Get(elem.PNum(3)).X();
  p2[1] = (double)points->Get(elem.PNum(3)).Y();
  p2[2] = (double)points->Get(elem.PNum(3)).Z();

  p3[0] = (double)points->Get(elem.PNum(4)).X();
  p3[1] = (double)points->Get(elem.PNum(4)).Y();
  p3[2] = (double)points->Get(elem.PNum(4)).Z();

  Minus(n1,p0,p1);
  A = Lenght(n1);
  Minus(n1,p1,p2);
  B = Lenght(n1);
  Minus(n1,p2,p0);
  C = Lenght(n1);
  Minus(n1,p0,p3);
  D = Lenght(n1);
  Minus(n1,p1,p3);
  E = Lenght(n1);
  Minus(n1,p2,p3);
  F = Lenght(n1);

  a = C * E;
  b = A * F;
  c = B * D;

  nn = (a+b+c)*(a+b-c)*(a-b+c)*(-a+b+c);
  if(nn<0.0)
    nn = (a+c+b)*(a+c-b)*(a-c+b)*(-a+c+b);
  if(nn<0.0)
    nn = (b+a+c)*(b+a-c)*(b-a+c)*(-b+a+c);
  if(nn<0.0)
    nn = (b+c+a)*(b+c-a)*(b-c+a)*(-b+c+a);
  if(nn<0.0)
    nn = (c+a+b)*(c+a-b)*(c-a+b)*(-c+a+b);
  if(nn<0.0)
    nn = (c+b+a)*(c+b-a)*(c-b+a)*(-c+b+a);

  out_circ = - sqrt(nn) / (24 * vol);

  return(out_circ);
}

static double VOL(Element elem)
{
  double x[4][3],vol;

  x[0][0] = points -> Get(elem.PNum(1)).X();
  x[0][1] = points -> Get(elem.PNum(1)).Y();
  x[0][2] = points -> Get(elem.PNum(1)).Z();
  x[1][0] = points -> Get(elem.PNum(2)).X();
  x[1][1] = points -> Get(elem.PNum(2)).Y();
  x[1][2] = points -> Get(elem.PNum(2)).Z();
  x[2][0] = points -> Get(elem.PNum(3)).X();
  x[2][1] = points -> Get(elem.PNum(3)).Y();
  x[2][2] = points -> Get(elem.PNum(3)).Z();
  x[3][0] = points -> Get(elem.PNum(4)).X();
  x[3][1] = points -> Get(elem.PNum(4)).Y();
  x[3][2] = points -> Get(elem.PNum(4)).Z();

  vol = (         (x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) * (x[3][2] - x[0][2])
                  -       (x[1][0] - x[0][0]) * (x[3][1] - x[0][1]) * (x[2][2] - x[0][2])
                  +       (x[2][0] - x[0][0]) * (x[3][1] - x[0][1]) * (x[1][2] - x[0][2])
                  -       (x[2][0] - x[0][0]) * (x[1][1] - x[0][1]) * (x[3][2] - x[0][2])
                  +       (x[3][0] - x[0][0]) * (x[1][1] - x[0][1]) * (x[2][2] - x[0][2])
                  -       (x[3][0] - x[0][0]) * (x[2][1] - x[0][1]) * (x[1][2] - x[0][2])
                  ) / 6;

  return(vol);
}

void my_meshing3 :: SaveElement (const Element & elem)
{
  double in_circ,out_circ,rho,vol,percent;
  int i,n,j;
  FILE *file;
  char name[10],buff[5];

  if (disp)
  {
    volelements -> Append (elem);
    //		printf("%s %d %d %d %d\n","element",elem.PNum(1)-1,elem.PNum(2)-1,elem.PNum(3)-1,elem.PNum(4)-1);

    percent = 100.0 * adfront->Volume() / vol0;

    if(elem.NP()==4)
    {
      vol = VOL(elem);
      in_circ = IN_CIRC(elem,vol);
      out_circ = OUT_CIRC(elem,vol);
      rho = 3 * in_circ / out_circ;
    }

    if(elem.NP()==4)
      UserWriteF(" ID(Tetrahedra)=%4d in_circ %6.2f out_circ %7.2f rho %4.2f vol %9.6f%\n",
                 volelements -> Size(),in_circ,out_circ,rho,percent);
    if(elem.NP()==5)
      UserWriteF(" ID(Pyramid)=%4d vol %9.6f%\n",
                 volelements -> Size(),percent);
    if(elem.NP()==6)
      UserWriteF(" ID(Prism)=%4d vol %9.6f%\n",
                 volelements -> Size(),percent);
  }
  /*	if(volelements->Size() % 50 == 0)
          {
                  name[0] = 'g';
                  name[1] = 'r';
                  name[2] = 'a';
                  name[3] = 'p';
                  name[4] = 'e';
                  sprintf(buff,"%d",volelements->Size());
                  name[5] = buff[0];
                  name[6] = buff[1];
                  name[7] = buff[2];
                  name[8] = buff[3];
                  name[9] = buff[4];
                  file = fopen(name,"w");
                  fprintf(file, "%s\n", "volmesh");
                  fprintf(file, "%d\n", volelements->Size());
                  for(i=1;i<=volelements->Size();i++)
                  {
                          fprintf(file, "%d\n",volelements->Get(i).NP());
                          for(j=1;j<=volelements->Get(i).NP()-1;j++)
                                  fprintf(file, "%d ",  volelements->Get(i).PNum(j)-1);
                          fprintf(file, "%d\n",  volelements->Get(i).PNum(volelements->Get(i).NP())-1);
                  }
                  fprintf(file, "%d\n", points->Size());
                  for(i=1;i<=points->Size();i++)
                          fprintf(file, "%f %f %f\n", points->Get(i).X(), points->Get(i).Y(), points->Get(i).Z());
                  fclose(file);
          }*/
}

void my_meshing3 :: Get_Local_h_3d(double *in,double *out)
{
  Get_h(in,out);
}


static my_meshing3 * meshing;


int AddSurfaceNode (int nodeid, double x, double y, double z)
{
  points -> Append (Point3d(x, y, z*1));
  meshing -> AddPoint (Point3d(x, y, z*1), nodeid+1);
  return 0;
}


int AddSurfaceTriangle (int node0, int node1, int node2,int prism_flag)
{
  Element elem(3);

  elem.PNum(1) = node0 + 1;
  elem.PNum(2) = node1 + 1;
  elem.PNum(3) = node2 + 1;
  elem.SetSurfaceIndex (1);
  meshing -> AddBoundaryElement (elem, 0, prism_flag);

  return 0;
}

static int openkey=0;

int InitNetgen (char * rulefilename)
{
  if(openkey==0) {testout = new ofstream("test.out"); openkey=1;}
  meshing = new my_meshing3(rulefilename);

  points = new ARRAY<Point3d>;
  volelements = new ARRAY<Element>;
  nbp = 0;

  return 0;
}

int StartNetgen (double h, int smooth, int display,int prism)
{
  int i;

  nbp = points->Size();
  disp = display;
  vol0 = meshing->adfront->Volume();

  meshing -> Mesh (h,prism);

  for (i=0; i<smooth; i++)
    meshing -> ImproveMesh (*points, *volelements, nbp, h);

  AllMemInnerPoints(points -> Size()-nbp);
  for (i = nbp + 1; i <= points -> Size(); i++)
  {
    AddInnerNode (  points -> Get(i).X(),
                    points -> Get(i).Y(),
                    points -> Get(i).Z()/1);
  }

  AllMemElements(volelements -> Size());
  for (i = 1; i <= volelements -> Size(); i++)
  {
    AddElement (volelements -> Get(i).NP(),
                volelements -> Get(i).PNum(1) - 1,
                volelements -> Get(i).PNum(2) - 1,
                volelements -> Get(i).PNum(3) - 1,
                volelements -> Get(i).PNum(4) - 1,
                volelements -> Get(i).PNum(5) - 1,
                volelements -> Get(i).PNum(6) - 1);
  }
  volelements->SetSize(0);
  points->SetSize(0);
  return (0);
}
