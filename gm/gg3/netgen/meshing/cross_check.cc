// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdio.h>
#include <math.h>
#include <math.h>

#include "cross_check.hh"


double ex[3], ey[3], ez[3], pp[3];

static int Check_One_point_common(      double p1[3][3], int *id1,
                                        double p2[3][3], int *id2,
                                        int i0, int i1, int i2, int j0, int j1, int j2,
                                        double *p10, double *p11, double *p12,
                                        double *p20, double *p21, double *p22);
static int Calculate_Point(double *n1, double *p1_0, double *n2, double *p2_0, double *sp, double *ep, double *p);

static double ABS(double a)
{
  if(a>=0.0)
    return(a);
  else
    return(-a);
}

int Element_Inside_Check(       double *p0, double *p1, double *p2,
                                double *p3, double *p4)
{
  int i;
  double a[3][3], rhs[3], sol[4];

  for(i=0; i<3; i++)
  {
    a[i][0] = p0[i] - p3[i];
    a[i][1] = p1[i] - p3[i];
    a[i][2] = p2[i] - p3[i];
    rhs[i] = p4[i] - p3[i];
  }

  InvMatMult3d(sol,rhs,a);

  sol[3] = 1 - sol[0] - sol[1] - sol[2];

  if( (0<sol[0])&&(sol[0]<1.0)
      &&      (0<sol[1])&&(sol[1]<1.0)
      &&      (0<sol[2])&&(sol[2]<1.0)
      &&      (0<sol[3])&&(sol[3]<1.0) )
    return(1);
  else
    return(0);
}

double Delaunay(double *p1, double *p2, double *p3, double *p4, double *p5)
{
  int i;
  double a[3][3], rhs[3], det, n[3], l1;

  rhs[0] = rhs[1] = rhs[2] = 0.0;
  for(i=0; i<3; i++)
  {
    a[0][i] = p1[i] - p2[i];
    a[1][i] = p1[i] - p3[i];
    a[2][i] = p1[i] - p4[i];
    rhs[0] = rhs[0] + (p1[i]*p1[i] - p2[i]*p2[i]) * 0.5;
    rhs[1] = rhs[1] + (p1[i]*p1[i] - p3[i]*p3[i]) * 0.5;
    rhs[2] = rhs[2] + (p1[i]*p1[i] - p4[i]*p4[i]) * 0.5;
  }
  det = Det3d(a);
  if(ABS(det)>SMALL)
  {
    InvMatMult3d(p5,rhs,a);
  }
  else
  {
    for(i=0; i<3; i++)
    {
      a[0][i] = p2[i] - p1[i];
      a[1][i] = p2[i] - p3[i];
      a[2][i] = p2[i] - p4[i];
      rhs[0] = rhs[0] + (p2[i]*p2[i] - p1[i]*p1[i]) * 0.5;
      rhs[1] = rhs[1] + (p2[i]*p2[i] - p3[i]*p3[i]) * 0.5;
      rhs[2] = rhs[2] + (p2[i]*p2[i] - p4[i]*p4[i]) * 0.5;
    }
    det = Det3d(a);
    if(ABS(det)>SMALL)
    {
      InvMatMult3d(p5,rhs,a);
    }
    else
    {
      for(i=0; i<3; i++)
      {
        a[0][i] = p3[i] - p1[i];
        a[1][i] = p3[i] - p2[i];
        a[2][i] = p3[i] - p4[i];
        rhs[0] = rhs[0] + (p3[i]*p3[i] - p1[i]*p1[i]) * 0.5;
        rhs[1] = rhs[1] + (p3[i]*p3[i] - p2[i]*p2[i]) * 0.5;
        rhs[2] = rhs[2] + (p3[i]*p3[i] - p4[i]*p4[i]) * 0.5;
      }
      det = Det3d(a);
      if(ABS(det)>SMALL)
      {
        InvMatMult3d(p5,rhs,a);
      }
      else
      {
        for(i=0; i<3; i++)
        {
          a[0][i] = p4[i] - p1[i];
          a[1][i] = p4[i] - p2[i];
          a[2][i] = p4[i] - p4[i];
          rhs[0] = rhs[0] + (p4[i]*p4[i] - p1[i]*p1[i]) * 0.5;
          rhs[1] = rhs[1] + (p4[i]*p4[i] - p2[i]*p2[i]) * 0.5;
          rhs[2] = rhs[2] + (p4[i]*p4[i] - p3[i]*p3[i]) * 0.5;
        }
        det = Det3d(a);
        if(ABS(det)>SMALL)
        {
          InvMatMult3d(p5,rhs,a);
        }
        else
        {
          //					printf("%s\n", "S L I V E R");
          return(-1.0);
        }
      }
    }
  }

  /* check distance */
  Minus(n,p1,p5);
  l1 = Lenght(n);

  return(l1);
}

double Calc_D_M(double *p1,double *p2,double *p3,double *D,double *M)
{
  int i;
  double dist,l,m[3];

  dist = 1000000000.0;
  for(i=0; i<3; i++)
    M[i] = ( p1[i] + p2[i] + p3[i] ) / 3;

  Minus(m,p1,M);
  l = Lenght(m);
  if(l<dist)
  {
    dist = l;
    VecCopy3d(D, p1);
  }
  Minus(m,p2,M);
  l = Lenght(m);
  if(l<dist)
  {
    dist = l;
    VecCopy3d(D, p2);
  }
  Minus(m,p3,M);
  l = Lenght(m);
  if(l<dist)
  {
    dist = l;
    VecCopy3d(D, p3);
  }

  return(0);
}

double Sphere_radius(double *Q,double *P,double *D,double *M)
{
  int i;
  double dist2,dist3,m[3],lb,lc,lambda,PP[3];

  lb =    ( D[0] - M[0] ) * ( D[0] - M[0] )
       +       ( D[1] - M[1] ) * ( D[1] - M[1] )
       +       ( D[2] - M[2] ) * ( D[2] - M[2] )
       -       ( Q[0] - M[0] ) * ( Q[0] - M[0] )
       -       ( Q[1] - M[1] ) * ( Q[1] - M[1] )
       -       ( Q[2] - M[2] ) * ( Q[2] - M[2] );
  lc =    ( P[0] - M[0] ) * ( D[0] - Q[0] )
       +       ( P[1] - M[1] ) * ( D[1] - Q[1] )
       +       ( P[2] - M[2] ) * ( D[2] - Q[2] );
  lambda = lb / (2 *lc);

  for(i=0; i<3; i++)
    PP[i] = M[i] + lambda * ( P[i] - M[i]);

  /* check */
  Minus(m,M,PP);
  dist2 = Lenght(m);
  Minus(m,D,PP);
  dist2 = Lenght(m);
  Minus(m,Q,PP);
  dist3 = Lenght(m);
  if((ABS(dist2-dist3)<SMALL))
    return(dist2);
  else
    return(100.0);
}


static double MAX2(double a,  double b)
{
  if(a>b)
    return(a);
  else
    return(b);
}

static double MAX3(double a,  double b, double c)
{
  return(MAX2(MAX2(a, b), c));
}

static double MIN2(double a,  double b)
{
  if(a<b)
    return(a);
  else
    return(b);
}

static double MIN3(double a,  double b, double c)
{
  return(MIN2(MIN2(a, b), c));
}

static int Bounding_Box_3d(     double *p1_0,
                                double *p1_1,
                                double *p1_2,
                                double *p2_0,
                                double *p2_1,
                                double *p2_2)
{
  double max1[3], min1[3], max2[3], min2[3];

  max1[0] = MAX3(p1_0[0], p1_1[0], p1_2[0]);
  max1[1] = MAX3(p1_0[1], p1_1[1], p1_2[1]);
  max1[2] = MAX3(p1_0[2], p1_1[2], p1_2[2]);

  min1[0] = MIN3(p1_0[0], p1_1[0], p1_2[0]);
  min1[1] = MIN3(p1_0[1], p1_1[1], p1_2[1]);
  min1[2] = MIN3(p1_0[2], p1_1[2], p1_2[2]);

  max2[0] = MAX3(p2_0[0], p2_1[0], p2_2[0]);
  max2[1] = MAX3(p2_0[1], p2_1[1], p2_2[1]);
  max2[2] = MAX3(p2_0[2], p2_1[2], p2_2[2]);

  min2[0] = MIN3(p2_0[0], p2_1[0], p2_2[0]);
  min2[1] = MIN3(p2_0[1], p2_1[1], p2_2[1]);
  min2[2] = MIN3(p2_0[2], p2_1[2], p2_2[2]);

  if( ( (max1[0]<min2[0]) || (max1[1]<min2[1]) || (max1[2]<min2[2]) )
      ||  ( (min1[0]>max2[0]) || (min1[1]>max2[1]) || (min1[2]>max2[2]) ) )
  {
    //		printf("%s\n", "Bounding_Box_3d: keine Ueberschneidung");
    return(0);                          /* keine Ueberschneidung der Bounding Boxes */
  }
  else
  {
    //		printf("%s\n", "Bounding_Box_3d: Ueberschneidung ist moeglich");
    return(1);
  }
}

int GetNormalVector(double *p0, double *p1, double *p2, double *n)
{
  double n1[3], n2[3], l;

  Minus(n1, p1, p0);
  l = Lenght(n1);
  Scale(n1,1/l);
  Minus(n2, p2, p0);
  l = Lenght(n2);
  Scale(n2,1/l);
  Cross(n, n1, n2);
  l = Lenght(n);
  Scale(n,1/l);

  return(1);
}

static int Faces_in_the_same_Plane(     double *p1_0, double *p1_1, double *p1_2,
                                        double *p2_0, double *p2_1, double *p2_2, int *flag)
{
  double n1[3], n2[3], sp, m[3], a;

  GetNormalVector(p1_0, p1_1, p1_2, n1);
  GetNormalVector(p2_0, p2_1, p2_2, n2);

  *flag = 0;
  sp = Mult(n1, n2);
  if((sp>1-SMALL)||(sp<-1+SMALL))
  {
    //		printf("%s\n", "Faces_in_the_same_Plane: Faces nearly planar");
    /* Check, if the triangles define the same plane */
    Minus(m, p1_0, p2_0);
    a = Mult(n1, m);
    if(ABS(a)<SMALL)
      *flag = 1;
    return(1);
  }
  else
  {
    //		printf("%s\n", "Faces_in_the_same_Plane: Faces not planar");
    return(0);
  }
}

static int Two_points_common(int *id1, int*id2)
{
  int i, j, f[3], g;

  f[0] = f[1] = f[2] = 0;
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      if(id1[i]==id2[j])
        f[i]++;

  g = 0;
  for(i=0; i<3; i++)
    if(f[i]==1)
      g++;
  if(g==2)
    return(1);
  else
    return(0);
}

static int One_point_common(double *p1_0, double *p1_1, double *p1_2, int *id1,
                            double *p2_0, double *p2_1, double *p2_2, int *id2,
                            double p1[3][3], double p2[3][3], int *flag1)
{
  int i, flag, comp;
  double n1[3], n2[3], m2[3], m3[3],b, c;
  double p10[3], p11[3], p12[3], p20[3], p21[3], p22[3];
  double pt1[3], pt2[3], m[3], l, max;

  for(i=0; i<3; i++)
  {
    p1[i][0] = p1_0[i]; p1[i][1] = p1_1[i]; p1[i][2] = p1_2[i];
    p2[i][0] = p2_0[i]; p2[i][1] = p2_1[i]; p2[i][2] = p2_2[i];
  }
  flag = 0;

  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 0, 2, 1, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 1, 0, 2, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 2, 1, 0, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 2, 0, 1, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 0, 1, 2, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 1, 2, 0, p10, p11, p12, p20, p21, p22)) flag = 1;else

  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 0, 2, 1, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 1, 0, 2, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 2, 1, 0, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 2, 0, 1, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 0, 1, 2, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 1, 2, 0, p10, p11, p12, p20, p21, p22)) flag = 1;else

  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 0, 2, 1, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 1, 0, 2, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 2, 1, 0, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 2, 0, 1, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 0, 1, 2, p10, p11, p12, p20, p21, p22)) flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 1, 2, 0, p10, p11, p12, p20, p21, p22)) flag = 1;

  if(flag)
  {
    GetNormalVector(p1_0, p1_1, p1_2, n1);
    for(i=0; i<3; i++)
    {
      m2[i] = p10[i] - p21[i];
      m3[i] = p10[i] - p22[i];
    }
    b = Mult(n1, m2);
    c = Mult(n1, m3);
    if( (b<-SMALL)&&(c<-SMALL) )
    {
      *flag1 = 0;
      return(1);
    }
    if( (b>SMALL)&&(c>SMALL) )
    {
      *flag1 = 0;
      return(1);
    }

    GetNormalVector(p2_0, p2_1, p2_2, n1);
    for(i=0; i<3; i++)
    {
      m2[i] = p20[i] - p11[i];
      m3[i] = p20[i] - p12[i];
    }
    b = Mult(n1, m2);
    c = Mult(n1, m3);
    if( (b<-SMALL)&&(c<-SMALL) )
    {
      *flag1 = 0;
      return(1);
    }
    if( (b>SMALL)&&(c>SMALL) )
    {
      *flag1 = 0;
      return(1);
    }

    GetNormalVector(p1_0, p1_1, p1_2, n1);
    GetNormalVector(p2_0, p2_1, p2_2, n2);
    Calculate_Point(n1, p1_0, n2, p2_0, p11, p12, pt1);
    Calculate_Point(n1, p1_0, n2, p2_0, p21, p22, pt2);

    Minus(m, pt2, pt1);
    l = Lenght(m);
    if(l>SMALL)
      Scale(m,1/l);

    max = 0.0;
    for(i=0; i<3; i++)
    {
      if(max<ABS(m[i]))
      {
        comp = i;
        max = ABS(m[i]);
      }
    }

    if( ((pt1[comp]<p10[comp])&&(p10[comp]<pt2[comp])) || ((pt1[comp]>p10[comp])&&(p10[comp]>pt2[comp])) )
      *flag1 = 0;                                       /* no intersection */
    else
      *flag1 = 1;

    return(1);
  }
  else
  {
    *flag1 = 0;
    return(0);
  }
}

static int Faces_on_the_same_side(      double *p1_0, double *p1_1, double *p1_2,
                                        double *p2_0, double *p2_1, double *p2_2)
{
  double n1[3], m1[3], m2[3], m3[3], a, b, c;

  GetNormalVector(p1_0, p1_1, p1_2, n1);
  Minus(m1, p1_0, p2_0);
  Minus(m2, p1_0, p2_1);
  Minus(m3, p1_0, p2_2);
  a = Mult(n1, m1);
  b = Mult(n1, m2);
  c = Mult(n1, m3);
  if( ( ( a<-SMALL )&&( b<-SMALL )&&( c<-SMALL ) )
      ||  ( ( a> SMALL )&&( b> SMALL )&&( b> SMALL ) ) )
    return(1);

  GetNormalVector(p2_0, p2_1, p2_2, n1);
  Minus(m1, p2_0, p1_0);
  Minus(m2, p2_0, p1_1);
  Minus(m3, p2_0, p1_2);
  a = Mult(n1, m1);
  b = Mult(n1, m2);
  c = Mult(n1, m3);
  if( ( ( a<-SMALL )&&( b<-SMALL )&&( c<-SMALL ) )
      ||  ( ( a> SMALL )&&( b> SMALL )&&( b> SMALL ) ) )
    return(1);

  return(0);
}

static int Faces_nearly_on_the_same_side(       double *p1_0, double *p1_1, double *p1_2, int *id1,
                                                double *p2_0, double *p2_1, double *p2_2, int *id2)
{
  double n1[3], m1[3], m2[3], m3[3], a, b, c;

  GetNormalVector(p1_0, p1_1, p1_2, n1);
  Minus(m1, p1_0, p2_0);
  Minus(m2, p1_0, p2_1);
  Minus(m3, p1_0, p2_2);
  a = Mult(n1, m1);
  b = Mult(n1, m2);
  c = Mult(n1, m3);

  if( ( (ABS(a)<SMALL)&&(b<-SMALL)&&(c<-SMALL) )
      || ( (a<-SMALL)&&(ABS(b)<SMALL)&&(c<-SMALL) )
      || ( (a<-SMALL)&&(b<-SMALL)&&(ABS(c)<SMALL) ) )
    return(1);
  if( ( (ABS(a)<SMALL)&&(b>SMALL)&&(c>SMALL) )
      || ( (a>SMALL)&&(ABS(b)<SMALL)&&(c>SMALL) )
      || ( (a>SMALL)&&(b>SMALL)&&(ABS(c)<SMALL) ) )
    return(1);
  if( ( (ABS(a)<SMALL)&&(ABS(b)<SMALL)&&(c<-SMALL) )
      || ( (ABS(a)<SMALL)&&(b<-SMALL)&&(ABS(c)<SMALL) )
      || ( (a<-SMALL)&&(ABS(b)<SMALL)&&(ABS(c)<SMALL) ) )
    return(1);
  if( ( (ABS(a)<SMALL)&&(ABS(b)<SMALL)&&(c>SMALL) )
      || ( (ABS(a)<SMALL)&&(b>SMALL)&&(ABS(c)<SMALL) )
      || ( (a>SMALL)&&(ABS(b)<SMALL)&&(ABS(c)<SMALL) ) )
    return(1);

  GetNormalVector(p2_0, p2_1, p2_2, n1);
  Minus(m1, p2_0, p1_0);
  Minus(m2, p2_0, p1_1);
  Minus(m3, p2_0, p1_2);
  a = Mult(n1, m1);
  b = Mult(n1, m2);
  c = Mult(n1, m3);

  if( ( (ABS(a)<SMALL)&&(b<-SMALL)&&(c<-SMALL) )
      || ( (a<-SMALL)&&(ABS(b)<SMALL)&&(c<-SMALL) )
      || ( (a<-SMALL)&&(b<-SMALL)&&(ABS(c)<SMALL) ) )
    return(1);
  if( ( (ABS(a)<SMALL)&&(b>SMALL)&&(c>SMALL) )
      || ( (a>SMALL)&&(ABS(b)<SMALL)&&(c>SMALL) )
      || ( (a>SMALL)&&(b>SMALL)&&(ABS(c)<SMALL) ) )
    return(1);
  if( ( (ABS(a)<SMALL)&&(ABS(b)<SMALL)&&(c<-SMALL) )
      || ( (ABS(a)<SMALL)&&(b<-SMALL)&&(ABS(c)<SMALL) )
      || ( (a<-SMALL)&&(ABS(b)<SMALL)&&(ABS(c)<SMALL) ) )
    return(1);
  if( ( (ABS(a)<SMALL)&&(ABS(b)<SMALL)&&(c>SMALL) )
      || ( (ABS(a)<SMALL)&&(b>SMALL)&&(ABS(c)<SMALL) )
      || ( (a>SMALL)&&(ABS(b)<SMALL)&&(ABS(c)<SMALL) ) )
    return(1);

  return(0);
}

static int Check_Two_points_common(     double p1[3][3], int *id1,
                                        double p2[3][3], int *id2,
                                        int i0, int i1, int i2, int j0, int j1, int j2,
                                        double *p10, double *p11, double *p12,
                                        double *p20, double *p21, double *p22)
{
  int i;

  if( (id1[i0]==id2[j0])&&(id1[i1]==id2[j1])&&(id1[i2]!=id2[j2]) )
  {
    for(i=0; i<3; i++)
    {
      p10[i]=p1[i][i0];
      p11[i]=p1[i][i1];
      p12[i]=p1[i][i2];

      p20[i]=p2[i][j0];
      p21[i]=p2[i][j1];
      p22[i]=p2[i][j2];
    }
    return(1);
  }
  else
    return(0);
}

static int Define_Tangential_Plane(double p[3][3])
{
  int i;
  double l, e1[3], e2[3];

  for(i=0; i<3; i++)
  {
    pp[i] = p[i][0];
    e1[i] = p[i][1] - p[i][0];
    e2[i] = p[i][2] - p[i][0];
    ex[i] = p[i][1] - p[i][0];
  }

  l = Lenght(e1);
  if(l>SMALL)
    Scale(e1,1/l);
  l = Lenght(e2);
  if(l>SMALL)
    Scale(e2,1/l);
  Cross(ez, e1, e2);
  l = Lenght(ez);
  if(l>SMALL)
    Scale(ez,1/l);

  l = Lenght(ex);
  if(l>SMALL)
    Scale(ex,1/l);

  l = Mult(ex, ez);
  for(i=0; i<3; i++)
    ex[i] = ex[i] - l * ez[i];
  l = Lenght(ex);
  if(l>SMALL)
    Scale(ex,1/l);

  Cross(ey, ez, ex);
  l = Lenght(ey);
  Scale(ey,1/l);

  return(1);
}

static int Project_Point_to_plane(double *p, double *p_plain)
{
  int i;
  double help[3], help1[3], help2[3], l;

  Minus(help, p, pp);
  l = Mult(help, ez);
  for(i=0; i<3; i++)
    help1[i] = ez[i] * l;

  Minus(help2, help, help1);

  p_plain[0] = Mult(help2, ex);
  p_plain[1] = Mult(help2, ey);

  return(0);
}

static int Check_Cross(double *ps, double *pe, double *p1, double *p2)
{
  double m[2], n[2], l, c, v1, v2;

  /* HNF der Geraden durch ps und pe */
  m[0] = pe[0] - ps[0];
  m[1] = pe[1] - ps[1];

  l = sqrt(m[0]*m[0]+m[1]*m[1]);
  m[0] = m[0] / l;
  m[1] = m[1] / l;

  n[0] = m[1];
  n[1] = -m[0];
  c = ps[0]*n[0] + ps[1]*n[1];

  v1 = n[0]*p1[0] + n[1]*p1[1] - c;
  v2 = n[0]*p2[0] + n[1]*p2[1] - c;

  if( ((v1>0.0) && (v2<0.0)) || ((v1<0.0) && (v2>0.0)) )
    return(0);                                                  /* kein Schnittpunkt */
  else
    return(1);
}

static int In_Plane_Two_points_common(  double *p1_0, double *p1_1, double *p1_2, int *id1,
                                        double *p2_0, double *p2_1, double *p2_2, int *id2, int *flag)
{
  int i;
  double p10[3], p11[3], p12[3], p20[3], p21[3], p22[3], p1[3][3], p2[3][3];
  double p_10[2], p_11[2], p_12[2], p_20[2], p_21[2], p_22[2];

  for(i=0; i<3; i++)
  {
    p1[i][0] = p1_0[i]; p1[i][1] = p1_1[i]; p1[i][2] = p1_2[i];
    p2[i][0] = p2_0[i]; p2[i][1] = p2_1[i]; p2[i][2] = p2_2[i];
  }
  *flag = 0;

  if(Check_Two_points_common(p1, id1, p2, id2, 0, 1, 2, 0, 2, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 0, 1, 2, 1, 0, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 0, 1, 2, 2, 1, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 0, 1, 2, 2, 0, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 0, 1, 2, 0, 1, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 0, 1, 2, 1, 2, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else

  if(Check_Two_points_common(p1, id1, p2, id2, 2, 0, 1, 0, 2, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 2, 0, 1, 1, 0, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 2, 0, 1, 2, 1, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 2, 0, 1, 2, 0, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 2, 0, 1, 0, 1, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 2, 0, 1, 1, 2, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else

  if(Check_Two_points_common(p1, id1, p2, id2, 1, 2, 0, 0, 2, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 1, 2, 0, 1, 0, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 1, 2, 0, 2, 1, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 1, 2, 0, 2, 0, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 1, 2, 0, 0, 1, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_Two_points_common(p1, id1, p2, id2, 1, 2, 0, 1, 2, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;

  if(*flag)
  {
    Define_Tangential_Plane(p1);
    Project_Point_to_plane(p10, p_10);
    Project_Point_to_plane(p11, p_11);
    Project_Point_to_plane(p12, p_12);
    Project_Point_to_plane(p20, p_20);
    Project_Point_to_plane(p21, p_21);
    Project_Point_to_plane(p22, p_22);
    //		printf("%s\n", "In_Plane_Two_points_common: 2 Punkte gemeinsam");

    if(Check_Cross(p_10, p_11, p_12, p_22))
    {
      //			printf("%s\n", "Schnittpunkt");
      return(1);
    }
    else
    {
      //			printf("%s\n", "kein Schnittpunkt");
      *flag = 0;
      return(1);
    }
  }

  return(0);
}

static int Check_One_point_common(      double p1[3][3], int *id1,
                                        double p2[3][3], int *id2,
                                        int i0, int i1, int i2, int j0, int j1, int j2,
                                        double *p10, double *p11, double *p12,
                                        double *p20, double *p21, double *p22)
{
  int i;

  if( (id1[i0]==id2[j0])&&(id1[i1]!=id2[j1])&&(id1[i2]!=id2[j2]) )
  {
    for(i=0; i<3; i++)
    {
      p10[i]=p1[i][i0];
      p11[i]=p1[i][i1];
      p12[i]=p1[i][i2];

      p20[i]=p2[i][j0];
      p21[i]=p2[i][j1];
      p22[i]=p2[i][j2];
    }
    return(1);
  }
  else
    return(0);
}

static double Get_Angle(double *n)
{
  double l, a, b, n1[2], n2[2], p2;

  p2 = PI/2;
  l = sqrt(n[0]*n[0]+n[1]*n[1]);

  n1[0] = 1.0;
  n1[1] = 0.0;
  n2[0] = 0.0;
  n2[1] = 1.0;

  a = acos( (n[0]*n1[0]+n[1]*n1[1]) / l);
  b = acos( (n[0]*n2[0]+n[1]*n2[1]) / l);

  if( (a>=0.0)&&(a<p2)&&(b>=0.0)&&(b<p2) )
    return(a);
  if( (a>=p2)&&(a<2*p2)&&(b>=0.0)&&(b<p2) )
    return(a);
  if( (a>=p2)&&(a<2*p2)&&(b>=p2)&&(b<2*p2) )
    return(4*p2-a);
  if( (a>=0)&&(a<p2)&&(b>=p2)&&(b<2*p2) )
    return(4*p2-a);
  return(0.0);
}

static int Check_Angle(double *p, double *p1, double *p2, double *p3, double *p4)
{
  int i;
  double n1[2], n2[2], n3[2], n4[2], a1, a2, b1, b2, a[2];

  for(i=0; i<2; i++)
  {
    n1[i] = p1[i] - p[i];
    n2[i] = p2[i] - p[i];
    n3[i] = p3[i] - p[i];
    n4[i] = p4[i] - p[i];
  }

  a[0] = Get_Angle(n1);
  a[1] = Get_Angle(n2);

  if(a[0]<a[1])
  {
    a1 = a[0];
    b1 = a[1];
  }
  else
  {
    a1 = a[1];
    b1 = a[0];
  }

  a[0] = Get_Angle(n3);
  a[1] = Get_Angle(n4);

  if(a[0]<a[1])
  {
    a2 = a[0];
    b2 = a[1];
  }
  else
  {
    a2 = a[1];
    b2 = a[0];
  }
  if( (b1-a1<PI)&&(a1<a2)&&(a2<b1)
      ||  (b1-a1<PI)&&(a1<b2)&&(b2<b1)
      ||  (b2-a2<PI)&&(a2<a1)&&(a1<b2)
      ||  (b2-a2<PI)&&(a2<b1)&&(b1<b2)
      ||  (b1-a1>=PI)&&!((a1<a2)&&(a2<b1))
      ||  (b1-a1>=PI)&&!((a1<b2)&&(b2<b1))
      ||  (b2-a2>=PI)&&!((a2<a1)&&(a1<b2))
      ||  (b2-a2>=PI)&&!((a2<b1)&&(b1<b2)) )
    return(1);
  else
    return(0);
}

static int In_Plane_One_point_common(   double *p1_0, double *p1_1, double *p1_2, int *id1,
                                        double *p2_0, double *p2_1, double *p2_2, int *id2, int *flag)
{
  int i;
  double p10[3], p11[3], p12[3], p20[3], p21[3], p22[3], p1[3][3], p2[3][3];
  double p_10[2], p_11[2], p_12[2], p_20[2], p_21[2], p_22[2];

  for(i=0; i<3; i++)
  {
    p1[i][0] = p1_0[i]; p1[i][1] = p1_1[i]; p1[i][2] = p1_2[i];
    p2[i][0] = p2_0[i]; p2[i][1] = p2_1[i]; p2[i][2] = p2_2[i];
  }
  *flag = 0;

  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 0, 2, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 1, 0, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 2, 1, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 2, 0, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 0, 1, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 0, 1, 2, 1, 2, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else

  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 0, 2, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 1, 0, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 2, 1, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 2, 0, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 0, 1, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 2, 0, 1, 1, 2, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else

  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 0, 2, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 1, 0, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 2, 1, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 2, 0, 1, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 0, 1, 2, p10, p11, p12, p20, p21, p22)) *flag = 1;else
  if(Check_One_point_common(p1, id1, p2, id2, 1, 2, 0, 1, 2, 0, p10, p11, p12, p20, p21, p22)) *flag = 1;

  if(*flag)
  {
    Define_Tangential_Plane(p1);
    Project_Point_to_plane(p10, p_10);
    Project_Point_to_plane(p11, p_11);
    Project_Point_to_plane(p12, p_12);
    Project_Point_to_plane(p20, p_20);
    Project_Point_to_plane(p21, p_21);
    Project_Point_to_plane(p22, p_22);
    //		printf("%s\n", "In_Plane_One_point_common: 1 Punkt gemeinsam");

    if(Check_Angle(p_10, p_11, p_12, p_21, p_22))
    {
      //			printf("%s\n", "Schnittpunkt");
      return(1);
    }
    else
    {
      //			printf("%s\n", "kein Schnittpunkt");
      *flag = 0;
      return(1);
    }
  }

  return(0);
}

static int Line_Cross_2d(double *p1_0, double *p1_1, double *p2_0, double *p2_1)
{
  double n1[2], n2[2], l, m, a[2][2], rhs[2], sol[2], p1[2], p2[2];

  n1[0] = p1_1[0] - p1_0[0];
  n1[1] = p1_1[1] - p1_0[1];

  n2[0] = p2_1[0] - p2_0[0];
  n2[1] = p2_1[1] - p2_0[1];

  a[0][0] = n1[0];
  a[0][1] = -n2[0];
  a[1][0] = n1[1];
  a[1][1] = -n2[1];

  rhs[0] = p2_0[0] - p1_0[0];
  rhs[1] = p2_0[1] - p1_0[1];

  if( (Det2d(a)>SMALL) || (Det2d(a)<-SMALL) )
  {
    InvMatMult2d(sol, rhs, a);

    l = sol[0];
    m = sol[1];

    p1[0] = p1_0[0] + l * n1[0];
    p1[1] = p1_0[1] + l * n1[1];
    p2[0] = p2_0[0] + m * n2[0];
    p2[1] = p2_0[1] + m * n2[1];

    if(sqrt( (p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]) )>SMALL)
      printf("%s\n", "error");

    if( (-SMALL<l)&&(l<1.0+SMALL)&&(-SMALL<m)&&(m<1.0+SMALL) )
      return(1);                                                        /* intersection */
    else
      return(0);
  }
  else
    return(0);
}

static int Cut_Test_2d( double *p1_0, double *p1_1, double *p1_2,
                        double *p2_0, double *p2_1, double *p2_2)
{
  int i;
  double p1[3][3];
  double p_10[2], p_11[2], p_12[2], p_20[2], p_21[2], p_22[2];

  for(i=0; i<3; i++)
  {
    p1[i][0] = p1_0[i]; p1[i][1] = p1_1[i]; p1[i][2] = p1_2[i];
  }

  Define_Tangential_Plane(p1);
  Project_Point_to_plane(p1_0, p_10);
  Project_Point_to_plane(p1_1, p_11);
  Project_Point_to_plane(p1_2, p_12);
  Project_Point_to_plane(p2_0, p_20);
  Project_Point_to_plane(p2_1, p_21);
  Project_Point_to_plane(p2_2, p_22);

  if(Line_Cross_2d(p_10, p_11, p_20, p_21)) return(1);else
  if(Line_Cross_2d(p_10, p_11, p_21, p_22)) return(1);else
  if(Line_Cross_2d(p_10, p_11, p_22, p_20)) return(1);else
  if(Line_Cross_2d(p_11, p_12, p_20, p_21)) return(1);else
  if(Line_Cross_2d(p_11, p_12, p_21, p_22)) return(1);else
  if(Line_Cross_2d(p_11, p_12, p_22, p_20)) return(1);else
  if(Line_Cross_2d(p_12, p_10, p_20, p_21)) return(1);else
  if(Line_Cross_2d(p_12, p_10, p_21, p_22)) return(1);else
  if(Line_Cross_2d(p_12, p_10, p_22, p_20)) return(1);else return(0);

}

static int Points_on_different_sides(double *n, double *p, double *p1, double *p2, double *sp, double *ep)
{
  double m1[3], m2[3], a, b;

  Minus(m1, p, p1);
  Minus(m2, p, p2);
  a = Mult(n, m1);
  b = Mult(n, m2);

  if( ((a<=0.0)&&(b>0)) || ((b<=0.0)&&(a>0)) )
  {
    /* endpoint on different sides of the face */
    sp[0] = p1[0];
    sp[1] = p1[1];
    sp[2] = p1[2];
    ep[0] = p2[0];
    ep[1] = p2[1];
    ep[2] = p2[2];
    return(1);
  }
  else
    return(0);
}

static int Calculate_Point(double *n1, double *p1_0, double *n2, double *p2_0, double *sp, double *ep, double *p)
{
  double n3[3], p1[3], p2[3], p3[3], det1, det2, det3, det, rhs1[3], rhs2[3], rhs3[3], help[3];
  double a1[3][3], a2[3][3], a3[3][3], l;

  p1[0] = p1[1] = p1[2] = 0.0;
  p2[0] = p2[1] = p2[2] = 0.0;
  p3[0] = p3[1] = p3[2] = 0.0;
  det1 = det2 = det3 = 0.0;

  n3[0] = ep[1] - sp[1];
  n3[1] = sp[0] - ep[0];
  n3[2] = 0.0;
  l = Lenght(n3);
  if(l>SMALL)
  {
    Scale(n3,1/l);
    help[0] = sp[0];
    help[1] = sp[1];
    help[2] = 0.0;
    rhs1[0] = Mult(n1, p1_0);
    rhs1[1] = Mult(n2, p2_0);
    rhs1[2] = Mult(n3, help);
    a1[0][0] = n1[0]; a1[0][1] = n1[1]; a1[0][2] = n1[2];
    a1[1][0] = n2[0]; a1[1][1] = n2[1]; a1[1][2] = n2[2];
    a1[2][0] = n3[0]; a1[2][1] = n3[1]; a1[2][2] = n3[2];
    det1 = ABS(Det3d(a1));
    if(det1>SMALL)
    {
      InvMatMult3d(p1,rhs1,a1);
    }
  }

  n3[0] = 0.0;
  n3[1] = ep[2] - sp[2];
  n3[2] = sp[1] - ep[1];
  l = Lenght(n3);
  if(l>SMALL)
  {
    Scale(n3,1/l);
    help[0] = 0.0;
    help[1] = sp[1];
    help[2] = sp[2];
    rhs2[0] = Mult(n1, p1_0);
    rhs2[1] = Mult(n2, p2_0);
    rhs2[2] = Mult(n3, help);
    a2[0][0] = n1[0]; a2[0][1] = n1[1]; a2[0][2] = n1[2];
    a2[1][0] = n2[0]; a2[1][1] = n2[1]; a2[1][2] = n2[2];
    a2[2][0] = n3[0]; a2[2][1] = n3[1]; a2[2][2] = n3[2];
    det2 = ABS(Det3d(a2));
    if(det2>SMALL)
    {
      InvMatMult3d(p2,rhs2,a2);
    }
  }


  n3[0] = sp[2] - ep[2];
  n3[1] = 0.0;
  n3[2] = ep[0] - sp[0];
  l = Lenght(n3);
  if(l>SMALL)
  {
    Scale(n3,1/l);
    help[0] = sp[0];
    help[1] = 0.0;
    help[2] = sp[2];
    rhs3[0] = Mult(n1, p1_0);
    rhs3[1] = Mult(n2, p2_0);
    rhs3[2] = Mult(n3, help);
    a3[0][0] = n1[0]; a3[0][1] = n1[1]; a3[0][2] = n1[2];
    a3[1][0] = n2[0]; a3[1][1] = n2[1]; a3[1][2] = n2[2];
    a3[2][0] = n3[0]; a3[2][1] = n3[1]; a3[2][2] = n3[2];
    det3 = ABS(Det3d(a3));
    if(det3>SMALL)
    {
      InvMatMult3d(p3,rhs3,a3);
    }
  }

  det = det1;
  VecCopy3d(p, p1);
  if(det2>det)
  {
    det = det2;
    VecCopy3d(p, p2);
  }
  if(det3>det)
  {
    det = det3;
    VecCopy3d(p, p3);
  }
  return(1);
}

static int Cut_Test_3d(double *p1_0, double *p1_1, double *p1_2,
                       double *p2_0, double *p2_1, double *p2_2)
{
  int i, count, comp;
  double n1[3], n2[3];
  double spl[3], epl[3], sp1[2][3], ep1[2][3], sp2[2][3], ep2[2][3];
  double p1[3], p2[3], p3[3], p4[3], sp[3], ep[3];
  double m[3], l, max, x11, x12, x21, x22;

  GetNormalVector(p1_0, p1_1, p1_2, n1);
  GetNormalVector(p2_0, p2_1, p2_2, n2);

  count = 0;
  if(Points_on_different_sides(n1, p1_0, p2_0, p2_1, spl, epl))
  {
    sp2[count][0] = spl[0]; sp2[count][1] = spl[1]; sp2[count][2] = spl[2];
    ep2[count][0] = epl[0]; ep2[count][1] = epl[1]; ep2[count][2] = epl[2];
    count++;
  }
  if(Points_on_different_sides(n1, p1_0, p2_1, p2_2, spl, epl))
  {
    sp2[count][0] = spl[0]; sp2[count][1] = spl[1]; sp2[count][2] = spl[2];
    ep2[count][0] = epl[0]; ep2[count][1] = epl[1]; ep2[count][2] = epl[2];
    count++;
  }
  if(Points_on_different_sides(n1, p1_0, p2_2, p2_0, spl, epl))
  {
    sp2[count][0] = spl[0]; sp2[count][1] = spl[1]; sp2[count][2] = spl[2];
    ep2[count][0] = epl[0]; ep2[count][1] = epl[1]; ep2[count][2] = epl[2];
    count++;
  }
  if(count!=2)
    printf("%s\n", "ERROR");

  count = 0;
  if(Points_on_different_sides(n2, p2_0, p1_0, p1_1, spl, epl))
  {
    sp1[count][0] = spl[0]; sp1[count][1] = spl[1]; sp1[count][2] = spl[2];
    ep1[count][0] = epl[0]; ep1[count][1] = epl[1]; ep1[count][2] = epl[2];
    count++;
  }
  if(Points_on_different_sides(n2, p2_0, p1_1, p1_2, spl, epl))
  {
    sp1[count][0] = spl[0]; sp1[count][1] = spl[1]; sp1[count][2] = spl[2];
    ep1[count][0] = epl[0]; ep1[count][1] = epl[1]; ep1[count][2] = epl[2];
    count++;
  }
  if(Points_on_different_sides(n2, p2_0, p1_2, p1_0, spl, epl))
  {
    sp1[count][0] = spl[0]; sp1[count][1] = spl[1]; sp1[count][2] = spl[2];
    ep1[count][0] = epl[0]; ep1[count][1] = epl[1]; ep1[count][2] = epl[2];
    count++;
  }
  if(count!=2)
    printf("%s\n", "ERROR");

  /* search crosspoint edge + 2 face planes */
  sp[0] = sp1[0][0];
  sp[1] = sp1[0][1];
  sp[2] = sp1[0][2];
  ep[0] = ep1[0][0];
  ep[1] = ep1[0][1];
  ep[2] = ep1[0][2];
  Calculate_Point(n1, p1_0, n2, p2_0, sp, ep, p1) ;

  sp[0] = sp1[1][0];
  sp[1] = sp1[1][1];
  sp[2] = sp1[1][2];
  ep[0] = ep1[1][0];
  ep[1] = ep1[1][1];
  ep[2] = ep1[1][2];
  Calculate_Point(n1, p1_0, n2, p2_0, sp, ep, p2) ;

  sp[0] = sp2[0][0];
  sp[1] = sp2[0][1];
  sp[2] = sp2[0][2];
  ep[0] = ep2[0][0];
  ep[1] = ep2[0][1];
  ep[2] = ep2[0][2];
  Calculate_Point(n1, p1_0, n2, p2_0, sp, ep, p3) ;

  sp[0] = sp2[1][0];
  sp[1] = sp2[1][1];
  sp[2] = sp2[1][2];
  ep[0] = ep2[1][0];
  ep[1] = ep2[1][1];
  ep[2] = ep2[1][2];
  Calculate_Point(n1, p1_0, n2, p2_0, sp, ep, p4) ;

  /* check, if intersection occur */
  Minus(m, p2, p1);
  l = Lenght(m);
  if(l>SMALL)
    Scale(m,1/l);

  max = 0.0;
  for(i=0; i<3; i++)
  {
    if(max<ABS(m[i]))
    {
      comp = i;
      max = ABS(m[i]);
    }
  }
  if(p1[comp]<p2[comp])
  {
    x11 = p1[comp];
    x12 = p2[comp];
  }
  else
  {
    x11 = p2[comp];
    x12 = p1[comp];
  }
  if(p3[comp]<p4[comp])
  {
    x21 = p3[comp];
    x22 = p4[comp];
  }
  else
  {
    x21 = p4[comp];
    x22 = p3[comp];
  }

  if( (x12<x21) || (x11>x22) )
    return(0);
  else
    return(1);
}

static int Three_Points_in_common(int *id1, int *id2)
{
  if( ( (id1[0]==id2[0])&&(id1[1]==id2[1])&&(id1[2]==id2[2]) )
      ||  ( (id1[0]==id2[1])&&(id1[1]==id2[2])&&(id1[2]==id2[0]) )
      ||  ( (id1[0]==id2[2])&&(id1[1]==id2[0])&&(id1[2]==id2[1]) )
      ||      ( (id1[1]==id2[0])&&(id1[0]==id2[1])&&(id1[2]==id2[2]) )
      ||  ( (id1[1]==id2[1])&&(id1[0]==id2[2])&&(id1[2]==id2[0]) )
      ||  ( (id1[1]==id2[2])&&(id1[0]==id2[0])&&(id1[2]==id2[1]) )    )
  {
    //		printf("%s\n", "Three_Points_in_common: gleiche Dreiecke");
    return(0);                          /* keine Ueberschneidung */
  }
  else
  {
    //		printf("%s\n", "Three_Points_in_common: Ueberschneidung ist moeglich");
    return(1);
  }
}

int Cross_Check(double *p1_0, double *p1_1, double *p1_2, int *id1,
                double *p2_0, double *p2_1, double *p2_2, int *id2)
{
  int flag, flag1;
  double p1[3][3], p2[3][3];

  flag = 0;
  if(Three_Points_in_common(id1, id2))
  {
    if(Bounding_Box_3d(p1_0, p1_1, p1_2, p2_0, p2_1, p2_2))
    {
      if(Faces_in_the_same_Plane(p1_0, p1_1, p1_2, p2_0, p2_1, p2_2, &flag1))
      {
        if(flag1==1)
        {
          if(In_Plane_Two_points_common(p1_0, p1_1, p1_2, id1, p2_0, p2_1, p2_2, id2, &flag1))
            flag = flag1;
          else
          if(In_Plane_One_point_common(p1_0, p1_1, p1_2, id1, p2_0, p2_1, p2_2, id2, &flag1))
            flag = flag1;
          else
          if(Cut_Test_2d(p1_0, p1_1, p1_2, p2_0, p2_1, p2_2))
            flag = 1;
          else
            flag = 0;
        }
        else
          flag = 0;
      }
      else
      {
        if(Two_points_common(id1, id2))
          flag = 0;
        else
        {
          if(One_point_common(p1_0, p1_1, p1_2, id1, p2_0, p2_1, p2_2, id2, p1, p2, &flag1))
            flag = flag1;
          else
          {
            if(Faces_on_the_same_side(p1_0, p1_1, p1_2, p2_0, p2_1, p2_2))
              flag = 0;
            else
            {
              if(Faces_nearly_on_the_same_side(p1_0, p1_1, p1_2, id1, p2_0, p2_1, p2_2, id2))
                flag = 0;
              else
              {
                /* richtiger Schnitt-Test */
                if(Cut_Test_3d(p1_0, p1_1, p1_2, p2_0, p2_1, p2_2))
                  flag = 1;
                else
                  flag = 0;
              }
            }
          }
        }
      }
    }
    else
      flag = 0;
  }
  return(flag);
}
