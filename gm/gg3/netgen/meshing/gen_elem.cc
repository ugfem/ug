// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <math.h>

#include <myadt.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>
#include <meshing/ruler3.hh>

const int npoint = 5;
const int aux_point = 0;

const double search_radius = 1.75;
const double min_angle = 0.53;
const double min_distance = 0.55;
static int DEBUG = 0;

double Sphere_radius(double *p,double *p_opt,double *D,double *M);
double Calc_D_M(double *p1,double *p2,double *p3,double *D,double *M);
extern double Delaunay(double *p1, double *p2, double *p3, double *p4, double *p5);
int Check_All_Points(ARRAY<Point3d> & points,Element element);
extern int Cross_Check( double *p1_0, double *p1_1, double *p1_2, int *id1,
                        double *p2_0, double *p2_1, double *p2_2, int *id2);


static double ABS(double a)
{
  if(a>=0.0)
    return(a);
  else
    return(-a);
}

static Vec3d Get_Normalvector_Triang(Point3d p1,Point3d p2,Point3d p3)
{
  Vec3d n,n1,n2;

  n1 = p2 - p1;
  n1 /= n1.Length();
  n2 = p3 - p1;
  n2 /= n2.Length();

  n = Cross(n1,n2);
  n /= n.Length();

  return(n);
}

static Vec3d Get_Normalvector_Quad(Point3d A,Point3d B,Point3d C,Point3d D)
{
  Vec3d a,b,n;

  a = (B - A) + (C - D);
  a /= a.Length();
  b = (D - A) + (C - B);
  b /= b.Length();
  n = Cross(a,b);
  n /= n.Length();
  return(n);
}

int Search_TriangFace(ARRAY<Element> &lfaces,int id1,int id2,int id3)
{
  int i, id;

  id = 0;
  for(i=1; i<=lfaces.Size(); i++)
  {
    if(lfaces.Get(i).NP()==3)
    {
      if( (lfaces.Get(i).PNum(1)==id1)&&(lfaces.Get(i).PNum(2)==id2)&&(lfaces.Get(i).PNum(3)==id3)
          || (lfaces.Get(i).PNum(1)==id2)&&(lfaces.Get(i).PNum(2)==id3)&&(lfaces.Get(i).PNum(3)==id1)
          || (lfaces.Get(i).PNum(1)==id3)&&(lfaces.Get(i).PNum(2)==id1)&&(lfaces.Get(i).PNum(3)==id2) )
        id = i;
    }
  }
  return(id);
}

int RemoveOrAdd_Triangle(       ARRAY<Element> & lfaces,
                                ARRAY<INDEX> & delfaces,
                                int id1,int id2,int id3)
{
  int i,flag;
  Element element;

  flag = Search_TriangFace(lfaces,id1,id2,id3);
  if(flag!=0)
  {
    delfaces.Append(flag);
    if(DEBUG) printf("%s %d %d %d\n","remove face",lfaces.Get(flag).PNum(1),lfaces.Get(flag).PNum(2),lfaces.Get(flag).PNum(3));
  }
  else
  {
    element.PNum(1) = id2;
    element.PNum(2) = id1;
    element.PNum(3) = id3;
    element.PNum(4) = 0;
    element.SetNP(3);
    lfaces.Append(element);
    if(DEBUG) printf("%s %d %d %d\n","add face",id1,id2,id3);
  }
  return(1);
}

int Search_QuadFace(ARRAY<Element> &lfaces,int id1,int id2,int id3,int id4)
{
  int i,id;

  id = 0;
  for(i=1; i<=lfaces.Size(); i++)
  {
    if(lfaces.Get(i).NP()==4)
    {
      if( (lfaces.Get(i).PNum(1)==id1)&&(lfaces.Get(i).PNum(2)==id2)&&(lfaces.Get(i).PNum(3)==id3)&&(lfaces.Get(i).PNum(4)==id4)
          || (lfaces.Get(i).PNum(1)==id2)&&(lfaces.Get(i).PNum(2)==id3)&&(lfaces.Get(i).PNum(3)==id4)&&(lfaces.Get(i).PNum(4)==id1)
          || (lfaces.Get(i).PNum(1)==id3)&&(lfaces.Get(i).PNum(2)==id4)&&(lfaces.Get(i).PNum(3)==id1)&&(lfaces.Get(i).PNum(4)==id2)
          || (lfaces.Get(i).PNum(1)==id4)&&(lfaces.Get(i).PNum(2)==id1)&&(lfaces.Get(i).PNum(3)==id2)&&(lfaces.Get(i).PNum(4)==id3) )
      {
        id = i;
      }
    }
  }
  return(id);
}

int RemoveOrAdd_Quadrilateral(  ARRAY<Element> & lfaces,
                                ARRAY<INDEX> & delfaces,
                                int id1,int id2,int id3,int id4)
{
  int i,flag;
  Element element;

  flag = Search_QuadFace(lfaces,id1,id2,id3,id4);
  if(flag!=0)
  {
    delfaces.Append(flag);
    if(DEBUG) printf("%s %d %d %d %d\n","remove face",lfaces.Get(flag).PNum(1),lfaces.Get(flag).PNum(2),lfaces.Get(flag).PNum(3),lfaces.Get(flag).PNum(4));
  }
  else
  {
    element.PNum(1) = id4;
    element.PNum(2) = id3;
    element.PNum(3) = id2;
    element.PNum(4) = id1;
    element.SetNP(4);
    lfaces.Append(element);
    if(DEBUG) printf("%s %d %d %d %d\n","add face",id1,id2,id3,id4);
  }
  return(1);
}

int Cut(int i1, int i2 ,int test_id, ARRAY<Point3d> & lpoints, ARRAY<Element> & lfaces, Point3d testpoint)
{
  int j,flag;
  FILE *file;
  int id1[3],id2[3];
  double p1_0[3], p1_1[3], p1_2[3], p2_0[3], p2_1[3], p2_2[3];

  flag = 0;
  id1[0] = lfaces.Get(1).PNum(i1);
  id1[1] = lfaces.Get(1).PNum(i2);
  id1[2] = test_id;
  p1_0[0] = lpoints[lfaces.Get(1).PNum(i1)].X();
  p1_0[1] = lpoints[lfaces.Get(1).PNum(i1)].Y();
  p1_0[2] = lpoints[lfaces.Get(1).PNum(i1)].Z();
  p1_1[0] = lpoints[lfaces.Get(1).PNum(i2)].X();
  p1_1[1] = lpoints[lfaces.Get(1).PNum(i2)].Y();
  p1_1[2] = lpoints[lfaces.Get(1).PNum(i2)].Z();
  p1_2[0] = testpoint.X();
  p1_2[1] = testpoint.Y();
  p1_2[2] = testpoint.Z();
  for(j=2; j<=lfaces.Size(); j++)
  {
    id2[0] = lfaces.Get(j).PNum(1);
    id2[1] = lfaces.Get(j).PNum(2);
    id2[2] = lfaces.Get(j).PNum(3);
    p2_0[0] = lpoints[id2[0]].X();
    p2_0[1] = lpoints[id2[0]].Y();
    p2_0[2] = lpoints[id2[0]].Z();
    p2_1[0] = lpoints[id2[1]].X();
    p2_1[1] = lpoints[id2[1]].Y();
    p2_1[2] = lpoints[id2[1]].Z();
    p2_2[0] = lpoints[id2[2]].X();
    p2_2[1] = lpoints[id2[2]].Y();
    p2_2[2] = lpoints[id2[2]].Z();

    if(Cross_Check(p1_0,p1_1,p1_2,id1,p2_0,p2_1,p2_2,id2))
      return(1);

    if(lfaces.Get(j).NP()==4)
    {
      id2[0] = lfaces.Get(j).PNum(1);
      id2[1] = lfaces.Get(j).PNum(3);
      id2[2] = lfaces.Get(j).PNum(4);
      p2_0[0] = lpoints[id2[0]].X();
      p2_0[1] = lpoints[id2[0]].Y();
      p2_0[2] = lpoints[id2[0]].Z();
      p2_1[0] = lpoints[id2[1]].X();
      p2_1[1] = lpoints[id2[1]].Y();
      p2_1[2] = lpoints[id2[1]].Z();
      p2_2[0] = lpoints[id2[2]].X();
      p2_2[1] = lpoints[id2[2]].Y();
      p2_2[2] = lpoints[id2[2]].Z();

      if(Cross_Check(p1_0,p1_1,p1_2,id1,p2_0,p2_1,p2_2,id2))
        return(1);
    }
  }

  return(flag);
}

int Save_local_Situation(       ARRAY<Point3d> & lpoints,
                                ARRAY<Element> & lfaces,
                                Element element)
{
  FILE *file;
  int i;

  file = fopen("grape", "w+");
  fprintf(file, "%d\n", lpoints.Size());
  for (i=1; i<=lpoints.Size(); i++)
    fprintf(file, "%lf %lf %lf\n", lpoints[i].X(), lpoints[i].Y(), lpoints[i].Z());
  fprintf(file, "%d\n", lfaces.Size());
  for (i=1; i<=lfaces.Size(); i++)
  {
    fprintf(file, "%d %d %d\n", lfaces.Get(i).PNum(1), lfaces.Get(i).PNum(2), lfaces.Get(i).PNum(3));
    if(lfaces.Get(i).NP()==4)
      fprintf(file, "%d %d %d\n", lfaces.Get(i).PNum(1), lfaces.Get(i).PNum(3), lfaces.Get(i).PNum(4));
  }
  if(element.NP()==4)
    fprintf(file, "%d %d %d %d\n",  element.PNum(1),
            element.PNum(2),
            element.PNum(3),
            element.PNum(4));


  fclose(file);
  return(0);
}



static Local_Out(ARRAY<Point3d> & lpoints,ARRAY<Element> & lfaces, ARRAY<int> & prism_flags)
{
  int i,j;

  printf("%s\n","points");
  for(i=1; i<=lpoints.Size(); i++)
    printf("%lf %lf %lf\n",lpoints[i].X(),lpoints[i].Y(),lpoints[i].Z());
  printf("%s\n","faces");
  for(i=1; i<=lfaces.Size(); i++)
  {
    for(j=1; j<=lfaces[i].NP()-1; j++)
      printf("%d %s",lfaces.Get(i).PNum(j),"- ");
    printf("%d %d\n",lfaces.Get(i).PNum(lfaces[i].NP()),prism_flags.Get(i));
  }
  return(1);
}

int Search_PrismTriangle(       ARRAY<Point3d> & lpoints,ARRAY<Element> & lfaces,
                                ARRAY<Point3d> & ASL,ARRAY<Point3d> & BSL,ARRAY<Point3d> & CSL,
                                ARRAY<int> & id_ASL,ARRAY<int> & id_BSL,ARRAY<int> & id_CSL,
                                ARRAY<int> & prism_flags,
                                int &id1,int &id2,int&id3)
{
  int i,id;

  id = 0;
  if( (ASL.Size()==1)&&(BSL.Size()==1)&&(CSL.Size()==1) )
  {
    for(i=1; i<=lfaces.Size(); i++)
    {
      if(lfaces.Get(i).NP()==3)
      {
        if( (lfaces.Get(i).PNum(1)==id_ASL.Get(1))&&(lfaces.Get(i).PNum(3)==id_BSL.Get(1))&&(lfaces.Get(i).PNum(2)==id_CSL.Get(1)) )
        {
          id = i;
          id1 = lfaces.Get(i).PNum(1);
          id2 = lfaces.Get(i).PNum(3);
          id3 = lfaces.Get(i).PNum(2);
        }
        if( (lfaces.Get(i).PNum(3)==id_ASL.Get(1))&&(lfaces.Get(i).PNum(2)==id_BSL.Get(1))&&(lfaces.Get(i).PNum(1)==id_CSL.Get(1)) )
        {
          id = i;
          id1 = lfaces.Get(i).PNum(3);
          id2 = lfaces.Get(i).PNum(2);
          id3 = lfaces.Get(i).PNum(1);
        }
        if( (lfaces.Get(i).PNum(2)==id_ASL.Get(1))&&(lfaces.Get(i).PNum(1)==id_BSL.Get(1))&&(lfaces.Get(i).PNum(3)==id_CSL.Get(1)) )
        {
          id = i;
          id1 = lfaces.Get(i).PNum(2);
          id2 = lfaces.Get(i).PNum(1);
          id3 = lfaces.Get(i).PNum(3);
        }
      }
    }
  }

  if(prism_flags.Get(1)==prism_flags.Get(id))
    return(id);                         // gefunden
  else
    return(-1);
}

int Free_Prism( ARRAY<Point3d> & lpoints,ARRAY<Element> & lfaces,Element & element,
                ARRAY<Point3d> & AL,ARRAY<Point3d> & BL,ARRAY<Point3d> & CL,
                ARRAY<int> & id_AL,ARRAY<int> & id_BL,ARRAY<int> & id_CL,
                ARRAY<Point3d> & ASL,ARRAY<Point3d> & BSL,ARRAY<Point3d> & CSL,
                ARRAY<int> & id_ASL,ARRAY<int> & id_BSL,ARRAY<int> & id_CSL,
                Point3d AN,Point3d BN,Point3d CN)
{
  int i;
  if( (AL.Size()==0)&&(BL.Size()==0)&&(CL.Size()==0) )
  {
    element.PNum(1) = 1;
    element.PNum(2) = 2;
    element.PNum(3) = 3;
    element.PNum(4) = lpoints.Size() + 1;
    element.PNum(5) = lpoints.Size() + 2;
    element.PNum(6) = lpoints.Size() + 3;
    element.SetNP(6);

    lpoints.Append(AN);
    lpoints.Append(BN);
    lpoints.Append(CN);

    return(1);
  }
  else
    return(0);
}

#define xy_dist 0.1

int Generate_Prism (ARRAY<Point3d> & lpoints, ARRAY<Element> & lfaces,
                    ARRAY<Element> & elements,
                    ARRAY<INDEX> & delfaces, ARRAY<int> & prism_flags)
{
  int i,j,flag,face_id, id1, id2, id3;
  double max_l,l,dA,dB,dC;
  Point3d A,B,C,AN,BN,CN,M,dummy;
  Vec3d a,b,c,m,n,dist_A,dist_B,dist_C;
  Element element;
  ARRAY<Point3d> AL;
  ARRAY<int> id_AL;
  ARRAY<Point3d> BL;
  ARRAY<int> id_BL;
  ARRAY<Point3d> CL;
  ARRAY<int> id_CL;
  ARRAY<Point3d> ASL;
  ARRAY<int> id_ASL;
  ARRAY<Point3d> BSL;
  ARRAY<int> id_BSL;
  ARRAY<Point3d> CSL;
  ARRAY<int> id_CSL;


  if(DEBUG) Local_Out(lpoints,lfaces,prism_flags);

  A = lpoints[lfaces.Get(1).PNum(1)];
  B = lpoints[lfaces.Get(1).PNum(2)];
  C = lpoints[lfaces.Get(1).PNum(3)];

  /* centroid of the face */
  M.X() = ( A.X() + B.X() + C.X() ) / 3;
  M.Y() = ( A.Y() + B.Y() + C.Y() ) / 3;
  M.Z() = ( A.Z() + B.Z() + C.Z() ) / 3;

  /* normalvector on the first face */
  n = Get_Normalvector_Triang(A,B,C);

  max_l = (B - A).Length();
  l =  (C - A).Length();
  if(l>max_l)
    max_l = l;
  l =  (C - B).Length();
  if(l>max_l)
    max_l = l;

  /* determine best new points */
  AN = A - n;
  BN = B - n;
  CN = C - n;

  AL.SetSize(0);
  id_AL.SetSize(0);
  BL.SetSize(0);
  id_BL.SetSize(0);
  BL.SetSize(0);
  id_BL.SetSize(0);

  for(i=1; i<=lpoints.Size(); i++)
  {
    if( (i!=lfaces.Get(1).PNum(1)) && (i!=lfaces.Get(1).PNum(2)) && (i!=lfaces.Get(1).PNum(3)) )
    {
      m = lpoints[i] - M;
      dummy = lpoints[i];
      m /= m.Length();
      if(m*n<-0.1)
      {
        dist_A = lpoints[i] - A;
        dist_B = lpoints[i] - B;
        dist_C = lpoints[i] - C;
        dA = dist_A.Length();
        dB = dist_B.Length();
        dC = dist_C.Length();
        if( (dA<dB)&&(dA<dC)&&(dA<search_radius) )
        {
          AL.Append(dummy);
          id_AL.Append(i);
        }
        else
        {
          if( (dB<dA)&&(dB<dC)&&(dB<search_radius) )
          {
            BL.Append(dummy);
            id_BL.Append(i);
          }
          else
          if( (dC<search_radius) )
          {
            CL.Append(dummy);
            id_CL.Append(i);
          }
        }
      }
    }
  }

  for(i=1; i<=AL.Size(); i++)
  {
    a = AL[i]-AN;
    a.Z() = 0.0;
    if(a.Length()<xy_dist)
    {
      ASL.Append(AL[i]);
      id_ASL.Append(id_AL[i]);
    }
  }
  for(i=1; i<=BL.Size(); i++)
  {
    b = BL[i]-BN;
    b.Z() = 0.0;
    if(b.Length()<xy_dist)
    {
      BSL.Append(BL[i]);
      id_BSL.Append(id_BL[i]);
    }
  }
  for(i=1; i<=CL.Size(); i++)
  {
    c = CL[i]-CN;
    c.Z() = 0.0;
    if(c.Length()<xy_dist)
    {
      CSL.Append(CL[i]);
      id_CSL.Append(id_CL[i]);
    }
  }

  /*	printf("%s\n","AL");
          for(i=1;i<=AL.Size();i++)
                  printf("%lf %lf %lf %d\n",AL[i].X(),AL[i].Y(),AL[i].Z(),id_AL[i]);
          printf("%s\n","BL");
          for(i=1;i<=BL.Size();i++)
                  printf("%lf %lf %lf %d\n",BL[i].X(),BL[i].Y(),BL[i].Z(),id_BL[i]);
          printf("%s\n","CL");
          for(i=1;i<=CL.Size();i++)
                  printf("%lf %lf %lf %d\n",CL[i].X(),CL[i].Y(),CL[i].Z(),id_CL[i]);
          printf("%s\n","################");
          printf("%s\n","ASL");
          for(i=1;i<=ASL.Size();i++)
                  printf("%lf %lf %lf %d\n",ASL[i].X(),ASL[i].Y(),ASL[i].Z(),id_ASL[i]);
          printf("%s\n","BSL");
          for(i=1;i<=BSL.Size();i++)
                  printf("%lf %lf %lf %d\n",BSL[i].X(),BSL[i].Y(),BSL[i].Z(),id_BSL[i]);
          printf("%s\n","CSL");
          for(i=1;i<=CSL.Size();i++)
                  printf("%lf %lf %lf %d\n",CSL[i].X(),CSL[i].Y(),CSL[i].Z(),id_CSL[i]);*/

  face_id = Search_PrismTriangle(lpoints,lfaces,ASL,BSL,CSL,id_ASL,id_BSL,id_CSL, prism_flags, id1, id2, id3);
  if(face_id==-1)
    return(0);
  else
  {
    // create element
    element.PNum(1) = 1;
    element.PNum(2) = 2;
    element.PNum(3) = 3;
    element.PNum(4) = id1;
    element.PNum(5) = id2;
    element.PNum(6) = id3;
    element.SetNP(6);

    // cross check

    // in_element check

    elements.Append (element);
    if(DEBUG) printf("%s %d %d %d %d %d %d\n","create element",element.PNum(1),element.PNum(2),element.PNum(3),element.PNum(4),
                     element.PNum(5),element.PNum(6));
    RemoveOrAdd_Triangle(lfaces,delfaces,element.PNum(5),element.PNum(4),element.PNum(6));
    RemoveOrAdd_Quadrilateral(lfaces,delfaces,element.PNum(4),element.PNum(5),element.PNum(2),element.PNum(1));
    RemoveOrAdd_Quadrilateral(lfaces,delfaces,element.PNum(5),element.PNum(6),element.PNum(3),element.PNum(2));
    RemoveOrAdd_Quadrilateral(lfaces,delfaces,element.PNum(6),element.PNum(4),element.PNum(1),element.PNum(3));
    delfaces.Append(1);
    return(1);

  }
}

int Order_List3(ARRAY<Point3d> & lpoints,
                ARRAY<Element> & lfaces,
                ARRAY<Point3d> &Q,
                Point3d P[1+npoint+aux_point],
                ARRAY<int> &id_Q,
                ARRAY<double> &distance,
                ARRAY<Point3d> &OL,
                ARRAY<int> &id_OL)
{
  int i,j,help_id,min_j;
  double p1[3],p2[3],p3[3],p4[3],p5[3],helpx,helpy,helpz,min,help,D[3],M[3],p_opt[3];
  Point3d dummy;
  double distance1[20];
  dummy = P[0];
  Q.Append(dummy);
  id_Q.Append(lpoints.Size()+1);

  p1[0] = lpoints[lfaces.Get(1).PNum(1)].X();
  p1[1] = lpoints[lfaces.Get(1).PNum(1)].Y();
  p1[2] = lpoints[lfaces.Get(1).PNum(1)].Z();

  p2[0] = lpoints[lfaces.Get(1).PNum(2)].X();
  p2[1] = lpoints[lfaces.Get(1).PNum(2)].Y();
  p2[2] = lpoints[lfaces.Get(1).PNum(2)].Z();

  p3[0] = lpoints[lfaces.Get(1).PNum(3)].X();
  p3[1] = lpoints[lfaces.Get(1).PNum(3)].Y();
  p3[2] = lpoints[lfaces.Get(1).PNum(3)].Z();

  Calc_D_M(p1,p2,p3,D,M);

  p_opt[0] = dummy.X();
  p_opt[1] = dummy.Y();
  p_opt[2] = dummy.Z();

  for(i=0; i<20; i++)
    distance1[i] = -1.0;
  for(i=1; i<=Q.Size(); i++)
  {
    p4[0] = Q[i].X();
    p4[1] = Q[i].Y();
    p4[2] = Q[i].Z();

    distance[i] = Sphere_radius(p4,p_opt,D,M);
    distance1[i] = Delaunay(p1,p2,p3,p4,p5);
    if(i==Q.Size())
      distance[i] = distance[i] * 1.65;
  }

  /* order  Q's */
  for(i=1; i<=Q.Size(); i++)
  {
    min = 1000000.0;
    for(j=i; j<=Q.Size(); j++)
    {
      if(min>distance[j])
      {
        min_j = j;
        min = distance[j];
      }
    }
    /* vertauschen */
    help = distance[i];
    help_id = id_Q[i];
    helpx = Q[i].X();
    helpy = Q[i].Y();
    helpz = Q[i].Z();

    distance[i] = distance[min_j];
    id_Q[i] = id_Q[min_j];
    Q[i].X() = Q[min_j].X();
    Q[i].Y() = Q[min_j].Y();
    Q[i].Z() = Q[min_j].Z();

    distance[min_j] = help;
    id_Q[min_j] = help_id;
    Q[min_j].X() = helpx;
    Q[min_j].Y() = helpy;
    Q[min_j].Z() = helpz;
  }

  for(i=1; i<npoint+aux_point; i++)
  {
    dummy = P[i];
    Q.Append(dummy);
    id_Q.Append(lpoints.Size()+1);
  }
  for(i=1; i<=Q.Size(); i++)
  {
    OL.Append(Q[i]);
    id_OL.Append(id_Q[i]);
  }

  if(DEBUG)
  {
    for(i=1; i<=OL.Size(); i++)
      printf("%lf %lf %lf %d %lf %lf\n",OL[i].X(),OL[i].Y(),OL[i].Z(),id_OL[i],distance[i],distance1[i]);
    printf("\n");
  }
  return(1);
}


int Generate_Pyramid (  ARRAY<Point3d> & lpoints, ARRAY<Element> & lfaces,
                        ARRAY<Element> & elements,
                        ARRAY<INDEX> & delfaces, ARRAY<int> & prism_flags)
{
  int i,j,flag1,flag2,flag3,flag4,new_id,ok,save,test_id;
  double l,max_l,dist;
  Point3d A,B,C,D,M,dummy,P[1+npoint+aux_point],testpoint;
  Vec3d a,b,n,m,dist_vec;
  Element element;
  ARRAY<Point3d> Q;
  ARRAY<int> id_Q;
  ARRAY<Point3d> OL;
  ARRAY<int> id_OL;
  ARRAY<Point3d> npoints;
  ARRAY<double> distance;
  ARRAY<Element> elem;

  if(DEBUG) Local_Out(lpoints,lfaces,prism_flags);

  A = lpoints[lfaces.Get(1).PNum(1)];
  B = lpoints[lfaces.Get(1).PNum(2)];
  C = lpoints[lfaces.Get(1).PNum(3)];
  D = lpoints[lfaces.Get(1).PNum(4)];

  /* centroid of the face */
  M.X() = ( A.X() + B.X() + C.X() + D.X() ) / 4;
  M.Y() = ( A.Y() + B.Y() + C.Y() + D.Y() ) / 4;
  M.Z() = ( A.Z() + B.Z() + C.Z() + D.Z() ) / 4;

  /* normalvector on the first face */
  n = Get_Normalvector_Quad(A,B,C,D);

  max_l = (B - A).Length();
  l =  (C - A).Length();
  if(l>max_l)
    max_l = l;
  l =  (C - B).Length();
  if(l>max_l)
    max_l = l;

  /* determine inner points */
  for(i=0; i<npoint; i++)
  {
    P[i].X() = M.X() - 1.0 * ((double)(npoint-i) / (double)(npoint)) * n.X()/1.0;
    P[i].Y() = M.Y() - 1.0 * ((double)(npoint-i) / (double)(npoint)) * n.Y()/1.0;
    P[i].Z() = M.Z() - 1.0 * ((double)(npoint-i) / (double)(npoint)) * n.Z()/1.0;
  }

  for(i=0; i<aux_point; i++)
  {
    P[npoint+i].X() = (M.X()+P[npoint-1+i].X())/2;
    P[npoint+i].Y() = (M.Y()+P[npoint-1+i].Y())/2;
    P[npoint+i].Z() = (M.Z()+P[npoint-1+i].Z())/2;
  }

  for(i=1; i<=lpoints.Size(); i++)
  {
    if( (i!=lfaces.Get(1).PNum(1)) && (i!=lfaces.Get(1).PNum(2)) && (i!=lfaces.Get(1).PNum(3)) && (i!=lfaces.Get(1).PNum(4)) )
    {
      m = lpoints[i] - M;
      dummy = lpoints[i];
      m /= m.Length();
      if(m*n<0.0)
      {
        dist_vec = lpoints[i] - P[0];
        dist = dist_vec.Length();
        if((dist<search_radius))
        {
          Q.Append(dummy);
          id_Q.Append(i);
        }
      }
    }
  }

  OL.SetSize(0);
  id_OL.SetSize(0);
  distance.SetSize(Q.Size()+npoint+aux_point+1);
  for(i=1; i<=Q.Size()+npoint+aux_point; i++)
    distance[i]=-1.0;

  Order_List3(lpoints,lfaces,Q,P,id_Q,distance,OL,id_OL);

  new_id = -1;
  ok = 0;
  elem.SetSize(0);

  for(i=1; i<=OL.Size(); i++)
  {
    test_id = id_OL[i];
    testpoint.X() = OL[i].X();
    testpoint.Y() = OL[i].Y();
    testpoint.Z() = OL[i].Z();
    flag1 = Cut(1,2,test_id,lpoints,lfaces,testpoint);
    flag2 = Cut(2,3,test_id,lpoints,lfaces,testpoint);
    flag3 = Cut(3,4,test_id,lpoints,lfaces,testpoint);
    flag4 = Cut(4,1,test_id,lpoints,lfaces,testpoint);
    save = 0;
    if( (flag1==0)&&(flag2==0)&&(flag3==0)&&(flag4==0) )
    {
      new_id = id_OL[i];
      for(j=1; j<=lpoints.Size(); j++)
        npoints.Append(lpoints[i]);
      if(new_id>lpoints.Size())
        npoints.Append(OL[i]);
      element.PNum(1) = lfaces.Get(1).PNum(1);
      element.PNum(2) = lfaces.Get(1).PNum(2);
      element.PNum(3) = lfaces.Get(1).PNum(3);
      element.PNum(4) = lfaces.Get(1).PNum(4);
      element.PNum(5) = new_id;
      element.SetNP(5);
      if(new_id>lpoints.Size())
        lpoints.Append(OL[i]);
      break;
    }
  }

  if(DEBUG)
    Save_local_Situation(lpoints,lfaces,element);

  /* create Element */
  if(new_id==-1)
  {
    return(0);
  }
  if(DEBUG) printf("%s %d %d %d %d %d\n","create element",        element.PNum(1),element.PNum(2),element.PNum(3),element.PNum(4),element.PNum(5));
  elements.Append (element);

  /* delfaces and lfaces */
  delfaces.Append(1);

  RemoveOrAdd_Triangle(lfaces,delfaces,element.PNum(5),element.PNum(2),element.PNum(1));
  RemoveOrAdd_Triangle(lfaces,delfaces,element.PNum(5),element.PNum(3),element.PNum(2));
  RemoveOrAdd_Triangle(lfaces,delfaces,element.PNum(5),element.PNum(4),element.PNum(3));
  RemoveOrAdd_Triangle(lfaces,delfaces,element.PNum(5),element.PNum(1),element.PNum(4));

  return(1);
}
