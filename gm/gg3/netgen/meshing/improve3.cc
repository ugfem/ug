// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include <myadt.hh>

#include <geom/geom2d.hh>
#include <geom/geom3d.hh>


#include <linalg/linalg.hh>

#include <meshing/global.hh>
#include <meshing/meshing3.hh>
// #include <meshing/meshtool.hh>

#ifdef SOLIDGEOM
#include <geom/solid.hh>
#include <meshing/specpoin.hh>
#include <meshing/edgeflw.hh>
#endif

#include <opti/opti.hh>


#ifdef SOLIDGEOM
extern ARRAY<Surface*> surfaces;
#endif


static double CalcVolume (const ARRAY<Point3d> & points,
                          const ARRAY<Element> & elements)
{
  double vol;
  int i;
  Vec3d v1, v2, v3;

  vol = 0;
  for (i = 1; i <= elements.Size(); i++)
  {
    v1 = points.Get(elements.Get(i).PNum(2)) -
         points.Get(elements.Get(i).PNum(1));
    v2 = points.Get(elements.Get(i).PNum(3)) -
         points.Get(elements.Get(i).PNum(1));
    v3 = points.Get(elements.Get(i).PNum(4)) -
         points.Get(elements.Get(i).PNum(1));
    vol -= (Cross (v1, v2) * v3) / 6;
  }
  return vol;
}

static double CalcTetBadness (const Point3d & p1, const Point3d & p2,
                              const Point3d & p3, const Point3d & p4, double h)
{
  double vol, l, l4, l5, l6;
  double err;

  Vec3d v1 = p2 - p1;
  Vec3d v2 = p3 - p1;
  Vec3d v3 = p4 - p1;

  vol = - (Cross (v1, v2) * v3) / 6;
  l4 = Dist (p2, p3);
  l5 = Dist (p2, p4);
  l6 = Dist (p3, p4);

  l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

  if (vol < 1e-12) return 1e14;
  err = (l*l*l) / (1832.82 * vol);    // 6^4 * sqrt(2)

  if (h > 0)
    err += l / h +
           h * ( 1 / v1.Length() + 1 / v2.Length() + 1 / v3.Length() +
                 1 / l4 + 1 / l5 + 1 / l6) - 12;
  return err;
}


class PointFunction
{
public:
  ARRAY<Point3d> & points;
  const ARRAY<Element> & elements;
  TABLE<INDEX> elementsonpoint;
  int actpind;
  double h;

public:
  PointFunction (ARRAY<Point3d> & apoints,
                 const ARRAY<Element> & aelements);

  void SetPointIndex (int aactpind) { actpind = aactpind; }
  void SetLocalH (double ah) { h = ah; }
  double PointFunctionValue (const Point3d & pp) const;
  double PointFunctionValueGrad (const Point3d & pp, Vector & grad) const;
};


PointFunction :: PointFunction (ARRAY<Point3d> & apoints,
                                const ARRAY<Element> & aelements)
  : points(apoints), elements(aelements), elementsonpoint(apoints.Size())
{
  INDEX i;
  int j;

  for (i = 1; i <= elements.Size(); i++)
    for (j = 1; j <= elements[i].NP(); j++)
      elementsonpoint.Add (elements[i].PNum(j), i);
}

double PointFunction :: PointFunctionValue (const Point3d & pp) const
{
  int j;
  INDEX eli;
  const Element * el;
  double badness;
  ARRAY<const Point3d*> p(4);
  Point3d hp;

  badness = 0;

  hp = points.Elem(actpind);
  points.Elem(actpind) = pp;

  for (j = 1; j <= elementsonpoint.EntrySize(actpind); j++)
  {
    eli = elementsonpoint.Get(actpind, j);
    el = &elements.Get(eli);
    badness += CalcTetBadness (points.Get(el->PNum(1)),
                               points.Get(el->PNum(2)), points.Get(el->PNum(3)),
                               points.Get(el->PNum(4)), h);
  }

  points.Elem(actpind) = hp;
  return badness;
}


double PointFunction :: PointFunctionValueGrad (const Point3d & pp, Vector & grad) const
{
  double f, fr, fl, delta = 1e-6;
  Point3d hpp;

  f = PointFunctionValue (pp);

  hpp = pp;
  hpp.X() = pp.X() + delta;
  fr = PointFunctionValue (hpp);
  hpp.X() = pp.X() - delta;
  fl = PointFunctionValue (hpp);
  grad.Elem(1) = (fr - fl) / (2 * delta);

  hpp = pp;
  hpp.Y() = pp.Y() + delta;
  fr = PointFunctionValue (hpp);
  hpp.Y() = pp.Y() - delta;
  fl = PointFunctionValue (hpp);
  grad.Elem(2) = (fr - fl) / (2 * delta);

  hpp = pp;
  hpp.Z() = pp.Z() + delta;
  fr = PointFunctionValue (hpp);
  hpp.Z() = pp.Z() - delta;
  fl = PointFunctionValue (hpp);
  grad.Elem(3) = (fr - fl) / (2 * delta);

  return f;
}




class Opti3FreeMinFunction : public MinFunction
{
  const PointFunction & pf;
  Point3d sp1;

public:
  Opti3FreeMinFunction (const PointFunction & apf);
  void SetPoint (const Point3d & asp1) { sp1 = asp1; }
  virtual double FuncGrad (const Vector & x, Vector & g) const;
};


Opti3FreeMinFunction :: Opti3FreeMinFunction (const PointFunction & apf)
  : pf(apf)
{
  ;
}

double Opti3FreeMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
{
  Point3d pp;
  pp.X() = sp1.X() + x.Get(1);
  pp.Y() = sp1.Y() + x.Get(2);
  pp.Z() = sp1.Z() + x.Get(3);

  return pf.PointFunctionValueGrad (pp, grad);
}




#ifdef SOLIDGEOM
class Opti3SurfaceMinFunction : public MinFunction
{
  const PointFunction & pf;
  Point3d sp1;
  int surf;
  Vec3d t1, t2;

public:
  Opti3SurfaceMinFunction (const PointFunction & apf);

  void SetPoint (int asurf, const Point3d & asp1);

  void CalcNewPoint (const Vector & x, Point3d & np) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
};


Opti3SurfaceMinFunction :: Opti3SurfaceMinFunction (const PointFunction & apf)
  : MinFunction(), pf(apf)
{
  ;
}

void Opti3SurfaceMinFunction :: SetPoint (int asurf, const Point3d & asp1)
{
  Vec3d n;
  sp1 = asp1;
  surf = asurf;

  surfaces[surf] -> GetNormalVector (sp1, n);
  n.GetNormal (t1);
  t1 /= t1.Length();
  t2 = Cross (n, t1);
}


void Opti3SurfaceMinFunction :: CalcNewPoint (const Vector & x,
                                              Point3d & np) const
{
  np.X() = sp1.X() + x.Get(1) * t1.X() + x.Get(2) * t2.X();
  np.Y() = sp1.Y() + x.Get(1) * t1.Y() + x.Get(2) * t2.Y();
  np.Z() = sp1.Z() + x.Get(1) * t1.Z() + x.Get(2) * t2.Z();

  surfaces[surf] -> Project (np);
}

double Opti3SurfaceMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
{
  Vec3d n, vgrad;
  Point3d pp1;
  double badness;
  static Vector freegrad(3);

  CalcNewPoint (x, pp1);

  badness = pf.PointFunctionValueGrad (pp1, freegrad);
  vgrad.X() = freegrad.Get(1);
  vgrad.Y() = freegrad.Get(2);
  vgrad.Z() = freegrad.Get(3);

  surfaces[surf] -> GetNormalVector (pp1, n);

  vgrad -= (vgrad * n) * n;

  grad.Elem(1) = vgrad * t1;
  grad.Elem(2) = vgrad * t2;

  return badness;
}
#endif








#ifdef SOLIDGEOM
class Opti3EdgeMinFunction : public MinFunction
{
  const PointFunction & pf;
  Point3d sp1;
  int surf1, surf2;
  Vec3d t1;

public:
  Opti3EdgeMinFunction (const PointFunction & apf);

  void SetPoint (int asurf1, int asurf2, const Point3d & asp1);
  void CalcNewPoint (const Vector & x, Point3d & np) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
};

Opti3EdgeMinFunction :: Opti3EdgeMinFunction (const PointFunction & apf)
  : MinFunction(), pf(apf)
{
  ;
}

void Opti3EdgeMinFunction :: SetPoint (int asurf1, int asurf2,
                                       const Point3d & asp1)
{
  Vec3d n1, n2;
  sp1 = asp1;
  surf1 = asurf1;
  surf2 = asurf2;

  surfaces[surf1] -> GetNormalVector (sp1, n1);
  surfaces[surf2] -> GetNormalVector (sp1, n2);
  t1 = Cross (n1, n2);
}

void Opti3EdgeMinFunction :: CalcNewPoint (const Vector & x,
                                           Point3d & np) const
{
  np.X() = sp1.X() + x.Get(1) * t1.X();
  np.Y() = sp1.Y() + x.Get(1) * t1.Y();
  np.Z() = sp1.Z() + x.Get(1) * t1.Z();
  ProjectToEdge ( surfaces[surf1], surfaces[surf2], np);
}


double Opti3EdgeMinFunction :: FuncGrad (const Vector & x, Vector & grad) const
{
  Vec3d n1, n2, v1, vgrad;
  Point3d pp1;
  double badness;
  static Vector freegrad(3);

  CalcNewPoint (x, pp1);


  badness = pf.PointFunctionValueGrad (pp1, freegrad);

  vgrad.X() = freegrad.Get(1);
  vgrad.Y() = freegrad.Get(2);
  vgrad.Z() = freegrad.Get(3);

  surfaces[surf1] -> GetNormalVector (pp1, n1);
  surfaces[surf2] -> GetNormalVector (pp1, n2);

  v1 = Cross (n1, n2);
  v1 /= v1.Length();

  grad.Elem(1) = (vgrad * v1) * (t1 * v1);
  return badness;
}
#endif






#ifdef SOLIDGEOM
void Meshing3 :: ImproveMesh (ARRAY<Point3d> & points,
                              const ARRAY<Element> & surfelements,
                              const ARRAY<Element> & elements, double h)
{
  INDEX i, eli;
  int j;

  TABLE<INDEX> surfelementsonpoint(points.Size());
  Vector x(3), xsurf(2), xedge(1);
  int surf, surf1, surf2, surf3;

  (*testout).precision(8);
  // cout << "ImproveMesh" << endl;
  // cout << "Vol = " << CalcVolume (points, elements) << endl;


  for (i = 1; i <= surfelements.Size(); i++)
    for (j = 1; j <= surfelements[i].NP(); j++)
      surfelementsonpoint.Add (surfelements[i].PNum(j), i);


  PointFunction pf(points, elements);
  pf.SetLocalH (h);

  Opti3FreeMinFunction freeminf(pf);
  Opti3SurfaceMinFunction surfminf(pf);
  Opti3EdgeMinFunction edgeminf(pf);

  for (i = 1; i <= points.Size(); i++)
  {
    // cout << ".";

    //    (*testout) << "Now point " << i << endl;
    //    (*testout) << "Old: " << points.Get(i) << endl;

    pf.SetPointIndex (i);
    surf1 = surf2 = surf3 = 0;

    for (j = 1; j <= surfelementsonpoint.EntrySize(i); j++)
    {
      eli = surfelementsonpoint.Get(i, j);
      surf = surfelements.Get(eli).SurfaceIndex();

      if (!surf1)
        surf1 = surf;
      else if (surf1 != surf)
      {
        if (!surf2)
          surf2 = surf;
        else if (surf2 != surf)
          surf3 = surf;
      }
    }


    if (surf2 && !surf3)
    {
      //      (*testout) << "On Edge" << endl;
      xedge = 0;
      edgeminf.SetPoint (surf1, surf2, points.Elem(i));
      BFGS (xedge, edgeminf);

      edgeminf.CalcNewPoint (xedge, points.Elem(i));
    }

    if (surf1 && !surf2)
    {
      //      (*testout) << "In Surface" << endl;
      xsurf = 0;
      surfminf.SetPoint (surf1, points.Get(i));
      BFGS (xsurf, surfminf);

      surfminf.CalcNewPoint (xsurf, points.Elem(i));
    }

    if (!surf1)
    {
      //      (*testout) << "In Volume" << endl;
      x = 0;
      freeminf.SetPoint (points.Elem(i));
      BFGS (x, freeminf);

      points.Elem(i).X() += x.Get(1);
      points.Elem(i).Y() += x.Get(2);
      points.Elem(i).Z() += x.Get(3);
    }

    //    (*testout) << "New Point: " << points.Elem(i) << endl << endl;

  }
  // cout << endl;

  // cout << "Vol = " << CalcVolume (points, elements) << endl;
}
#endif


void Meshing3 :: ImproveMesh (ARRAY<Point3d> & points,
                              const ARRAY<Element> & elements,
                              int nboundnodes, double h)
{
  INDEX i;

  (*testout) << "Improve Mesh" << endl;
  // cout << "ImproveMesh" << endl;

  Vector x(3);

  (*testout).precision(8);

  PointFunction pf(points, elements);
  pf.SetLocalH (h);

  Opti3FreeMinFunction freeminf(pf);

  for (i = nboundnodes+1; i <= points.Size(); i++)
  {
    // cout << "." << flush;
    //    if (elementsonpoint.EntrySize(i) == 0) continue;

    freeminf.SetPoint (points.Elem(i));
    pf.SetPointIndex (i);

    x = 0;
    BFGS (x, freeminf);

    points.Elem(i).X() += x.Get(1);
    points.Elem(i).Y() += x.Get(2);
    points.Elem(i).Z() += x.Get(3);
  }
  // cout << endl;
}



static double CalcBad (const ARRAY<Point3d> & points,
                       const Element & elem, double h)
{
  return CalcTetBadness (
           points.Get(elem.PNum(1)),
           points.Get(elem.PNum(2)),
           points.Get(elem.PNum(3)),
           points.Get(elem.PNum(4)), h);
}


void Meshing3 :: CombineImprove (ARRAY<Point3d> & points,
                                 ARRAY<Element> & elements,
                                 int nboundnodes, double /* h */)
{
  int i, j, k, l;
  Point3d p1, p2, pnew;
  int pi1, pi2;
  Element * elem;
  double bad1, bad2;

  TABLE<int> elementsonnode(points.Size());
  ARRAY<Element*> hasonepoint;
  ARRAY<Element*> hasbothpoints;

  // cout << "CombineImprove" << endl;
  // cout << "nbp = " << nboundnodes << endl;

  // cout << "Vol = " << CalcVolume (points, elements) << endl;

  for (i = 1; i <= elements.Size(); i++)
    for (j = 1; j <= elements.Get(i).NP(); j++)
      elementsonnode.Add (elements.Get(i).PNum(j), i);

  bad1 = 0;
  for (i = 1; i <= elements.Size(); i++)
    bad1 += CalcBad (points, elements.Get(i), 0);
  // cout << "Total badness = " << bad1 << endl;

  for (i = 1; i <= elements.Size(); i++)
  {
    for (j = 1; j <= 6; j++)
    {
      elem = &elements.Elem(i);
      if (elem->PNum(1) == 0) continue;

      switch (j)
      {
      case 1 : pi1 = elem->PNum(1); pi2 = elem->PNum(2); break;
      case 2 : pi1 = elem->PNum(1); pi2 = elem->PNum(3); break;
      case 3 : pi1 = elem->PNum(1); pi2 = elem->PNum(4); break;
      case 4 : pi1 = elem->PNum(2); pi2 = elem->PNum(3); break;
      case 5 : pi1 = elem->PNum(2); pi2 = elem->PNum(4); break;
      case 6 : pi1 = elem->PNum(3); pi2 = elem->PNum(4); break;
      }

      if (pi2 < pi1) swap (pi1, pi2);

      hasonepoint.SetSize(0);
      hasbothpoints.SetSize(0);

      for (k = 1; k <= elementsonnode.EntrySize(pi1); k++)
      {
        elem = &elements[elementsonnode.Get (pi1, k)];
        if (elem->PNum(1) > 0)
        {
          if (elem->PNum(1) == pi2 || elem->PNum(2) == pi2 ||
              elem->PNum(3) == pi2 || elem->PNum(4) == pi2)
            hasbothpoints.Append (elem);
          else
            hasonepoint.Append (elem);
        }
      }

      for (k = 1; k <= elementsonnode.EntrySize(pi2); k++)
      {
        elem = &elements[elementsonnode.Get (pi2, k)];
        if (elem->PNum(1) > 0)
        {
          if (elem->PNum(1) == pi1 || elem->PNum(2) == pi1 ||
              elem->PNum(3) == pi1 || elem->PNum(4) == pi1)
            ;
          else
            hasonepoint.Append (elem);
        }
      }

      bad1 = 0;
      for (k = 1; k <= hasonepoint.Size(); k++)
      {
        elem = hasonepoint.Get(k);
        bad1 += CalcTetBadness (
          points.Get(elem->PNum(1)),
          points.Get(elem->PNum(2)),
          points.Get(elem->PNum(3)),
          points.Get(elem->PNum(4)), 0);
      }

      for (k = 1; k <= hasbothpoints.Size(); k++)
      {
        elem = hasbothpoints.Get(k);
        bad1 += CalcTetBadness (
          points.Get(elem->PNum(1)),
          points.Get(elem->PNum(2)),
          points.Get(elem->PNum(3)),
          points.Get(elem->PNum(4)), 0);
      }


      p1 = points[pi1];
      p2 = points[pi2];

      if (pi2 <= nboundnodes)
        continue;

      if (pi1 <= nboundnodes)
        pnew = p1;
      else
        pnew = Center (p1, p2);

      points[pi1] = pnew;
      points[pi2] = pnew;

      bad2 = 0;
      for (k = 1; k <= hasonepoint.Size(); k++)
      {
        elem = hasonepoint.Get(k);
        bad2 += CalcTetBadness (
          points.Get(elem->PNum(1)),
          points.Get(elem->PNum(2)),
          points.Get(elem->PNum(3)),
          points.Get(elem->PNum(4)), 0);
      }

      points[pi1] = p1;
      points[pi2] = p2;

      if (bad2 / hasonepoint.Size()  <
          bad1 / (hasonepoint.Size()+hasbothpoints.Size())  )
      {
        points[pi1] = pnew;
        // cout << "Connect point " << pi2 << " to " << pi1 << endl;
        for (k = 1; k <= elementsonnode.EntrySize(pi2); k++)
        {
          elem = &elements[elementsonnode.Get (pi2, k)];
          elementsonnode.Add (pi1, elementsonnode.Get (pi2, k));
          for (l = 1; l <= elem->NP(); l++)
            if (elem->PNum(l) == pi2)
              elem->PNum(l) = pi1;
        }

        for (k = 1; k <= hasbothpoints.Size(); k++)
          for (l = 1; l <= 4; l++)
            hasbothpoints[k] -> PNum(l) = 0;

        // cout << "bad1 = " << bad1 << " bad2 = " << bad2 << endl;
      }
    }
  }

  // cout << "Now deleting elements" << endl;

  for (k = elements.Size(); k >= 1; k--)
    if (elements[k].PNum(1) == 0)
      elements.DeleteElement(k);

  // cout << "Elements now: " << elements.Size() << endl;
  // cout << "Vol = " << CalcVolume (points, elements) << endl;
}


void Meshing3 :: SplitImprove (ARRAY<Point3d> & points,
                               ARRAY<Element> & elements, ARRAY<Element> & surfelements,
                               double /* h */)
{
  int i, j, k, l;
  INDEX_2 i2;
  Point3d p1, p2, pnew;
  int elnr;
  int pi1, pi2;
  Element * elem;
  double bad1, bad2;
  int has1, has2;

  TABLE<int> elementsonnode(points.Size());
  ARRAY<int> hasbothpoints;
  INDEX_2_HASHTABLE<int> boundaryedges(surfelements.Size());


  // cout << "SplitImprove" << endl;
  // cout << "Vol = " << CalcVolume (points, elements) << endl;

  bad1 = 0;
  for (i = 1; i <= elements.Size(); i++)
    bad1 += CalcBad (points, elements.Get(i), 0);
  // cout << "Total badness = " << bad1 << endl;

  for (i = 1; i <= surfelements.Size(); i++)
    for (j = 1; j <= 3; j++)
    {
      i2.I1() = surfelements.Get(i).PNumMod(j);
      i2.I2() = surfelements.Get(i).PNumMod(j+1);
      if (i2.I1() > i2.I2()) swap (i2.I1(), i2.I2());
      boundaryedges.Set (i2, 1);
    }

  for (i = 1; i <= elements.Size(); i++)
    for (j = 1; j <= elements.Get(i).NP(); j++)
      elementsonnode.Add (elements.Get(i).PNum(j), i);

  for (i = 1; i <= elements.Size(); i++)
  {
    for (j = 1; j <= 6; j++)
    {
      elem = &elements.Elem(i);
      switch (j)
      {
      case 1 : pi1 = elem->PNum(1); pi2 = elem->PNum(2); break;
      case 2 : pi1 = elem->PNum(1); pi2 = elem->PNum(3); break;
      case 3 : pi1 = elem->PNum(1); pi2 = elem->PNum(4); break;
      case 4 : pi1 = elem->PNum(2); pi2 = elem->PNum(3); break;
      case 5 : pi1 = elem->PNum(2); pi2 = elem->PNum(4); break;
      case 6 : pi1 = elem->PNum(3); pi2 = elem->PNum(4); break;
      }

      if (pi2 < pi1) swap (pi1, pi2);

      if (pi2 > elementsonnode.Size()) continue;

      i2.I1() = pi1;
      i2.I2() = pi2;
      if (boundaryedges.Used (i2)) continue;


      hasbothpoints.SetSize (0);
      for (k = 1; k <= elementsonnode.EntrySize(pi1); k++)
      {
        has1 = has2 = 0;
        elnr = elementsonnode.Get(pi1, k);
        elem = &elements.Elem (elnr);
        for (l = 1; l <= elem->NP(); l++)
        {
          if (elem->PNum(l) == pi1) has1 = 1;
          if (elem->PNum(l) == pi2) has2 = 1;
        }
        if (has1 && has2)
        {   // only once
          for (l = 1; l <= hasbothpoints.Size(); l++)
            if (hasbothpoints.Get(l) == elnr)
              has1 = 0;

          if (has1)
            hasbothpoints.Append (elnr);
        }
      }


      bad1 = 0;
      for (k = 1; k <= hasbothpoints.Size(); k++)
        bad1 += CalcBad (points, elements.Get(hasbothpoints.Get(k)), 0);

      p1 = points.Get(pi1);
      p2 = points.Get(pi2);
      pnew = Center (p1, p2);

      points.Elem(pi1) = pnew;
      bad2 = 0;
      for (k = 1; k <= hasbothpoints.Size(); k++)
        bad2 += CalcBad (points, elements.Get(hasbothpoints.Get(k)), 0);
      points.Elem(pi1) = p1;

      points.Elem(pi2) = pnew;
      for (k = 1; k <= hasbothpoints.Size(); k++)
        bad2 += CalcBad (points, elements.Get(hasbothpoints.Get(k)), 0);
      points.Elem(pi2) = p2;

      if (bad2 < bad1)
      {
        // cout << "Split Edge" << endl;
      }
    }
  }
  // cout << "Splitt - Improve done" << endl;
}





void Meshing3 :: SwapImprove (ARRAY<Point3d> & points,
                              ARRAY<Element> & elements, ARRAY<Element> & surfelements,
                              double /* h */)
{
  int i, j, k, l;
  INDEX_2 i2;
  //  Point3d p1, p2, pnew;
  int elnr;
  int pi1, pi2, pi3, pi4, pi5, pi6;
  Element el21(4), el22(4), el31(4), el32(4), el33(4);
  Element el1(4), el2(4), el3(4), el4(4);
  Element el1b(4), el2b(4), el3b(4), el4b(4);
  Element * elem;
  double bad1, bad2, bad3;
  int has1, has2;

  TABLE<int> elementsonnode(points.Size());
  ARRAY<int> hasbothpoints;
  INDEX_2_HASHTABLE<int> boundaryedges(surfelements.Size());


  // cout << "SwapImprove" << endl;
  // cout << "Vol = " << CalcVolume (points, elements) << endl;

  bad1 = 0;
  for (i = 1; i <= elements.Size(); i++)
    bad1 += CalcBad (points, elements.Get(i), 0);
  // cout << "Total badness = " << bad1 << endl;

  for (i = 1; i <= surfelements.Size(); i++)
    for (j = 1; j <= 3; j++)
    {
      i2.I1() = surfelements.Get(i).PNumMod(j);
      i2.I2() = surfelements.Get(i).PNumMod(j+1);
      if (i2.I1() > i2.I2()) swap (i2.I1(), i2.I2());
      boundaryedges.Set (i2, 1);
    }

  for (i = 1; i <= elements.Size(); i++)
    for (j = 1; j <= elements.Get(i).NP(); j++)
      elementsonnode.Add (elements.Get(i).PNum(j), i);

  for (i = 1; i <= elements.Size(); i++)
  {
    for (j = 1; j <= 6; j++)
    {
      elem = &elements.Elem(i);
      if (elem->PNum(1) == 0) continue;

      switch (j)
      {
      case 1 : pi1 = elem->PNum(1); pi2 = elem->PNum(2); break;
      case 2 : pi1 = elem->PNum(1); pi2 = elem->PNum(3); break;
      case 3 : pi1 = elem->PNum(1); pi2 = elem->PNum(4); break;
      case 4 : pi1 = elem->PNum(2); pi2 = elem->PNum(3); break;
      case 5 : pi1 = elem->PNum(2); pi2 = elem->PNum(4); break;
      case 6 : pi1 = elem->PNum(3); pi2 = elem->PNum(4); break;
      }

      if (pi2 < pi1) swap (pi1, pi2);

      i2.I1() = pi1;
      i2.I2() = pi2;
      if (boundaryedges.Used (i2)) continue;
      //      if (pi2 <= nboundnodes) continue;

      hasbothpoints.SetSize (0);
      for (k = 1; k <= elementsonnode.EntrySize(pi1); k++)
      {
        has1 = has2 = 0;
        elnr = elementsonnode.Get(pi1, k);
        elem = &elements.Elem (elnr);
        for (l = 1; l <= elem->NP(); l++)
        {
          if (elem->PNum(l) == pi1) has1 = 1;
          if (elem->PNum(l) == pi2) has2 = 1;
        }
        if (has1 && has2)
        {   // only once
          for (l = 1; l <= hasbothpoints.Size(); l++)
            if (hasbothpoints.Get(l) == elnr)
              has1 = 0;

          if (has1)
            hasbothpoints.Append (elnr);
        }
      }


      if (hasbothpoints.Size() == 3)
      {
        elem = &elements.Elem (hasbothpoints.Get(1));
        for (l = 1; l <= 4; l++)
          if (elem->PNum(l) != pi1 && elem->PNum(l) != pi2)
          {
            pi4 = pi3;
            pi3 = elem->PNum(l);
          }

        el31.PNum(1) = pi1;
        el31.PNum(2) = pi2;
        el31.PNum(3) = pi3;
        el31.PNum(4) = pi4;
        if (CalcBad (points, el31, 0) > 1e8)
        {
          swap (pi3, pi4);
          el31.PNum(3) = pi3;
          el31.PNum(4) = pi4;
        }

        pi5 = 0;
        for (k = 2; k <= 3; k++)
        {
          elem = &elements.Elem(hasbothpoints.Get(k));
          has1 = 0;
          for (l = 1; l <= 4; l++)
            if (elem->PNum(l) == pi4)
              has1 = 1;
          if (has1)
          {
            for (l = 1; l <=4; l++)
              if (elem->PNum(l) != pi1 && elem->PNum(l) != pi2 &&
                  elem->PNum(l) != pi4)
                pi5 = elem->PNum(l);
          }
        }

        el32.PNum(1) = pi1;
        el32.PNum(2) = pi2;
        el32.PNum(3) = pi4;
        el32.PNum(4) = pi5;

        el33.PNum(1) = pi1;
        el33.PNum(2) = pi2;
        el33.PNum(3) = pi5;
        el33.PNum(4) = pi3;

        elementsonnode.Add (pi4, hasbothpoints.Elem(2));
        elementsonnode.Add (pi3, hasbothpoints.Elem(3));

        bad1 = CalcBad (points, el31, 0) + CalcBad (points, el32, 0) +
               CalcBad (points, el33, 0);

        el21.PNum(1) = pi3;
        el21.PNum(2) = pi4;
        el21.PNum(3) = pi5;
        el21.PNum(4) = pi2;

        el22.PNum(1) = pi5;
        el22.PNum(2) = pi4;
        el22.PNum(3) = pi3;
        el22.PNum(4) = pi1;

        bad2 = CalcBad (points, el21, 0) + CalcBad (points, el22, 0);

        if (bad2 < bad1)
        {
          // cout << "3->2 conversion" << endl;
          elements.Elem (hasbothpoints.Get(1)) = el21;
          elements.Elem (hasbothpoints.Get(2)) = el22;
          for (l = 1; l <= 4; l++)
            elements.Elem (hasbothpoints.Get(3)).PNum(l) = 0;

          for (k = 1; k <= 2; k++)
            for (l = 1; l <= 4; l++)
              elementsonnode.Add (
                elements.Get(hasbothpoints.Get(k)).PNum(l),
                hasbothpoints.Get(k));
        }
      }



      if (hasbothpoints.Size() == 4)
      {
        elem = &elements.Elem (hasbothpoints.Get(1));
        for (l = 1; l <= 4; l++)
          if (elem->PNum(l) != pi1 && elem->PNum(l) != pi2)
          {
            pi4 = pi3;
            pi3 = elem->PNum(l);
          }

        el1.PNum(1) = pi1;
        el1.PNum(2) = pi2;
        el1.PNum(3) = pi3;
        el1.PNum(4) = pi4;
        if (CalcBad (points, el1, 0) > 1e8)
        {
          swap (pi3, pi4);
          el1.PNum(3) = pi3;
          el1.PNum(4) = pi4;
        }

        pi5 = 0;
        for (k = 2; k <= 4; k++)
        {
          elem = &elements.Elem(hasbothpoints.Get(k));
          has1 = 0;
          for (l = 1; l <= 4; l++)
            if (elem->PNum(l) == pi4)
              has1 = 1;
          if (has1)
          {
            for (l = 1; l <=4; l++)
              if (elem->PNum(l) != pi1 && elem->PNum(l) != pi2 &&
                  elem->PNum(l) != pi4)
                pi5 = elem->PNum(l);
          }
        }

        pi6 = 0;
        for (k = 2; k <= 4; k++)
        {
          elem = &elements.Elem(hasbothpoints.Get(k));
          has1 = 0;
          for (l = 1; l <= 4; l++)
            if (elem->PNum(l) == pi3)
              has1 = 1;
          if (has1)
          {
            for (l = 1; l <=4; l++)
              if (elem->PNum(l) != pi1 && elem->PNum(l) != pi2 &&
                  elem->PNum(l) != pi3)
                pi6 = elem->PNum(l);
          }
        }


        el1.PNum(1) = pi1;
        el1.PNum(2) = pi2;
        el1.PNum(3) = pi3;
        el1.PNum(4) = pi4;

        el2.PNum(1) = pi1;
        el2.PNum(2) = pi2;
        el2.PNum(3) = pi4;
        el2.PNum(4) = pi5;

        el3.PNum(1) = pi1;
        el3.PNum(2) = pi2;
        el3.PNum(3) = pi5;
        el3.PNum(4) = pi6;

        el4.PNum(1) = pi1;
        el4.PNum(2) = pi2;
        el4.PNum(3) = pi6;
        el4.PNum(4) = pi3;

        //        elementsonnode.Add (pi4, hasbothpoints.Elem(2));
        //        elementsonnode.Add (pi3, hasbothpoints.Elem(3));

        bad1 = CalcBad (points, el1, 0) + CalcBad (points, el2, 0) +
               CalcBad (points, el3, 0) + CalcBad (points, el4, 0);


        el1.PNum(1) = pi3;
        el1.PNum(2) = pi5;
        el1.PNum(3) = pi2;
        el1.PNum(4) = pi4;

        el2.PNum(1) = pi3;
        el2.PNum(2) = pi5;
        el2.PNum(3) = pi4;
        el2.PNum(4) = pi1;

        el3.PNum(1) = pi3;
        el3.PNum(2) = pi5;
        el3.PNum(3) = pi1;
        el3.PNum(4) = pi6;

        el4.PNum(1) = pi3;
        el4.PNum(2) = pi5;
        el4.PNum(3) = pi6;
        el4.PNum(4) = pi2;

        bad2 = CalcBad (points, el1, 0) + CalcBad (points, el2, 0) +
               CalcBad (points, el3, 0) + CalcBad (points, el4, 0);

        el1b.PNum(1) = pi4;
        el1b.PNum(2) = pi6;
        el1b.PNum(3) = pi3;
        el1b.PNum(4) = pi2;

        el2b.PNum(1) = pi4;
        el2b.PNum(2) = pi6;
        el2b.PNum(3) = pi2;
        el2b.PNum(4) = pi5;

        el3b.PNum(1) = pi4;
        el3b.PNum(2) = pi6;
        el3b.PNum(3) = pi5;
        el3b.PNum(4) = pi1;

        el4b.PNum(1) = pi4;
        el4b.PNum(2) = pi6;
        el4b.PNum(3) = pi1;
        el4b.PNum(4) = pi3;

        bad3 = CalcBad (points, el1b, 0) + CalcBad (points, el2b, 0) +
               CalcBad (points, el3b, 0) + CalcBad (points, el4b, 0);


        if (bad2 < bad1 && bad2 < bad3)
        {
          // cout << "bad1 = " << bad1 << " bad2 = " << bad2 << endl;

          elements.Elem (hasbothpoints.Get(1)) = el1;
          elements.Elem (hasbothpoints.Get(2)) = el2;
          elements.Elem (hasbothpoints.Get(3)) = el3;
          elements.Elem (hasbothpoints.Get(4)) = el4;

          for (k = 1; k <= 4; k++)
            for (l = 1; l <= 4; l++)
              elementsonnode.Add (
                elements.Get(hasbothpoints.Get(k)).PNum(l),
                hasbothpoints.Get(k));
        }
        else if (bad3 < bad1)
        {
          // cout << "bad1 = " << bad1 << " bad3 = " << bad3 << endl;

          elements.Elem (hasbothpoints.Get(1)) = el1b;
          elements.Elem (hasbothpoints.Get(2)) = el2b;
          elements.Elem (hasbothpoints.Get(3)) = el3b;
          elements.Elem (hasbothpoints.Get(4)) = el4b;

          for (k = 1; k <= 4; k++)
            for (l = 1; l <= 4; l++)
              elementsonnode.Add (
                elements.Get(hasbothpoints.Get(k)).PNum(l),
                hasbothpoints.Get(k));
        }
      }




    }
  }

  for (k = elements.Size(); k >= 1; k--)
    if (elements[k].PNum(1) == 0)
      elements.DeleteElement(k);

  // cout << "Vol = " << CalcVolume (points, elements) << endl;
}
