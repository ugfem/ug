// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   meshing3.cc                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <new.h>
#include <ctype.h>
#include <string.h>

#include <template.hh>
#include <array.hh>
#include <table.hh>

#include <geom/geom2d.hh>
#include <geom/geom3d.hh>
#include <geom/reftrans.hh>

#ifdef MYGRAPH
#include <geom/rot3d.hh>
#include <graphics/mygraph.hh>
#endif

#include <linalg/linalg.hh>

#include <meshing/global.hh>
#include <meshing/meshing3.hh>

#ifdef SOLIDGEOM
#include <geom/solid.hh>
#include <meshing/specpoin.hh>
#include <meshing/edgeflw.hh>
#endif

extern void BFGS (VECTOR & x, double (*f)(const VECTOR & x, VECTOR & g));

#ifdef SOLIDGEOM
extern ARRAY<SURFACE*> surfaces;
#endif

// global variables for plotting:

static ARRAY<POINT3D> locpoints;
static ARRAY<ELEMENT> locfaces;
static ARRAY<INDEX> pindex, findex;
static ARRAY<int> style;
static int drawrad = 100;
static int cntsucc = 0;
static long cnttrials = 0;
static int cntelem = 0;
static char buf[100];
static int qualclass;
static float vol0, err, h;
static int problemindex = 1;


static meshing3 * meshing;



// for useful point
static ARRAY<POINT3D> grouppoints;
static ARRAY<ELEMENT> groupfaces;
static ARRAY<INDEX> grouppindex, groupfindex;

extern int FindInnerPoint (const ARRAY<POINT3D> & grouppoints,
                           const ARRAY<ELEMENT> & groupfaces,
                           POINT3D & p);


meshing3 :: meshing3 (char * rulefilename)
{
  int i;

  LoadRules (rulefilename);
  adfront = new ADFRONT3;

  ruleused.SetSize (rules.Size());
  for (i = 1; i <= ruleused.Size(); i++)
    ruleused[i] = 0;

  problems.SetSize (rules.Size());
  for (i = 1; i <= problems.Size(); i++)
    problems[i] = new char[255];
}

meshing3 :: ~meshing3 ()
{}

void meshing3 :: AddPoint (const POINT3D & p, INDEX globind)
{
  adfront->AddPoint (p, globind);
}

void meshing3 :: AddBoundaryElement (const ELEMENT & elem, int inverse)
{
  ELEMENT helem;
  INDEX hi;

  helem = elem;
  if (inverse)
  {
    hi = helem.PNum(1);
    helem.PNum(1) = helem.PNum(2);
    helem.PNum(2) = hi;
  }

  adfront->AddFace (helem);
}






void DrawPnF(const ARRAY<POINT3D> & pointl, const ARRAY<ELEMENT> & facel,
             ARRAY<int> & style, float xh , const ROT3D & r)
{
#ifdef MYGRAPH
  INDEX i;
  int j;
  const ELEMENT * face;
  int points2d[8];
  int nold = 0, nnew = 0;

  POINT3D c, c2;

  MySetColor (LIGHTBLUE);

  c = Center (pointl[facel[1].PNum(1)], pointl[facel[1].PNum(2)],
              pointl[facel[1].PNum(3)]);

  for (i = 1; i <= facel.Size(); i++)
  {
    c2 = Center (pointl[facel[i].PNum(1)], pointl[facel[i].PNum(2)],
                 pointl[facel[i].PNum(3)]);

    if (Dist(c, c2) > xh) continue;

    face = &facel[i];
    for (j = 0; j <= face->NP(); j++)
    {
      points2d[2*j  ] = int (r.X(pointl[face->PNumMod(j+1)]));
      points2d[2*j+1] = int (r.Y(pointl[face->PNumMod(j+1)]));
    }

    switch (style[i])
    {
    case 1 : case 2 :
    {
      MySetColor (RED);
      MyFillPoly (4, points2d, nold % 4);
      nold++;
      break;
    }
    case 3 :
    {
      MySetColor (BLUE);
      MyFillPoly (4, points2d, nnew % 4 + 4);
      nnew++;
      break;
    }
    default :
    {
      MySetColor (LIGHTBLUE);
      MyLine3D (pointl[facel[i].PNum(1)], pointl[facel[i].PNum(2)], r);
      MyLine3D (pointl[facel[i].PNum(1)], pointl[facel[i].PNum(3)], r);
      MyLine3D (pointl[facel[i].PNum(2)], pointl[facel[i].PNum(3)], r);
    }
    }
  }
#endif
}




void PlotVolMesh (const ROT3D & r, char key)
{
#ifdef MYGRAPH
  int i, y;
  char cbuf[20];

  MyClearDevice();

  MySetColor (BLUE);
  y = 30;

  //  sprintf (buf, "Trials:   %ld", cnttrials);
  //  MyOutTextXY (20, y, buf);
  //  y += 15;
  //  sprintf (buf, "Success:  %d", cntsucc);
  //  MyOutTextXY (20, y, buf);
  //  y += 15;
  sprintf (buf, "Elements: %d", cntelem);
  MyOutTextXY (20, y, buf);
  y += 15;
  sprintf (buf, "Quality:  %d", qualclass);
  MyOutTextXY (20, y, buf);
  y += 15;
  sprintf (buf, "Volume:   %4.1f %%", float(100 * meshing->adfront->Volume() / vol0));
  MyOutTextXY (20, y, buf);
  y += 15;
  //  sprintf (buf, "Error:    %f", err);
  //  MyOutTextXY (20, y, buf);
  //  y += 15;
  sprintf (buf, "Time:     %d", GetTime());
  MyOutTextXY (20, y, buf);
  y += 15;

  DrawPnF(locpoints, locfaces, style, drawrad, r);


  for (i = 1; i <= locpoints.Size(); i++)
  {
    MySetColor(BLUE);
    MyPoint3D (locpoints[i], r);
    if (testmode && i <= pindex.Size())
    {
      MySetColor(LIGHTBLUE);
      sprintf (cbuf, "%d", pindex[i]);
      MyOutTextXY (r.X(locpoints[i]), r.Y(locpoints[i]), cbuf);
    }
  }

  if (testmode)
  {
    MySetColor (BLACK);
    sprintf (buf, "Problem with rule %d: %s",
             problemindex, meshing->rules[problemindex]->Name());
    MyOutTextXY (120, 30, buf);
    MyOutTextXY (120, 45, meshing->problems[problemindex]);
  }

  if (key == 's')
    meshing -> adfront -> SaveSurface ("surf.out", h);
#endif
}



void meshing3 :: Mesh (double ah)
{
  ARRAY<POINT3D> plainpoints;
  ARRAY<int> delpoints, delfaces;
  ARRAY<int> changepoints;
  ARRAY<ELEMENT> locelements;
  int i, j, oldnp, oldnf;
  int found;
  referencetransform trans;
  int rotind;
  INDEX globind;
  POINT3D inp;

  float minerr;
  int hasfound;
  double tetvol;

  ARRAY<POINT3D> tempnewpoints;
  ARRAY<ELEMENT> tempnewfaces;
  ARRAY<int> tempdelfaces;
  ARRAY<ELEMENT> templocelements;

  int pause = 1, redraw = 1, shouldredraw = 1;


  h = ah;
  meshing = this;
  adfront->SetStartFront ();


  found = 0;
  vol0 = adfront -> Volume();
  tetvol = 0;

  //  adfront -> GetGroup (1, grouppoints, groupfaces);
  //  FindInnerPoint (grouppoints, groupfaces, inp);

  //  grouppoints.Append (inp);
  //  adfront->AddPoint (inp, 10000);

  //  Plot3D (PlotGroup, 1, 1);


  while (!adfront -> Empty())

  {
    locpoints.SetSize(0);
    locfaces.SetSize(0);
    pindex.SetSize(0);
    findex.SetSize(0);

    qualclass =
      adfront -> GetLocals (locpoints, locfaces, pindex, findex, 3 * h);

    oldnp = locpoints.Size();
    oldnf = locfaces.Size();

    if (qualclass >= 4)
    {
      grouppoints.SetSize (0);
      groupfaces.SetSize (0);
      grouppindex.SetSize (0);
      groupfindex.SetSize (0);

      adfront -> GetGroup (findex[1], grouppoints, groupfaces, grouppindex, groupfindex);
      if (groupfaces.Size() <= 20 && FindInnerPoint (grouppoints, groupfaces, inp))
      {
        (*testout) << "inner point found" << endl;

        for (i = 1; i <= groupfaces.Size(); i++)
          adfront -> DeleteFace (groupfindex[i]);

        for (i = 1; i <= groupfaces.Size(); i++)
          for (j = 1; j <= locfaces.Size(); j++)
            if (findex[j] == groupfindex[i])
              delfaces.Append (j);

        style.SetSize (locfaces.Size());
        for (i = 1; i <= style.Size(); i++)
          style[i] = 0;
        for (i = 1; i <= delfaces.Size(); i++)
          style[delfaces[i]] = 2;
        style[1] = 1;

#ifdef MYGRAPH
        Plot3D (PlotVolMesh, 1, pause);
#endif
        pause = 0;

        delfaces.SetSize (0);

        INDEX npi;
        ELEMENT newel;

        npi = SavePoint (inp);
        newel.SetNP(4);
        newel.PNum(4) = npi;

        for (i = 1; i <= groupfaces.Size(); i++)
        {
          for (j = 1; j <= 3; j++)
            newel.PNum(j) = adfront->GetGlobalIndex (groupfaces.Get(i).PNum(j));
          SaveElement (newel);
        }
        continue;
      }
    }



    found = 0;
    hasfound = 0;
    minerr = 1000;

    for (rotind = 1; rotind <= 3; rotind++)
    {
      trans.Set (locpoints[locfaces[1].PNumMod(1+rotind)],
                 locpoints[locfaces[1].PNumMod(2+rotind)],
                 locpoints[locfaces[1].PNumMod(3+rotind)], h);

      trans.ToPlain (locpoints, plainpoints);

      cnttrials++;

      found = ApplyVRules (rules, tolfak, plainpoints, locfaces, locelements,
                           delfaces, qualclass, rotind, err, problems);


      if (found) cntsucc++;


      locpoints.SetSize (plainpoints.Size());
      for (i = oldnp+1; i <= plainpoints.Size(); i++)
        trans.FromPlain (plainpoints[i], locpoints[i]);

      style.SetSize (locfaces.Size());
      for (i = 1; i <= oldnf; i++)
        style[i] = 0;
      for (i = oldnf+1; i <= locfaces.Size(); i++)
        style[i] = 3;
      for (i = 1; i <= delfaces.Size(); i++)
        style[delfaces[i]] = 2;
      style[1] = 1;

      if (found)
      {
        ruleused[found]++;

        ARRAY<POINT3D> pts;
        POINT3D hp;

        pts.SetSize (rules[found]->GetTransFreeZone().Size());

        for (i = 1; i <= pts.Size(); i++)
        {
          trans.FromPlain (rules[found]->GetTransFreeZone().Get(i), hp);
          pts.Set(i, hp);
        }
      }

      if (testmode && found)
      {
        pause = 1;
        redraw = 1;
      }
      if (found && shouldredraw)
        redraw = 1;

      if (qualclass > 30)
        pause = 1;

      if (found && rules[found]->TestFlag ('p'))
      {
        redraw = 1;
        pause = 1;
      }

#ifdef MYGRAPH
      Plot3D (PlotVolMesh, redraw, pause);
#endif
      redraw = 0;
      pause = 0;


      if (found && (!hasfound || err < minerr) )
      {

        if (testmode)
        {
          (*testout) << "testmode found" << endl;
          for (i = 1; i <= plainpoints.Size(); i++)
          {
            (*testout) << "p";
            if (i <= pindex.Size())
              (*testout) << pindex[i] << ": ";
            else
              (*testout) << "new: ";
            (*testout) << plainpoints[i] << endl;
          }
        }



        hasfound = found;
        minerr = err;

        tempnewpoints.SetSize (0);
        for (i = oldnp+1; i <= locpoints.Size(); i++)
          tempnewpoints.Append (locpoints[i]);

        tempnewfaces.SetSize (0);
        for (i = oldnf+1; i <= locfaces.Size(); i++)
          tempnewfaces.Append (locfaces[i]);

        tempdelfaces.SetSize (0);
        for (i = 1; i <= delfaces.Size(); i++)
          tempdelfaces.Append (delfaces[i]);

        templocelements.SetSize (0);
        for (i = 1; i <= locelements.Size(); i++)
          templocelements.Append (locelements[i]);
      }

      locpoints.SetSize (oldnp);
      locfaces.SetSize (oldnf);
      delfaces.SetSize (0);
      locelements.SetSize (0);
    }


    if (hasfound)
    {
      for (i = 1; i <= tempnewpoints.Size(); i++)
        locpoints.Append (tempnewpoints[i]);
      for (i = 1; i <= tempnewfaces.Size(); i++)
        locfaces.Append (tempnewfaces[i]);
      for (i = 1; i <= tempdelfaces.Size(); i++)
        delfaces.Append (tempdelfaces[i]);
      for (i = 1; i <= templocelements.Size(); i++)
        locelements.Append (templocelements[i]);


      /*
            changepoints.SetSize (locpoints.Size());
            for (i = 1; i <= changepoints.Size(); i++)
              changepoints[i] = 0;

            for (i = oldnf + 1; i <= locfaces.Size(); i++)
              for (j = 1; j <= locfaces[i].NP(); j++)
                changepoints[locfaces[i].PNum(j)] = 1;

            for (i = 1; i <= oldnf; i++)
              {
              hc = 0;
              for (j = 1; j <= locfaces[i].NP(); j++)
                if (changepoints[locfaces[i].PNum(j)])
                  hc = 1;

              if (hc)
                adfront->ResetClass (findex[i]);
              }
       */

      if (testmode)
      {
        (*testout) << "testmode locpoints" << endl;
        for (i = 1; i <= locpoints.Size(); i++)
        {
          (*testout) << "p";
          if (i <= pindex.Size())
            (*testout) << pindex[i] << ": ";
          else
            (*testout) << "new: ";
          (*testout) << locpoints[i] << endl;
        }
      }



      pindex.SetSize(locpoints.Size());

      for (i = oldnp+1; i <= locpoints.Size(); i++)
      {
        globind = SavePoint (locpoints[i]);
        pindex[i] = adfront -> AddPoint (locpoints[i], globind);
      }

      for (i = 1; i <= locelements.Size(); i++)
      {
        POINT3D * hp1, * hp2, * hp3, * hp4;
        hp1 = &locpoints[locelements[i].PNum(1)];
        hp2 = &locpoints[locelements[i].PNum(2)];
        hp3 = &locpoints[locelements[i].PNum(3)];
        hp4 = &locpoints[locelements[i].PNum(4)];

        tetvol += (1.0 / 6.0) * ( Cross ( *hp2 - *hp1, *hp3 - *hp1) * (*hp4 - *hp1) );

        for (j = 1; j <= locelements[i].NP(); j++)
          locelements[i].PNum(j) =
            adfront -> GetGlobalIndex (pindex[locelements[i].PNum(j)]);

        SaveElement (locelements[i]);
        cntelem++;
      }

      for (i = oldnf+1; i <= locfaces.Size(); i++)
      {
        for (j = 1; j <= locfaces[i].NP(); j++)
          locfaces[i].PNum(j) = pindex[locfaces[i].PNum(j)];

        adfront->AddFace (locfaces[i]);
      }

      for (i = 1; i <= delfaces.Size(); i++)
      {
        adfront->DeleteFace (findex[delfaces[i]]);
      }

      //      testout << "rule: " << rules[hasfound] -> Name() << endl;
      //      testout << "Vol from surface: " << (vol0 - adfront -> Volume()) << endl;
      //      testout << "Tetra volume:     " << tetvol << endl << endl;
    }
    else
    {
      adfront->IncrementClass (findex[1]);
    }

    locelements.SetSize (0);
    delpoints.SetSize(0);
    delfaces.SetSize(0);
  }



  for (i = 1; i <= ruleused.Size(); i++)
    (*testout) << setw(4) << ruleused[i]
               << " times used rule " << rules[i] -> Name() << endl;
}





static POINT3D sp1;
static VEC3D n, t1, t2;
static ARRAY<INDEX> locelements(0);
static ARRAY<int> locrots(0);
static ARRAY<POINT3D> * optpoints;
static const ARRAY<ELEMENT> * optelements;
static int locerr2;
static double loch;
static int surfi, surfi2;



static double CalcTetBadness (const POINT3D & p1, const POINT3D & p2,
                              const POINT3D & p3, const POINT3D & p4)
{
  double vol, l, l4, l5, l6;
  double err1, err2;

  VEC3D v1 = p2 - p1;
  VEC3D v2 = p3 - p1;
  VEC3D v3 = p4 - p1;

  vol = - (Cross (v1, v2) * v3) / 6;
  l4 = Dist (p2, p3);
  l5 = Dist (p2, p4);
  l6 = Dist (p3, p4);

  l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

  if (vol < 1e-8) return 1e10;
  err1 = (l*l*l) / (1832.82 * vol);    // 6^4 * sqrt(2)
  err2 = l / (6 * loch) + loch * ( 1 / v1.Length() + 1 / v2.Length() + 1 / v3.Length() + 1 / l4 + 1 / l5 + 1 / l6);
  return err1 + err2;
}


double PointFunctionValue (const POINT3D & pp)
{
  int j, k;
  INDEX eli;
  const ELEMENT * el;
  double badness;
  ARRAY<const POINT3D*> p(4);

  badness = 0;

  for (j = 1; j <= locelements.Size(); j++)
  {
    eli = locelements.Get(j);
    el = &optelements->Get(eli);

    for (k = 1; k <= 4; k++)
      p.Elem(k) = &optpoints->Get(el->PNum(k));
    p.Elem(locrots.Get(j)) = &pp;

    badness += CalcTetBadness (*p.Get(1), *p.Get(2), *p.Get(3), *p.Get(4));
  }
  return badness;
}


double PointFunctionValueGrad (const POINT3D & pp, VECTOR & grad)
{
  double f, fr, fl, delta = 1e-6;
  POINT3D hpp;

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

double Opti3FunctionValueGrad (const VECTOR & x, VECTOR & grad)
{
  POINT3D pp;
  pp.X() = sp1.X() + x.Get(1);
  pp.Y() = sp1.Y() + x.Get(2);
  pp.Z() = sp1.Z() + x.Get(3);

  return PointFunctionValueGrad (pp, grad);
}


#ifdef SOLIDGEOM
double Opti3SurfFunctionValueGrad (const VECTOR & x, VECTOR & grad)
{
  VEC3D n, v1, v2, e1, e2, vgrad;
  POINT3D pp1;
  double badness;
  static VECTOR freegrad(3);

  pp1 = sp1 + (x.Get(1) * t1) + (x.Get(2) * t2);
  surfaces[surfi] -> Project (pp1);

  badness = PointFunctionValueGrad (pp1, freegrad);
  vgrad.X() = freegrad.Get(1);
  vgrad.Y() = freegrad.Get(2);
  vgrad.Z() = freegrad.Get(3);

  surfaces[surfi] -> GetNormalVector (pp1, n);

  vgrad -= (vgrad * n) * n;

  grad.Elem(1) = vgrad * t1;
  grad.Elem(2) = vgrad * t2;
  return badness;
}
#endif


#ifdef SOLIDGEOM
double Opti3EdgeFunctionValueGrad (const VECTOR & x, VECTOR & grad)
{
  VEC3D n1, n2, v1, vgrad;
  POINT3D pp1;
  double badness;
  static VECTOR freegrad(3);

  pp1 = sp1 + x.Get(1) * t1;
  ProjectToEdge ( surfaces[surfi], surfaces[surfi2], pp1);

  badness = PointFunctionValueGrad (pp1, freegrad);

  vgrad.X() = freegrad.Get(1);
  vgrad.Y() = freegrad.Get(2);
  vgrad.Z() = freegrad.Get(3);

  surfaces[surfi] -> GetNormalVector (pp1, n1);
  surfaces[surfi2] -> GetNormalVector (pp1, n2);

  v1 = Cross (n1, n2);
  v1 /= v1.Length();

  grad.Elem(1) = (vgrad * v1) * (t1 * v1);
  return badness;
}
#endif




#ifdef SOLIDGEOM
void meshing3 :: ImproveMesh (ARRAY<POINT3D> & points,
                              const ARRAY<ELEMENT> & surfelements,
                              const ARRAY<ELEMENT> & elements, double h)
{
  INDEX i, eli;
  int j, k, it, ito, rot, surfi3;
  const ELEMENT * el;

  VEC3D n1, n2;
  TABLE<INDEX> elementsonpoint(points.Size());
  TABLE<INDEX> surfelementsonpoint(points.Size());
  VECTOR x(3), xsurf(2), xedge(1);

  (*testout).precision(8);

  for (i = 1; i <= elements.Size(); i++)
    for (j = 1; j <= elements[i].NP(); j++)
      elementsonpoint.Add (elements[i].PNum(j), i);
  for (i = 1; i <= surfelements.Size(); i++)
    for (j = 1; j <= surfelements[i].NP(); j++)
      surfelementsonpoint.Add (surfelements[i].PNum(j), i);

  loch = h;
  optpoints = &points;
  optelements = &elements;

  for (i = 1; i <= points.Size(); i++)
  {
    //    testout << "i = " << i << endl;
    sp1 = points.Elem(i);

    locelements.SetSize(0);
    locrots.SetSize (0);

    surfi = surfi2 = surfi3 = 0;

    for (j = 1; j <= surfelementsonpoint.EntrySize(i); j++)
    {
      eli = surfelementsonpoint.Get(i, j);
      el = &surfelements.Get(eli);

      if (!surfi)
        surfi = el->SurfaceIndex();
      else if (surfi != el->SurfaceIndex())
      {
        if (!surfi2)
          surfi2 = el->SurfaceIndex();
        else if (surfi2 != el->SurfaceIndex())
          surfi3 = el->SurfaceIndex();
      }
    }

    for (j = 1; j <= elementsonpoint.EntrySize(i); j++)
    {
      eli = elementsonpoint.Get(i, j);
      el = &elements.Get(eli);

      locelements.Append (eli);

      for (k = 1; k <= el->NP(); k++)
        if (el->PNum(k) == i)
          locrots.Append (k);
    }

    if (surfi2 && !surfi3)
    {
      (*testout) << "Edgepoint" << endl;
      surfaces[surfi] -> GetNormalVector (sp1, n1);
      surfaces[surfi2] -> GetNormalVector (sp1, n2);
      t1 = Cross (n1, n2);

      xedge = 0;
      BFGS (xedge, Opti3EdgeFunctionValueGrad);

      points.Elem(i).X() += xedge.Get(1) * t1.X();
      points.Elem(i).Y() += xedge.Get(1) * t1.Y();
      points.Elem(i).Z() += xedge.Get(1) * t1.Z();
      ProjectToEdge (surfaces[surfi], surfaces[surfi2], points.Elem(i));
    }

    if (surfi && !surfi2)
    {
      //      testout << "Surfacepoint" << endl;
      surfaces[surfi] -> GetNormalVector (sp1, n);
      if (fabs (n.X()) > 0.3)
      {
        t1.X() = n.Y();
        t1.Y() = -n.X();
        t1.Z() = 0;
        t1 /= t1.Length();
      }
      else
      {
        t1.X() = 0;
        t1.Y() = -n.Z();
        t1.Z() = n.Y();
        t1 /= t1.Length();
      }
      t2 = Cross (n, t1);

      xsurf = 0;
      BFGS (xsurf, Opti3SurfFunctionValueGrad);

      points.Elem(i).X() += (xsurf.Get(1) * t1.X() + xsurf.Get(2) * t2.X());
      points.Elem(i).Y() += (xsurf.Get(1) * t1.Y() + xsurf.Get(2) * t2.Y());
      points.Elem(i).Z() += (xsurf.Get(1) * t1.Z() + xsurf.Get(2) * t2.Z());
      surfaces[surfi]->Project (points.Elem(i));
    }

    if (!surfi)
    {
      //      testout << "Volumepoint" << endl;
      x = 0;
      BFGS (x, Opti3FunctionValueGrad);

      points.Elem(i).X() += x.Get(1);
      points.Elem(i).Y() += x.Get(2);
      points.Elem(i).Z() += x.Get(3);
    }
  }
}
#endif


void meshing3 :: ImproveMesh (ARRAY<POINT3D> & points,
                              const ARRAY<ELEMENT> & elements,
                              int nboundnodes, double h)
{
  INDEX i, eli;
  int j, k, it, ito, rot, surfi3;
  const ELEMENT * el;

  (*testout) << "Improve Mesh" << endl;

  VEC3D n1, n2;
  TABLE<INDEX> elementsonpoint(points.Size());
  VECTOR x(3);

  (*testout).precision(8);

  for (i = 1; i <= elements.Size(); i++)
    for (j = 1; j <= elements[i].NP(); j++)
      elementsonpoint.Add (elements[i].PNum(j), i);

  loch = h;
  optpoints = &points;
  optelements = &elements;



  for (i = nboundnodes+1; i <= points.Size(); i++)
  {
    (*testout) << "Node " << i << endl;
    sp1 = points.Elem(i);

    locelements.SetSize(0);
    locrots.SetSize (0);

    for (j = 1; j <= elementsonpoint.EntrySize(i); j++)
    {
      eli = elementsonpoint.Get(i, j);
      el = &elements.Get(eli);

      locelements.Append (eli);

      for (k = 1; k <= el->NP(); k++)
        if (el->PNum(k) == i)
          locrots.Append (k);
    }

    x = 0;
    (*testout) << "Start BFGS" << endl;

    BFGS (x, Opti3FunctionValueGrad);

    points.Elem(i).X() += x.Get(1);
    points.Elem(i).Y() += x.Get(2);
    points.Elem(i).Z() += x.Get(3);
  }
}
