// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <new.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include <myadt.hh>

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


// global variables for plotting:

static ARRAY<Point3d> locpoints;
static ARRAY<Element> locfaces;
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


static Meshing3 * meshing;



// for useful point
static ARRAY<Point3d> grouppoints;
static ARRAY<Element> groupfaces;
static ARRAY<INDEX> grouppindex, groupfindex;

extern int FindInnerPoint (const ARRAY<Point3d> & grouppoints,
                           const ARRAY<Element> & groupfaces,
                           Point3d & p);


Meshing3 :: Meshing3 (char * rulefilename)
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

Meshing3 :: ~Meshing3 ()
{}

void Meshing3 :: AddPoint (const Point3d & p, INDEX globind)
{
  adfront->AddPoint (p, globind);
}

void Meshing3 :: AddBoundaryElement (const Element & elem, int inverse)
{
  Element helem;
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






void DrawPnF(const ARRAY<Point3d> & pointl, const ARRAY<Element> & facel,
             ARRAY<int> & style, float xh , const ROT3D & r)
{
#ifdef MYGRAPH
  INDEX i;
  int j;
  const Element * face;
  int points2d[8];
  int nold = 0, nnew = 0;

  Point3d c, c2;

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



void Meshing3 :: Mesh (double ah)
{
  ARRAY<Point3d> plainpoints;
  ARRAY<int> delpoints, delfaces;
  ARRAY<int> changepoints;
  ARRAY<Element> locelements;
  int i, j, oldnp, oldnf;
  int found;
  referencetransform trans;
  int rotind;
  INDEX globind;
  Point3d inp;

  float minerr;
  int hasfound;
  double tetvol;

  ARRAY<Point3d> tempnewpoints;
  ARRAY<Element> tempnewfaces;
  ARRAY<int> tempdelfaces;
  ARRAY<Element> templocelements;

  int pause = 1, redraw = 1, shouldredraw = 1;


  h = ah;
  meshing = this;
  adfront->SetStartFront ();


  found = 0;
  vol0 = adfront -> Volume();
  tetvol = 0;



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
        Element newel;

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

        ARRAY<Point3d> pts;
        Point3d hp;

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
        Point3d * hp1, * hp2, * hp3, * hp4;
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
