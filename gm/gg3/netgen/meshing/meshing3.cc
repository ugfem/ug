// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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

/*#include "../../../include/gginterface.h"*/

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
static float vol0, err;
static int problemindex = 1;
static int MESH_DEBUG = 0;
double h;

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

void Meshing3 :: AddBoundaryElement (const Element & elem, int inverse, int prism_flag)
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

  adfront->AddFace (helem,prism_flag);
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
  //  if (GSF_DEBUG) MyOutTextXY (20, y, buf);
  //  y += 15;
  //  sprintf (buf, "Success:  %d", cntsucc);
  //  if (GSF_DEBUG) MyOutTextXY (20, y, buf);
  //  y += 15;
  sprintf (buf, "Elements: %d", cntelem);
  if (GSF_DEBUG) MyOutTextXY (20, y, buf);
  y += 15;
  sprintf (buf, "Quality:  %d", qualclass);
  if (GSF_DEBUG) MyOutTextXY (20, y, buf);
  y += 15;
  sprintf (buf, "Volume:   %4.1f %%", float(100 * meshing->adfront->Volume() / vol0));
  if (GSF_DEBUG) MyOutTextXY (20, y, buf);
  y += 15;
  //  sprintf (buf, "Error:    %f", err);
  //  if (GSF_DEBUG) MyOutTextXY (20, y, buf);
  //  y += 15;
  sprintf (buf, "Time:     %d", GetTime());
  if (GSF_DEBUG) MyOutTextXY (20, y, buf);
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
      if (GSF_DEBUG) MyOutTextXY (r.X(locpoints[i]), r.Y(locpoints[i]), cbuf);
    }
  }

  if (testmode)
  {
    MySetColor (BLACK);
    sprintf (buf, "Problem with rule %d: %s",
             problemindex, meshing->rules[problemindex]->Name());
    if (GSF_DEBUG) MyOutTextXY (120, 30, buf);
    if (GSF_DEBUG) MyOutTextXY (120, 45, meshing->problems[problemindex]);
  }

  if (key == 's')
    meshing -> adfront -> SaveSurface ("surf.out", h);
#endif
}


void Meshing3 :: Mesh (double ah,int make_prism)
{
  ARRAY<Point3d> plainpoints;
  ARRAY<int> delpoints, delfaces;
  ARRAY<int> changepoints;
  ARRAY<Element> locelements;
  int i, j, oldnp, oldnf;
  int found,rotind;
  referencetransform trans;
  int gen_prism, gen_pyramid, gen_tetrahedra, prism, pyramid,flag;
  INDEX globind;
  Point3d inp,bemp,bemp1,bemp2,bemp3;
  Point3d *hp1, *hp2, *hp3, *hp4, *hp5, *hp6;
  double in[5];
  float minerr;
  int hasfound;
  double tetvol,prismvol,pyramidvol;

  ARRAY<Point3d> tempnewpoints;
  ARRAY<Element> tempnewfaces;
  ARRAY<int> tempdelfaces;
  ARRAY<Element> templocelements;
  ARRAY<int> prism_flags;

  int pause = 1, redraw = 1, shouldredraw = 1;

  //  ah = 15;

  h = ah;
  meshing = this;
  adfront->SetStartFront ();

  found = 0;
  vol0 = adfront -> Volume();
  tetvol = 0.0;
  prismvol = 0.0;
  pyramidvol = 0.0;

  gen_pyramid = 0;
  qualclass = 1;
  //	gen_prism = adfront->Prism();
  gen_prism = make_prism;
  if(gen_prism)
    gen_tetrahedra = 0;
  else
    gen_tetrahedra = 1;

  // generate prism if possible
  if(gen_prism)
    while (!adfront -> Empty())
    {
      //adfront->Print();
      h = ah;
      locpoints.SetSize(0);
      locfaces.SetSize(0);
      pindex.SetSize(0);
      findex.SetSize(0);

      prism = adfront->GetLocals_Prism(locpoints, locfaces, pindex, findex, 3 * h, prism_flags);
      if(prism==-1)
      {
        gen_prism = 0;
        gen_pyramid = 1;
        break;
      }

      oldnp = locpoints.Size();
      oldnf = locfaces.Size();
      found = 0;
      hasfound = 0;
      minerr = 1000;

      trans.Set (     locpoints[locfaces[1].PNumMod(1)],
                      locpoints[locfaces[1].PNumMod(2)],
                      locpoints[locfaces[1].PNumMod(3)], h);

      trans.ToPlain (locpoints, plainpoints);

      cnttrials++;

      found = Generate_Prism (plainpoints, locfaces, locelements,
                              delfaces, prism_flags);

      if (found) cntsucc++;

      locpoints.SetSize (plainpoints.Size());
      for (i = oldnp+1; i <= plainpoints.Size(); i++)
        trans.FromPlain (plainpoints[i], locpoints[i]);

      tempnewfaces.SetSize (0);
      for (i = oldnf+1; i <= locfaces.Size(); i++)
        tempnewfaces.Append (locfaces[i]);

      tempdelfaces.SetSize (0);
      for (i = 1; i <= delfaces.Size(); i++)
        tempdelfaces.Append (delfaces[i]);

      templocelements.SetSize (0);
      for (i = 1; i <= locelements.Size(); i++)
        templocelements.Append (locelements[i]);

      locpoints.SetSize (oldnp);
      locfaces.SetSize (oldnf);
      delfaces.SetSize (0);
      locelements.SetSize (0);

      if (found)
      {
        for (i = 1; i <= tempnewfaces.Size(); i++)
          locfaces.Append (tempnewfaces[i]);
        for (i = 1; i <= tempdelfaces.Size(); i++)
          delfaces.Append (tempdelfaces[i]);
        for (i = 1; i <= templocelements.Size(); i++)
          locelements.Append (templocelements[i]);

        pindex.SetSize(locpoints.Size());

        for (i = 1; i <= locelements.Size(); i++)
        {
          hp1 = &locpoints[locelements[i].PNum(1)];
          hp2 = &locpoints[locelements[i].PNum(2)];
          hp3 = &locpoints[locelements[i].PNum(3)];
          hp4 = &locpoints[locelements[i].PNum(4)];
          hp5 = &locpoints[locelements[i].PNum(5)];
          hp6 = &locpoints[locelements[i].PNum(6)];

          if(MESH_DEBUG)
          {
            cout << locelements[i].PNum(1) <<  " (" <<      locpoints[locelements[i].PNum(1)].X() << "," <<
            locpoints[locelements[i].PNum(1)].Y() << "," <<
            locpoints[locelements[i].PNum(1)].Z() << ")" << endl;
            cout << locelements[i].PNum(2) <<  " (" <<      locpoints[locelements[i].PNum(2)].X() << "," <<
            locpoints[locelements[i].PNum(2)].Y() << "," <<
            locpoints[locelements[i].PNum(2)].Z() << ")" << endl;
            cout << locelements[i].PNum(3) <<  " (" <<      locpoints[locelements[i].PNum(3)].X() << "," <<
            locpoints[locelements[i].PNum(3)].Y() << "," <<
            locpoints[locelements[i].PNum(3)].Z() << ")" << endl;
            cout << locelements[i].PNum(4) <<  " (" <<      locpoints[locelements[i].PNum(4)].X() << "," <<
            locpoints[locelements[i].PNum(4)].Y() << "," <<
            locpoints[locelements[i].PNum(4)].Z() << ")" << endl;
            cout << locelements[i].PNum(5) <<  " (" <<      locpoints[locelements[i].PNum(5)].X() << "," <<
            locpoints[locelements[i].PNum(5)].Y() << "," <<
            locpoints[locelements[i].PNum(5)].Z() << ")" << endl;
            cout << locelements[i].PNum(6) <<  " (" <<      locpoints[locelements[i].PNum(6)].X() << "," <<
            locpoints[locelements[i].PNum(6)].Y() << "," <<
            locpoints[locelements[i].PNum(6)].Z() << ")" << endl;

            printf("%d %d %d %d %d %d\n",   locelements[i].PNum(1),
                   locelements[i].PNum(2),
                   locelements[i].PNum(3),
                   locelements[i].PNum(4),
                   locelements[i].PNum(5),
                   locelements[i].PNum(6)
                   );
            printf("%s %d\n","nff:",adfront->NFF());
          }
          prismvol = prismvol
                     + (1.0 / 6.0) * ( Cross ( *hp2 - *hp1, *hp3 - *hp1) * (*hp4 - *hp1) )
                     + (1.0 / 6.0) * ( Cross ( *hp4 - *hp3, *hp5 - *hp3) * (*hp6 - *hp3) )
                     + (1.0 / 6.0) * ( Cross ( *hp3 - *hp2, *hp4 - *hp2) * (*hp5 - *hp2) );

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
          adfront->AddFace (locfaces[i],-1);
        }

        for (i = 1; i <= delfaces.Size(); i++)
          adfront->DeleteFace (findex[delfaces[i]]);
      }
      else
        adfront->ResetPrism (findex[1]);

      locelements.SetSize (0);
      delpoints.SetSize(0);
      delfaces.SetSize(0);
      prism_flags.SetSize(0);
    }

  //adfront->Print();

  // close with pyramids, if nessesary
  if(gen_pyramid)
  {
    while (!adfront -> Empty())
    {
      //adfront->Print();
      h = ah;
      locpoints.SetSize(0);
      locfaces.SetSize(0);
      pindex.SetSize(0);
      findex.SetSize(0);

      pyramid = adfront->GetLocals_Pyramid(locpoints, locfaces, pindex, findex, 3 * h);
      if(pyramid==-1)
      {
        gen_pyramid = 0;
        gen_tetrahedra = 1;
        break;
      }

      oldnp = locpoints.Size();
      oldnf = locfaces.Size();
      found = 0;
      hasfound = 0;
      minerr = 1000;

      trans.Set (     locpoints[locfaces[1].PNumMod(1)],
                      locpoints[locfaces[1].PNumMod(2)],
                      locpoints[locfaces[1].PNumMod(3)], h);

      trans.ToPlain (locpoints, plainpoints);

      cnttrials++;

      found = Generate_Pyramid (      plainpoints, locfaces, locelements,
                                      delfaces, prism_flags);

      if (found) cntsucc++;

      locpoints.SetSize (plainpoints.Size());
      for (i = oldnp+1; i <= plainpoints.Size(); i++)
        trans.FromPlain (plainpoints[i], locpoints[i]);

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

      locpoints.SetSize (oldnp);
      locfaces.SetSize (oldnf);
      delfaces.SetSize (0);
      locelements.SetSize (0);


      if (found)
      {
        for (i = 1; i <= tempnewpoints.Size(); i++)
          locpoints.Append (tempnewpoints[i]);
        for (i = 1; i <= tempnewfaces.Size(); i++)
          locfaces.Append (tempnewfaces[i]);
        for (i = 1; i <= tempdelfaces.Size(); i++)
          delfaces.Append (tempdelfaces[i]);
        for (i = 1; i <= templocelements.Size(); i++)
          locelements.Append (templocelements[i]);

        pindex.SetSize(locpoints.Size());

        for (i = oldnp+1; i <= locpoints.Size(); i++)
        {
          globind = SavePoint (locpoints[i]);
          pindex[i] = adfront -> AddPoint (locpoints[i], globind);
        }

        for (i = 1; i <= locelements.Size(); i++)
        {
          hp1 = &locpoints[locelements[i].PNum(1)];
          hp2 = &locpoints[locelements[i].PNum(2)];
          hp3 = &locpoints[locelements[i].PNum(3)];
          hp4 = &locpoints[locelements[i].PNum(4)];
          hp5 = &locpoints[locelements[i].PNum(5)];

          if(MESH_DEBUG)
          {
            cout << locelements[i].PNum(1) <<  " (" <<      locpoints[locelements[i].PNum(1)].X() << "," <<
            locpoints[locelements[i].PNum(1)].Y() << "," <<
            locpoints[locelements[i].PNum(1)].Z() << ")" << endl;
            cout << locelements[i].PNum(2) <<  " (" <<      locpoints[locelements[i].PNum(2)].X() << "," <<
            locpoints[locelements[i].PNum(2)].Y() << "," <<
            locpoints[locelements[i].PNum(2)].Z() << ")" << endl;
            cout << locelements[i].PNum(3) <<  " (" <<      locpoints[locelements[i].PNum(3)].X() << "," <<
            locpoints[locelements[i].PNum(3)].Y() << "," <<
            locpoints[locelements[i].PNum(3)].Z() << ")" << endl;
            cout << locelements[i].PNum(4) <<  " (" <<      locpoints[locelements[i].PNum(4)].X() << "," <<
            locpoints[locelements[i].PNum(4)].Y() << "," <<
            locpoints[locelements[i].PNum(4)].Z() << ")" << endl;
            cout << locelements[i].PNum(5) <<  " (" <<      locpoints[locelements[i].PNum(5)].X() << "," <<
            locpoints[locelements[i].PNum(5)].Y() << "," <<
            locpoints[locelements[i].PNum(5)].Z() << ")" << endl;

            printf("%d %d %d %d %d\n",      locelements[i].PNum(1),
                   locelements[i].PNum(2),
                   locelements[i].PNum(3),
                   locelements[i].PNum(4),
                   locelements[i].PNum(5)
                   );
            printf("%s %d\n","nff:",adfront->NFF());
          }
          pyramidvol = pyramidvol
                       + (1.0 / 6.0) * ( Cross ( *hp2 - *hp1, *hp4 - *hp1) * (*hp5 - *hp1) )
                       + (1.0 / 6.0) * ( Cross ( *hp3 - *hp2, *hp4 - *hp2) * (*hp5 - *hp2) );

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

          adfront->AddFace (locfaces[i],-1);
        }

        for (i = 1; i <= delfaces.Size(); i++)
          adfront->DeleteFace (findex[delfaces[i]]);
      }
      else
      {
        printf("%s\n", "Can not create pyramid");
        assert(0);
      }
      locelements.SetSize (0);
      delpoints.SetSize(0);
      delfaces.SetSize(0);
    }
  }
  if(make_prism)
    for(i=1; i<=adfront->NFF(); i++)
      adfront->SetClass(i,1);

  // and now the tetrahedra
  while (!adfront -> Empty())
  {
    //adfront->Print();
    locpoints.SetSize(0);
    locfaces.SetSize(0);
    pindex.SetSize(0);
    findex.SetSize(0);

    h = ah;
    if(h>0.0)
      qualclass = adfront -> GetLocals_Tetrahedra (locpoints, locfaces, pindex, findex, 3 * h);
    else
    {
      qualclass = adfront -> GetLocals_Tetrahedra (locpoints, locfaces, pindex, findex, -3 * h);

      bemp1.X() = locpoints[locfaces[1].PNum(1)].X();
      bemp1.Y() = locpoints[locfaces[1].PNum(1)].Y();
      bemp1.Z() = locpoints[locfaces[1].PNum(1)].Z();
      bemp2.X() = locpoints[locfaces[1].PNum(2)].X();
      bemp2.Y() = locpoints[locfaces[1].PNum(2)].Y();
      bemp2.Z() = locpoints[locfaces[1].PNum(2)].Z();
      bemp3.X() = locpoints[locfaces[1].PNum(3)].X();
      bemp3.Y() = locpoints[locfaces[1].PNum(3)].Y();
      bemp3.Z() = locpoints[locfaces[1].PNum(3)].Z();

      bemp.X() = ( bemp1.X() + bemp2.X() + bemp3.X() ) / 3;
      bemp.Y() = ( bemp1.Y() + bemp2.Y() + bemp3.Y() ) / 3;
      bemp.Z() = ( bemp1.Z() + bemp2.Z() + bemp3.Z() ) / 3;

      locpoints.SetSize(0);
      locfaces.SetSize(0);
      pindex.SetSize(0);
      findex.SetSize(0);

      in[0] = bemp.X();
      in[1] = bemp.Y();
      in[2] = bemp.Z();
      in[3] = ah;
      Get_Local_h_3d(in,&h);
      qualclass = adfront -> GetLocals_Tetrahedra (locpoints, locfaces, pindex, findex, 3 * h);
    }

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
      trans.Set (     locpoints[locfaces[1].PNumMod(1+rotind)],
                      locpoints[locfaces[1].PNumMod(2+rotind)],
                      locpoints[locfaces[1].PNumMod(3+rotind)], h);

      trans.ToPlain (locpoints, plainpoints);

      cnttrials++;

      found = ApplyVRules (   rules, tolfak, plainpoints, locfaces, locelements,
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
          locelements[i].PNum(j) = adfront -> GetGlobalIndex (pindex[locelements[i].PNum(j)]);
        //			if(MESH_DEBUG)
        printf("%s %d\n","NFF:",adfront->NFF());

        SaveElement (locelements[i]);
        cntelem++;
      }

      for (i = oldnf+1; i <= locfaces.Size(); i++)
      {
        for (j = 1; j <= locfaces[i].NP(); j++)
          locfaces[i].PNum(j) = pindex[locfaces[i].PNum(j)];

        adfront->AddFace (locfaces[i], -1);
      }

      for (i = 1; i <= delfaces.Size(); i++)
      {
        adfront->DeleteFace (findex[delfaces[i]]);
      }

    }
    else
    {
      adfront->IncrementClass (findex[1]);
    }

    locelements.SetSize (0);
    delpoints.SetSize(0);
    delfaces.SetSize(0);
  }

}
