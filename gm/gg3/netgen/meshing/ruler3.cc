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

static double CalcElementBadness (const ARRAY<Point3d> & points,
                                  const Element & elem)
{
  double vol, l, l4, l5, l6;

  Vec3d v1 = points[elem.PNum(2)] - points[elem.PNum(1)];
  Vec3d v2 = points[elem.PNum(3)] - points[elem.PNum(1)];
  Vec3d v3 = points[elem.PNum(4)] - points[elem.PNum(1)];

  vol = - (Cross(v1, v2) * v3);
  l4 = Dist (points[elem.PNum(2)], points[elem.PNum(3)]);
  l5 = Dist (points[elem.PNum(2)], points[elem.PNum(4)]);
  l6 = Dist (points[elem.PNum(3)], points[elem.PNum(4)]);

  l = v1.Length() + v2.Length() + v3.Length() + l4 + l5 + l6;

  if (vol < 1e-8)
    return 1e10;

  return (l*l*l/vol)/500;
}

int ApplyVRules (const ARRAY<vnetrule*> & rules, double tolfak,
                 ARRAY<Point3d> & lpoints, ARRAY<Element> & lfaces,
                 ARRAY<Element> & elements,
                 ARRAY<INDEX> & delfaces, int tolerance, int rotind1,
                 float & retminerr, ARRAY<char*> & problems)
{
  return(0);
}

int ApplyVRules_new
(
  const ARRAY<vnetrule_new*> & rules, double tolfak,
  ARRAY<Point3d> & lpoints,    // in: local points, out: old+new local points
  ARRAY<Element> & lfaces,  // in: local faces, out: old+new local faces
  INDEX lfacesplit,            // for local faces in outer radius
  ARRAY<Element> & elements,   // out: new elements
  ARRAY<INDEX> & delfaces,     // out: face indices of faces to delete
  int tolerance,               // quality class: 1 best
  int rotind1,                 // how to rotate base element
  float & retminerr,           // element error
  ARRAY<char*> & problems
)

{
  const int rotsym = 3;
  int i, j, k, ri, nfok, npok, incnpok, refpi, locpi, locfi, locfr;
  int hi, minn, hpi;
  float hf, err, minerr, teterr, minteterr;
  char ok, found, hc;
  vnetrule_new * rule;
  Vector oldu, newu, newu1, newu2;
  Vec3d ui;
  Point3d np;
  int oldnp, noldlp, noldlf;
  const Element * locface = NULL;
  int loktestmode;
  Element own;

  static ARRAY<int> pused;       // point is already mapped
  static ARRAY<int> fused;       // face is already mapped
  static ARRAY<int> pmap;        // map of reference point to local point
  static ARRAY<int> pfixed;      // ???
  static ARRAY<int> fmapi;       // face in reference is mapped to face nr ...
  static ARRAY<int> fmapr;       // face in reference is rotated to map
  static ARRAY<Point3d> transfreezone;  // transformed free-zone
  static INDEX cnt = 0;
  INDEX_2_HASHTABLE<int> ledges(lfaces.Size());  // edges in local environment

  static ARRAY<Point3d> tempnewpoints;
  static ARRAY<Element> tempnewfaces;
  static ARRAY<int> tempdelfaces;
  static ARRAY<Element> tempelements;

  static ARRAY<double> triminx, trimaxx, triminy, trimaxy, triminz, trimaxz;

  //  static ARRAY<int> notnewnewfaces;

  static ARRAY<int> pnearness;
  static ARRAY<int> fnearness;

  cnt++;

  delfaces.SetSize (0);
  elements.SetSize (0);


  // determine topological distance of faces and points to
  // base element

  pnearness.SetSize (lpoints.Size());
  fnearness.SetSize (lfacesplit);

  for (i = 1; i <= pnearness.Size(); i++)
    pnearness.Set(i, INT_MAX/4);

  for (j = 1; j <= 3; j++)
    pnearness.Set(lfaces.Get(1).PNum(j), 0);

  do
  {
    ok = 1;

    for (i = 1; i <= lfacesplit; i++)
    //      for (i = 1; i <= lfaces.Size(); i++)
    {
      const Element & hface = lfaces.Get(i);

      minn = INT_MAX-1;
      for (j = 1; j <= 3; j++)
      {
        hi = pnearness.Get(hface.PNum(j));
        if (hi < minn) minn = hi;
      }

      for (j = 1; j <= 3; j++)
      {
        hpi = hface.PNum(j);
        if (pnearness.Get(hpi) > minn+1)
        {
          ok = 0;
          pnearness.Set(hpi, minn+1);
        }
      }
    }
  }
  while (!ok);

  for (i = 1; i <= fnearness.Size(); i++)
  {
    fnearness.Set(i, 0);
    for (j = 1; j <= 3; j++)
      fnearness.Elem(i) += pnearness.Get(lfaces.Get(i).PNum(j));
  }


  triminx.SetSize (lfaces.Size());
  trimaxx.SetSize (lfaces.Size());
  triminy.SetSize (lfaces.Size());
  trimaxy.SetSize (lfaces.Size());
  triminz.SetSize (lfaces.Size());
  trimaxz.SetSize (lfaces.Size());

  for (i = 1; i <= lfaces.Size(); i++)
  {
    const Element & face = lfaces.Get(i);
    const Point3d & p1 = lpoints.Get(face.PNum(1));
    const Point3d & p2 = lpoints.Get(face.PNum(2));
    const Point3d & p3 = lpoints.Get(face.PNum(3));
    triminx.Elem(i) = min (p1.X(), p2.X(), p3.X());
    trimaxx.Elem(i) = max (p1.X(), p2.X(), p3.X());
    triminy.Elem(i) = min (p1.Y(), p2.Y(), p3.Y());
    trimaxy.Elem(i) = max (p1.Y(), p2.Y(), p3.Y());
    triminz.Elem(i) = min (p1.Z(), p2.Z(), p3.Z());
    trimaxz.Elem(i) = max (p1.Z(), p2.Z(), p3.Z());
  }



  for (j = 1; j <= lfacesplit; j++)
  {
    INDEX_2 i2;
    own = lfaces.Get(j);
    i2.I1() = lfaces.Get(j).PNumMod(1);
    i2.I2() = lfaces.Get(j).PNumMod(2);
    ledges.Set (i2, 1);
    i2.I1() = lfaces.Get(j).PNumMod(2);
    i2.I2() = lfaces.Get(j).PNumMod(3);
    ledges.Set (i2, 1);
    i2.I1() = lfaces.Get(j).PNumMod(3);
    i2.I2() = lfaces.Get(j).PNumMod(1);
    ledges.Set (i2, 1);
  }




  pused.SetSize (lpoints.Size());
  fused.SetSize (lfaces.Size());



  found = 0;
  minerr = tolfak * tolerance * tolerance;
  minteterr = 500 * tolerance * tolerance;

  // check each rule:

  for (ri = 1; ri <= rules.Size(); ri++)
  {   // loop ri


    rule = rules.Get(ri);

    if (rule->GetQuality() > tolerance)
    {
      if (testmode)
        sprintf (problems[ri], "Quality not ok");
      continue;
    }

    if (testmode)
      sprintf (problems[ri], "no mapping found");

    loktestmode = testmode || rule->TestFlag ('t');

    pmap.SetSize (rule->GetNP());
    fmapi.SetSize (rule->GetNF());
    fmapr.SetSize (rule->GetNF());

    for (i = 1; i <= lfaces.Size(); i++)
      fused.Set (i, 0);
    for (i = 1; i <= lpoints.Size(); i++)
      pused.Set (i, 0);
    for (i = 1; i <= pmap.Size(); i++)
      pmap.Set(i, 0);
    for (i = 1; i <= fmapi.Size(); i++)
      fmapi.Set(i, 0);
    for (i = 1; i <= fmapr.Size(); i++)
      fmapr.Set(i, 3);

    fused.Set (1, 1);

    fmapi.Set (1, 1);
    fmapr.Set (1, rotind1);


    for (j = 1; j <= 3; j++)
    {
      locpi = lfaces.Get(1).PNumMod (j+rotind1);
      pmap.Set (rule->GetPointNr (1, j), locpi);
      pused.Elem(locpi)++;
    }



    /*
       map all faces
       nfok .. first nfok-1 faces are mapped properly
     */

    nfok = 2;
    while (nfok >= 2)
    {

      if (nfok <= rule->GetNOldF())
      {
        // not all faces mapped

        ok = 0;
        locfi = fmapi.Get(nfok);
        locfr = fmapr.Get(nfok);

        while (!ok)
        {
          locfr++;
          if (locfr == rotsym + 1)
          {
            locfr = 1;
            locfi++;
            if (locfi > lfacesplit) break;
          }


          if (fnearness.Get(locfi) > rule->GetFNearness (nfok) ||
              fused.Get(locfi) )
          {
            // face not feasible in any rotation

            locfr = rotsym;
          }
          else
          {

            ok = 1;

            locface = &lfaces.Get(locfi);

            // reference point already mapped differently ?
            for (j = 1; j <= 3 && ok; j++)
            {
              locpi = pmap.Get(rule->GetPointNr (nfok, j));

              if (locpi && locpi != locface->PNumMod(j+locfr))
                ok = 0;
            }

            // local point already used or point outside tolerance ?
            for (j = 1; j <= 3 && ok; j++)
            {
              refpi = rule->GetPointNr (nfok, j);

              if (pmap.Get(refpi) == 0)
              {
                locpi = locface->PNumMod (j + locfr);

                /*
                   if (pused.Get(locpi) ||
                    rule->CalcPointDist (refpi, lpoints.Get(locpi)) > minerr)
                   ok = 0;
                 */

                if (pused.Get(locpi))
                  ok = 0;
                else
                {
                  const Point3d & lp = lpoints.Get(locpi);
                  const Point3d & rp = rule->GetPoint(refpi);
                  if ( ( (lp.X()-rp.X())*(lp.X()-rp.X()) +
                         (lp.Y()-rp.Y())*(lp.Y()-rp.Y()) +
                         (lp.Z()-rp.Z())*(lp.Z()-rp.Z()) )
                       * rule->PointDistFactor(refpi) > minerr)
                    //				  if (Dist2(lp, rp) * rule->PointDistFactor(refpi) > minerr)
                    ok = 0;
                }
              }
            }
          }
        }


        if (ok)
        {
          // map face nfok

          fmapi.Set (nfok, locfi);
          fmapr.Set (nfok, locfr);
          fused.Set (locfi, 1);

          for (j = 1; j <= rule->GetNP (nfok); j++)
          {
            locpi = locface->PNumMod(j+locfr);

            if (rule->GetPointNr (nfok, j) <= 3 &&
                pmap.Get(rule->GetPointNr(nfok, j)) != locpi)
              cout << "change face1 point, mark1" << endl;

            pmap.Set(rule->GetPointNr (nfok, j), locpi);
            pused.Elem(locpi)++;
          }

          nfok++;
        }
        else
        {
          // backtrack one face
          fmapi.Set (nfok, 0);
          fmapr.Set (nfok, 3);
          nfok--;

          fused.Set (fmapi.Get(nfok), 0);
          for (j = 1; j <= rule->GetNP (nfok); j++)
          {
            refpi = rule->GetPointNr (nfok, j);
            pused.Elem(pmap.Get(refpi))--;

            if (pused.Get(pmap.Get(refpi)) == 0)
            {
              pmap.Set(refpi, 0);
            }
          }
        }
      }

      else

      {

        // all faces are mapped
        // now map all isolated points:

        if (loktestmode)
        {
          sprintf (problems[ri], "Faces Ok");
        }

        npok = 1;
        incnpok = 1;

        pfixed.SetSize (pmap.Size());
        for (i = 1; i <= pmap.Size(); i++)
          pfixed.Set(i, (pmap.Get(i) != 0) );

        while (npok >= 1)
        {

          if (npok <= rule->GetNOldP())
          {

            if (pfixed.Get(npok))

            {
              if (incnpok)
                npok++;
              else
                npok--;
            }

            else

            {
              locpi = pmap.Elem(npok);
              ok = 0;

              if (locpi)
                pused.Elem(locpi)--;

              while (!ok && locpi < lpoints.Size())
              {
                ok = 1;
                locpi++;

                if (pused.Get(locpi))
                {
                  ok = 0;
                }
                else
                {
                  // if (rule->CalcPointDist (npok, lpoints.Get(locpi)) > minerr)
                  //    ok = 0;

                  const Point3d & lp = lpoints.Get(locpi);
                  const Point3d & rp = rule->GetPoint(npok);
                  if ( ( (lp.X()-rp.X())*(lp.X()-rp.X()) +
                         (lp.Y()-rp.Y())*(lp.Y()-rp.Y()) +
                         (lp.Z()-rp.Z())*(lp.Z()-rp.Z()) )
                       * rule->PointDistFactor(npok) > minerr)
                    ok = 0;
                }
              }


              if (ok)
              {
                pmap.Set (npok, locpi);

                pused.Elem(locpi)++;
                npok++;
                incnpok = 1;
              }

              else

              {
                pmap.Set (npok, 0);

                npok--;
                incnpok = 0;
              }
            }
          }

          else

          {

            // all points are mapped

            if (loktestmode)
            {
              sprintf (problems[ri], "mapping found");
            }

            ok = 1;


            // check mapedges:



            for (i = 1; i <= rule->GetNEd(); i++)
            {

              int i1, i2;
              i1 = pmap.Get(rule->GetEdge(i).i1);
              i2 = pmap.Get(rule->GetEdge(i).i2);

              INDEX_2 in2(i1, i2);
              if (!ledges.Used (in2)) ok = 0;
            }


            for (i = rule->GetNOldF() + 1; i <= rule->GetNF(); i++)
              fmapi.Set(i, 0);



            // deviation of existing points

            oldu.SetLength (3 * rule->GetNOldP());

            for (i = 1; i <= rule->GetNOldP(); i++)
            {
              //			  ui = lpoints.Get(pmap.Get(i)) - rule->GetPoint(i);
              const Point3d & lp = lpoints.Get(pmap.Get(i));
              const Point3d & rp = rule->GetPoint(i);
              oldu.Set (3*i-2, lp.X()-rp.X());
              oldu.Set (3*i-1, lp.Y()-rp.Y());
              oldu.Set (3*i  , lp.Z()-rp.Z());
            }


            //		      rule->GetOldUToFreeZone().Mult (oldu, newu);
            //		      rule->SetFreeZoneTransformation (newu);
            rule->SetFreeZoneTransformation (oldu, tolerance);

            if (!rule->ConvexFreeZone())
            {
              if (loktestmode)
              {
                sprintf (problems[ri], "Freezone not convex");
              }
              ok = 0;
            }


            // check freezone:

            for (i = 1; i <= lpoints.Size(); i++)
            {
              if ( !pused.Get(i) )
              {
                const Point3d & lp = lpoints.Get(i);

                if (lp.X() >= rule->fzminx && lp.X() <= rule->fzmaxx &&
                    lp.Y() >= rule->fzminy && lp.Y() <= rule->fzmaxy &&
                    lp.Z() >= rule->fzminz && lp.Z() <= rule->fzmaxz)
                {
                  if (rule->IsInFreeZone(lpoints.Get(i)))
                  {
                    if (loktestmode)
                    {
                      sprintf (problems[ri],
                               "locpoint %d in Freezone", i);
                    }
                    ok = 0;
                    break;
                  }
                }
              }
            }


            for (i = 1; i <= lfaces.Size() && ok; i++)
            {
              if (!fused.Get(i))
              {
                int triin;

                if (triminx.Elem(i) > rule->fzmaxx ||
                    trimaxx.Elem(i) < rule->fzminx ||
                    triminy.Elem(i) > rule->fzmaxy ||
                    trimaxy.Elem(i) < rule->fzminy ||
                    triminz.Elem(i) > rule->fzmaxz ||
                    trimaxz.Elem(i) < rule->fzminz)
                  triin = 0;
                else
                  triin = rule->IsTriangleInFreeZone
                          (
                    lpoints.Get(lfaces.Get(i).PNum(1)),
                    lpoints.Get(lfaces.Get(i).PNum(2)),
                    lpoints.Get(lfaces.Get(i).PNum(3))
                          );

                if (triin == -1)
                {
                  ok = 0;
                }

                if (triin == 1)
                {
                  hc = 0;
                  for (k = rule->GetNOldF() + 1; k <= rule->GetNF(); k++)
                  {
                    if (rule->GetPointNr(k, 1) <= rule->GetNOldP() &&
                        rule->GetPointNr(k, 2) <= rule->GetNOldP() &&
                        rule->GetPointNr(k, 3) <= rule->GetNOldP())
                    {
                      for (j = 1; j <= 3; j++)
                        if (lfaces.Get(i).PNumMod(j  ) == pmap.Get(rule->GetPointNr(k, 1)) &&
                            lfaces.Get(i).PNumMod(j+1) == pmap.Get(rule->GetPointNr(k, 3)) &&
                            lfaces.Get(i).PNumMod(j+2) == pmap.Get(rule->GetPointNr(k, 2)))
                        {
                          fmapi.Elem(k) = i;
                          hc = 1;
                        }
                    }
                  }

                  if (!hc)
                  {
                    if (loktestmode)
                    {
                      sprintf (problems[ri], "triangle (%d, %d, %d) in Freezone",
                               lfaces.Get(i).PNum(1), lfaces.Get(i).PNum(2),
                               lfaces.Get(i).PNum(3));
                    }
                    ok = 0;
                  }
                }
              }
            }


            if (ok)
            {

              err = 0;
              for (i = 1; i <= rule->GetNOldP(); i++)
              {
                hf = rule->CalcPointDist (i, lpoints[pmap[i]]);
                if (hf > err) err = hf;
              }

              if (loktestmode)
              {
                sprintf (problems[ri], "Rule ok, err = %f", err);
              }

              newu = rule->GetOldUToNewU() * oldu;

              // set new points:

              oldnp = rule->GetNOldP();
              noldlp = lpoints.Size();
              noldlf = lfaces.Size();


              for (i = oldnp + 1; i <= rule->GetNP(); i++)
              {
                np = rule->GetPoint(i);
                np.X() += newu (3 * (i-oldnp) - 2);
                np.Y() += newu (3 * (i-oldnp) - 1);
                np.Z() += newu (3 * (i-oldnp));

                pmap[i] = lpoints.Append (np);
              }

              // Set new Faces:

              for (i = rule->GetNOldF() + 1; i <= rule->GetNF(); i++)
                if (!fmapi[i])
                {
                  Element nface;
                  for (j = 1; j <= 3; j++)
                    nface.PNum(j) = pmap[rule->GetPointNr (i, j)];
                  nface.SetNP(3);
                  lfaces.Append (nface);
                }


              // Delete old Faces:

              for (i = 1; i <= rule->GetNDelF(); i++)
                delfaces.Append (fmapi.Get(rule->GetDelFace(i)));
              for (i = rule->GetNOldF()+1; i <= rule->GetNF(); i++)
                if (fmapi.Get(i))
                {
                  delfaces.Append (fmapi[i]);
                  fmapi[i] = 0;
                }


              for (i = 1; i <= rule->GetNO() && ok; i++)
              {
                const fourint * fouri;
                double v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;

                fouri = &rule->GetOrientation(i);

                v1x = lpoints.Get(pmap.Get(fouri->i2)).X() -
                      lpoints.Get(pmap.Get(fouri->i1)).X();
                v1y = lpoints.Get(pmap.Get(fouri->i2)).Y() -
                      lpoints.Get(pmap.Get(fouri->i1)).Y();
                v1z = lpoints.Get(pmap.Get(fouri->i2)).Z() -
                      lpoints.Get(pmap.Get(fouri->i1)).Z();

                v2x = lpoints.Get(pmap.Get(fouri->i3)).X() -
                      lpoints.Get(pmap.Get(fouri->i1)).X();
                v2y = lpoints.Get(pmap.Get(fouri->i3)).Y() -
                      lpoints.Get(pmap.Get(fouri->i1)).Y();
                v2z = lpoints.Get(pmap.Get(fouri->i3)).Z() -
                      lpoints.Get(pmap.Get(fouri->i1)).Z();

                v3x = lpoints.Get(pmap.Get(fouri->i4)).X() -
                      lpoints.Get(pmap.Get(fouri->i1)).X();
                v3y = lpoints.Get(pmap.Get(fouri->i4)).Y() -
                      lpoints.Get(pmap.Get(fouri->i1)).Y();
                v3z = lpoints.Get(pmap.Get(fouri->i4)).Z() -
                      lpoints.Get(pmap.Get(fouri->i1)).Z();


                if (v1x * v2y * v3z +
                    v1y * v2z * v3x +
                    v1z * v2x * v3y -
                    v1x * v2z * v3y -
                    v1y * v2x * v3z -
                    v1z * v2y * v3x > -1E-7)
                {
                  if (loktestmode)
                  {
                    sprintf (problems[ri], "Orientation wrong");
                  }
                  ok = 0;
                }
              }


              // new points in free-zone ?
              for (i = rule->GetNOldP() + 1; i <= rule->GetNP() && ok; i++)
                if (!rule->IsInFreeZone (lpoints[pmap[i]]))
                {
                  if (loktestmode)
                  {
                    sprintf (problems[ri], "newpoint outside convex hull");
                  }
                  ok = 0;

                }

              // insert new elements

              for (i = 1; i <= rule->GetNE(); i++)
              {
                elements.Append (rule->GetElement(i));
                for (j = 1; j <= elements[i].NP(); j++)
                  elements[i].PNum(j) = pmap[elements[i].PNum(j)];
              }


              // Calculate Element badness

              teterr = 0;
              for (i = 1; i <= elements.Size(); i++)
              {
                hf = CalcElementBadness (lpoints, elements[i]);
                if (hf > teterr) teterr = hf;
              }

              if (ok && teterr < 1e6 &&
                  (rule->TestFlag('b') || tolerance > 10) )
              {
                //  cout << "Reset teterr "
                //   << rule->Name()
                //   << " err = " << teterr
                //   << endl;
                printf("Reset teterr %s err = %g\n",rule->Name(),teterr);
                teterr = 1;
              }

              // compare edgelength
              if (rule->TestFlag('l'))
              {
                double oldlen = 0;
                double newlen = 0;

                for (i = 1; i <= rule->GetNDelF(); i++)
                {
                  const Element & face =
                    rule->GetFace (rule->GetDelFace(i));
                  for (j = 1; j <= 3; j++)
                  {
                    const Point3d & p1 =
                      lpoints.Get(pmap.Get(face.PNumMod(j)));
                    const Point3d & p2 =
                      lpoints.Get(pmap.Get(face.PNumMod(j+1)));
                    oldlen += Dist(p1, p2);
                  }
                }

                for (i = rule->GetNOldF()+1; i <= rule->GetNF(); i++)
                {
                  const Element & face = rule->GetFace (i);
                  for (j = 1; j <= 3; j++)
                  {
                    const Point3d & p1 =
                      lpoints.Get(pmap.Get(face.PNumMod(j)));
                    const Point3d & p2 =
                      lpoints.Get(pmap.Get(face.PNumMod(j+1)));
                    newlen += Dist(p1, p2);
                  }
                }

                if (oldlen < newlen)
                {
                  ok = 0;
                  if (loktestmode)
                    sprintf (problems[ri], "oldlen < newlen");
                }
              }

              if (ok && teterr < minteterr)
              {
                found = ri;
                minteterr = teterr;

                if (testmode)
                {}

                tempnewpoints.SetSize (0);
                for (i = noldlp+1; i <= lpoints.Size(); i++)
                  tempnewpoints.Append (lpoints[i]);

                tempnewfaces.SetSize (0);
                for (i = noldlf+1; i <= lfaces.Size(); i++)
                  tempnewfaces.Append (lfaces[i]);

                tempdelfaces.SetSize (0);
                for (i = 1; i <= delfaces.Size(); i++)
                  tempdelfaces.Append (delfaces[i]);

                tempelements.SetSize (0);
                for (i = 1; i <= elements.Size(); i++)
                  tempelements.Append (elements[i]);
              }

              lpoints.SetSize (noldlp);
              lfaces.SetSize (noldlf);
              delfaces.SetSize (0);
              elements.SetSize (0);
            }

            npok = rule->GetNOldP();
            incnpok = 0;
          }
        }

        nfok = rule->GetNOldF();

        for (j = 1; j <= rule->GetNP (nfok); j++)
        {
          refpi = rule->GetPointNr (nfok, j);
          pused.Elem(pmap.Get(refpi))--;

          if (pused.Get(pmap.Get(refpi)) == 0)
          {
            pmap.Set(refpi, 0);
          }
        }

      }
    }
  }


  // if successfull, reload best choice

  if (found)
  {

#ifdef debug
    // if face in advancing front ???
    for (i = 1; i <= tempnewfaces.Size(); i++)
    {
      hc = 1;
      for (k = 1; k <= lfaces.Size() && hc; k++)
        for (j = 1; j <= 3 && hc; j++)
          if (tempnewfaces.Elem(i).PNumMod(j  ) == lfaces.Get(k).PNum(1) &&
              tempnewfaces.Elem(i).PNumMod(j+1) == lfaces.Get(k).PNum(3) &&
              tempnewfaces.Elem(i).PNumMod(j+2) == lfaces.Get(k).PNum(2))
          {
            tempdelfaces.Append(k);
            tempnewfaces.Elem(i).PNum(1) = 0;
            hc = 0;
            cerr << "Ruler-reload necessary" << endl;
          }
    }
#endif

    for (i = 1; i <= tempnewpoints.Size(); i++)
      lpoints.Append (tempnewpoints.Get(i));
    for (i = 1; i <= tempnewfaces.Size(); i++)
      if (tempnewfaces.Get(i).PNum(1))
        lfaces.Append (tempnewfaces.Get(i));
    for (i = 1; i <= tempdelfaces.Size(); i++)
      delfaces.Append (tempdelfaces.Get(i));
    for (i = 1; i <= tempelements.Size(); i++)
      elements.Append (tempelements.Get(i));
  }

  retminerr = minerr;
  return found;
}
