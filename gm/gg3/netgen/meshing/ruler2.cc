// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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


static double CalcElementBadness_old (const ARRAY<Point2d> & points,
                                      const Element & elem)
{
  // badness = sqrt(3) /36 * circumference^2 / area - 1 +
  //           h / li + li / h - 2

  Vec2d v12, v13, v23;
  double l12, l13, l23, cir, area;
  static const double c = sqrt(3) / 36;

  v12 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(1));
  v13 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(1));
  v23 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(2));

  l12 = v12.Length();
  l13 = v13.Length();
  l23 = v23.Length();

  cir = l12 + l13 + l23;
  area = 0.5 * (v12.X() * v13.Y() - v12.Y() * v13.X());
  if (area < 1e-6)
  {
    return 1e8;
  }

  return c * cir * cir / area - 1
         + 1/l12 + l12 + 1/l13 + l13 + 1/l23 + l23 - 6;
}

static double CalcElementBadness (const ARRAY<Point2d> & points,
                                  const Element & elem)
{
  Vec2d v12, v13, v23, v21, v31, v32,v;
  double l12, l13, l23, l21, l31, l32, cir, area;
  double alpha1, alpha2, alpha3, maxalpha, minalpha,max_l,l;
  double v1x,v1y,v1z,v2x,v2y,v2z;
  Vec3d n;

  v12 = points.Get(elem.PNum(1)) - points.Get(elem.PNum(2));
  v13 = points.Get(elem.PNum(1)) - points.Get(elem.PNum(3));
  v23 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(3));

  v21 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(1));
  v31 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(1));
  v32 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(2));

  l12 = v12.Length();
  l13 = v13.Length();
  l23 = v23.Length();

  l21 = v21.Length();
  l31 = v31.Length();
  l32 = v32.Length();

  alpha1 = acos( v21*v31 / (l21*l31) );
  alpha2 = acos( v12*v32 / (l12*l32) );
  alpha3 = acos( v23*v13 / (l23*l13) );

  maxalpha = alpha1;
  if(maxalpha<alpha2)
    maxalpha = alpha2;
  if(maxalpha<alpha3)
    maxalpha = alpha3;

  minalpha = alpha1;
  if(minalpha>alpha2)
    minalpha = alpha2;
  if(minalpha>alpha3)
    minalpha = alpha3;

  if(minalpha<1e-6)
    return(1e10);
  return( 0.5*(0.5-cos(maxalpha)) + 0.5*(max_l/l-0.707)/2 );

  max_l = l12;
  if(max_l<l13)
    max_l = l13;
  if(max_l<l23)
    max_l = l23;

  v = points.Get(1) - points.Get(2);
  l = v.Length();

  v12 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(1));
  v13 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(1));
  v23 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(2));

  area = 0.5 * (v12.X() * v13.Y() - v12.Y() * v13.X());
  if (area < 1e-6)
  {
    return 1e8;
  }

  //	return((0.5-cos(maxalpha))  );
  //	return( 0.5*(0.5-cos(maxalpha)) + 0.5*(max_l/l-0.707)/4 );
  return( 0.5*(0.5-cos(maxalpha)) + 0.5*(max_l/l-0.707)/2 );
}


int ApplyRules ( const ARRAY<netrule*> & rules,
                 ARRAY<Point2d> & lpoints, ARRAY<ILINE> & llines,
                 ARRAY<Element> & elements,
                 ARRAY<INDEX> & dellines, int tolerance)
{
  int i, j, ri, nlok, npok, incnpok, refpi, locli;
  float maxerr = 0.3 * tolerance;
  char ok, found;
  netrule * rule;
  Vector oldu, newu;
  Vec2d ui;
  Point2d np;
  Vec2d linevec;
  int oldnp;
  ILINE loclin;
  double hf, elerr, minelerr;
  int noldlp, noldll;
  int loctestmode;

  static ARRAY<int> pused, lused;
  static ARRAY<int> pmap, lmap, pfixed;


  static ARRAY<Point2d> tempnewpoints;
  static ARRAY<ILINE> tempnewlines;
  static ARRAY<int> tempdellines;
  static ARRAY<Element> tempelements;



  elements.SetSize (0);
  dellines.SetSize (0);

  noldlp = lpoints.Size();
  noldll = llines.Size();

  pused.SetSize (lpoints.Size());
  lused.SetSize (llines.Size());

  loctestmode = (testmode && (tolerance % 50 == 0));
  if (loctestmode)
  {
    (*testout) << "tolerance = " << tolerance << endl;
    for (i = 1; i <= lpoints.Size(); i++)
      (*testout) << lpoints[i] << " ";
    (*testout) << endl;
    for (i = 1; i <= llines.Size(); i++)
      (*testout) << "(" << llines[i].I1() << "-" << llines[i].I2() << ")" << endl;
  }

  // check every rule

  found = 0;
  minelerr = 0.1 * tolerance * tolerance;

  for (ri = 1; ri <= rules.Size(); ri++)
  {
    rule = rules[ri];

    if (loctestmode)
      (*testout) << "Rule " << rule->Name() << endl;

    if (rule->GetQuality() > tolerance) continue;

    pmap.SetSize (rule->GetNP());
    lmap.SetSize (rule->GetNL());

    for (i = 1; i <= llines.Size(); i++)
      lused.Set (i, 0);
    for (i = 1; i <= lpoints.Size(); i++)
      pused.Set (i, 0);
    for (i = 1; i <= pmap.Size(); i++)
      pmap.Set(i, 0);
    for (i = 1; i <= lmap.Size(); i++)
      lmap.Set(i, 0);

    lused.Set (1, 1);
    lmap.Set (1, 1);

    for (j = 1; j <= 2; j++)
    {
      pmap[rule->GetPointNr (1, j)] = llines[1].I(j);
      pused[llines[1].I(j)]++;
    }


    nlok = 2;

    while (nlok >= 2)
    {

      if (nlok <= rule->GetNOldL())

      {
        ok = 0;
        while (!ok && lmap.Get(nlok) < llines.Size())
        {
          lmap.Elem(nlok)++;
          locli = lmap.Get(nlok);

          if (!lused.Get(locli))
          {
            ok = 1;

            loclin = llines.Get(locli);

            linevec.X() = lpoints.Get (loclin.I2()).X() -
                          lpoints.Get (loclin.I1()).X();
            linevec.Y() = lpoints.Get (loclin.I2()).Y() -
                          lpoints.Get (loclin.I1()).Y();

            if (rule->CalcLineError (nlok, linevec) > maxerr)
              ok = 0;

            for (j = 1; j <= 2 && ok; j++)
            {
              refpi = rule->GetPointNr (nlok, j);

              if (pmap.Get(refpi) != 0)
              {
                if (pmap.Get(refpi) != loclin.I(j))
                  ok = 0;
              }
              else
              {
                if (rule->CalcPointDist (refpi, lpoints.Get(loclin.I(j))) > maxerr
                    || pused[loclin.I(j)])
                  ok = 0;
              }
            }
          }
        }

        if (ok)
        {
          lused.Set (locli, 1);
          for (j = 1; j <= 2; j++)
          {
            pmap.Set(rule->GetPointNr (nlok, j), loclin.I(j));
            pused.Elem(loclin.I(j))++;
          }

          nlok++;
        }
        else
        {
          lmap[nlok] = 0;
          nlok--;

          lused.Set (lmap.Get(nlok), 0);
          for (j = 1; j <= 2; j++)
          {
            pused.Elem(llines.Get(lmap.Get(nlok)).I(j)) --;
            if (! pused.Get (llines.Get (lmap.Get (nlok)).I(j)))
              pmap.Set (rule->GetPointNr (nlok, j), 0);
          }
        }
      }

      else

      {

        // all lines are mapped !!

        // map also all points:

        npok = 1;
        incnpok = 1;

        pfixed.SetSize (pmap.Size());
        for (i = 1; i <= pmap.Size(); i++)
          pfixed[i] = (pmap[i] >= 1);

        while (npok >= 1)
        {

          if (npok <= rule->GetNOldP())

          {
            if (pfixed[npok])

            {
              if (incnpok)
                npok++;
              else
                npok--;
            }

            else

            {
              ok = 0;

              if (pmap[npok])
                pused[pmap[npok]]--;

              while (!ok && pmap[npok] < lpoints.Size())
              {
                ok = 1;

                pmap[npok]++;

                if (pused[pmap[npok]])
                {
                  ok = 0;
                }
                else
                {
                  if (rule->CalcPointDist (npok, lpoints[pmap[npok]]) > maxerr)
                    ok = 0;
                }
              }

              if (ok)
              {
                pused[pmap[npok]]++;
                npok++;
                incnpok = 1;
              }

              else

              {
                pmap[npok] = 0;
                npok--;
                incnpok = 0;
              }
            }
          }

          else

          {
            if (loctestmode)
              (*testout) << "lines and points mapped" << endl;

            oldu.SetLength (2 * rule->GetNOldP());

            for (i = 1; i <= rule->GetNOldP(); i++)
            {
              ui = lpoints.Get(pmap.Get(i)) - rule->GetPoint(i);
              oldu.Set (2*i-1, ui.X());
              oldu.Set (2*i  , ui.Y());
            }

            rule->GetOldUToFreeArea().Mult (oldu, newu);
            rule -> SetFreeZoneTransformation (newu);


            ok = 1;

            if (!rule->ConvexFreeZone())
            {
              ok = 0;
            }

            /*
                        rule->GetFreeArea (transfreearea);

                        for (i = 1; i <= transfreearea.Size(); i++)
                          {
                          transfreearea.Elem(i).X() += newu.Get(2*i - 1);
                          transfreearea.Elem(i).Y() += newu.Get(2*i);
                          }
             */

            // check freezone:

            for (i = 1; i <= lpoints.Size() && ok; i++)
            {
              if ( !pused.Get(i) &&
                   rule->IsInFreeZone (lpoints.Get(i)) )
                ok = 0;
            }


            for (i = 1; i <= llines.Size() && ok; i++)
            {
              if (!lused.Get(i) && rule->IsLineInFreeZone (
                    lpoints.Get(llines.Get(i).I1()),
                    lpoints.Get(llines.Get(i).I2())))

                ok = 0;
            }




            // check orientations

            for (i = 1; i <= rule->GetNOrientations() && ok; i++)
            {
              if (CW (lpoints[pmap[rule->GetOrientation(i).i1]],
                      lpoints[pmap[rule->GetOrientation(i).i2]],
                      lpoints[pmap[rule->GetOrientation(i).i3]]) ) ok = 0;
            }


            if (ok)
            {
              if (loctestmode)
                (*testout) << "rule ok" << endl;

              newu = rule->GetOldUToNewU() * oldu;

              // Setze neue Punkte:

              oldnp = rule->GetNOldP();

              for (i = oldnp + 1; i <= rule->GetNP(); i++)
              {
                np = rule->GetPoint(i);
                np.X() += newu (2 * (i-oldnp) - 1);
                np.Y() += newu (2 * (i-oldnp));

                pmap[i] = lpoints.Append (np);
              }

              // Setze neue Linien:

              for (i = rule->GetNOldL() + 1; i <= rule->GetNL(); i++)
              {
                llines.Append (ILINE (pmap[rule->GetPointNr (i, 1)],
                                      pmap[rule->GetPointNr (i, 2)]));
              }


              // delete old lines:

              for (i = 1; i <= rule->GetNDelL(); i++)
                dellines.Append (lmap[rule->GetDelLine(i)]);

              // insert new elements:

              for (i = 1; i <= rule->GetNE(); i++)
              {
                elements.Append (rule->GetElement(i));
                for (j = 1; j <= elements[i].NP(); j++)
                  elements[i].PNum(j) = pmap[elements[i].PNum(j)];
              }


              elerr = 0;
              for (i = 1; i <= elements.Size(); i++)
              {
                hf = CalcElementBadness (lpoints, elements[i]);
                if (hf > elerr) elerr = hf;
              }

              if (loctestmode)
                (*testout) << "error = " << elerr;


              if (elerr < minelerr)
              {

                if (testmode)
                {
                  (*testout) << "rule = " << rule->Name() << endl;
                  (*testout) << "class = " << tolerance << endl;
                  (*testout) << "lpoints: " << endl;
                  for (i = 1; i <= lpoints.Size(); i++)
                    (*testout) << lpoints[i] << endl;
                  (*testout) << "llines: " << endl;
                  for (i = 1; i <= llines.Size(); i++)
                    (*testout) << llines[i].I1() << " " << llines[i].I2() << endl;

                  (*testout) << "Freezone: ";
                  for (i = 1; i <= rule -> GetTransFreeZone().Size(); i++)
                    (*testout) << rule->GetTransFreeZone()[i] << endl;
                }



                minelerr = elerr;
                found = ri;

                tempnewpoints.SetSize (0);
                for (i = noldlp+1; i <= lpoints.Size(); i++)
                  tempnewpoints.Append (lpoints[i]);

                tempnewlines.SetSize (0);
                for (i = noldll+1; i <= llines.Size(); i++)
                  tempnewlines.Append (llines[i]);

                tempdellines.SetSize (0);
                for (i = 1; i <= dellines.Size(); i++)
                  tempdellines.Append (dellines[i]);

                tempelements.SetSize (0);
                for (i = 1; i <= elements.Size(); i++)
                  tempelements.Append (elements[i]);
              }

              lpoints.SetSize (noldlp);
              llines.SetSize (noldll);
              dellines.SetSize (0);
              elements.SetSize (0);
              ok = 0;

            }


            npok = rule->GetNOldP();
            incnpok = 0;
          }
        }

        nlok = rule->GetNOldL();

        lused.Set (lmap.Get(nlok), 0);

        for (j = 1; j <= 2; j++)
        {
          refpi = rule->GetPointNr (nlok, j);
          pused.Elem(pmap.Get(refpi))--;

          if (pused.Get(pmap.Get(refpi)) == 0)
          {
            pmap.Set(refpi, 0);
          }
        }
      }
    }
  }


  if (found)
  {
    for (i = 1; i <= tempnewpoints.Size(); i++)
      lpoints.Append (tempnewpoints[i]);
    for (i = 1; i <= tempnewlines.Size(); i++)
      llines.Append (tempnewlines[i]);
    for (i = 1; i <= tempdellines.Size(); i++)
      dellines.Append (tempdellines[i]);
    for (i = 1; i <= tempelements.Size(); i++)
      elements.Append (tempelements[i]);

    //    (*testout) << "minelerr = " << minelerr << endl;
  }


  return found;
}
