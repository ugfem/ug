// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
// #include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <string.h>

#include <myadt.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>

#include <meshing/ruler3.hh>
#include <meshing/meshing3.hh>


void LoadVMatrixLine (istream & ist, DenseMatrix & m, int line)
{
  char ch;
  int pnum;
  float f;

  ist >> ch;
  while (ch != '}')
  {
    ist.putback (ch);
    ist >> f;
    ist >> ch;
    ist >> pnum;

    if (ch == 'x' || ch == 'X')
      m(line, 3 * pnum - 2) = f;
    if (ch == 'y' || ch == 'Y')
      m(line, 3 * pnum - 1) = f;
    if (ch == 'z' || ch == 'Z')
      m(line, 3 * pnum    ) = f;

    ist >> ch;
    if (ch == ',')
      ist >> ch;
  }
}







void vnetrule :: LoadRule (istream & ist)
{
  char buf[256];
  char ch, ok;
  Point3d p;
  Element face;
  int i, j, i1, i2, i3, fs, ii, ii1, ii2, ii3;
  twoint edge;
  DenseMatrix tempoldutonewu(30, 20), tempoldutofreezone(30, 20);

  tempoldutonewu = 0;
  tempoldutofreezone = 0;

  noldp = 0;
  noldf = 0;

  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);
  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);

  name = new char[strlen (buf) + 1];
  strcpy (name, buf);
  //  cout << "Rule " << name << " found." << endl;

  do
  {
    ist >> buf;

    if (strcmp (buf, "quality") == 0)

    {
      ist >> quality;
    }

    else if (strcmp (buf, "flags") == 0)
    {
      ist >> ch;
      while (ch != ';')
      {
        flags.Append (ch);
        ist >> ch;
      }
    }

    else if (strcmp (buf, "mappoints") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> p.X();
        ist >> ch;    // ','
        ist >> p.Y();
        ist >> ch;    // ','
        ist >> p.Z();
        ist >> ch;    // ')'

        points.Append (p);
        noldp++;

        tolerances.SetSize (noldp);
        tolerances[noldp] = 1;

        ist >> ch;
        while (ch != ';')
        {
          if (ch == '{')
          {
            ist >> tolerances[noldp];
            ist >> ch;  // '}'
          }

          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }


    else if (strcmp (buf, "mapfaces") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        face.SetNP(3);
        ist >> face.PNum(1);
        ist >> ch;    // ','
        ist >> face.PNum(2);
        ist >> ch;    // ','
        ist >> face.PNum(3);
        ist >> ch;    // ')'

        faces.Append (face);
        noldf++;

        ist >> ch;
        while (ch != ';')
        {
          if (ch == 'd')
          {
            delfaces.Append (noldf);
            ist >> ch; // 'e'
            ist >> ch; // 'l'
          }

          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }

    else if (strcmp (buf, "mapedges") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> edge.i1;
        ist >> ch;    // ','
        ist >> edge.i2;
        ist >> ch;    // ')'

        edges.Append (edge);

        ist >> ch;
        while (ch != ';')
        {
          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }


    else if (strcmp (buf, "newpoints") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> p.X();
        ist >> ch;    // ','
        ist >> p.Y();
        ist >> ch;    // ','
        ist >> p.Z();
        ist >> ch;    // ')'

        points.Append (p);

        ist >> ch;
        while (ch != ';')
        {
          if (ch == '{')
          {
            LoadVMatrixLine (ist, tempoldutonewu,
                             3 * (points.Size()-noldp) - 2);

            ist >> ch; // '{'
            LoadVMatrixLine (ist, tempoldutonewu,
                             3 * (points.Size()-noldp) - 1);

            ist >> ch; // '{'
            LoadVMatrixLine (ist, tempoldutonewu,
                             3 * (points.Size()-noldp)    );
          }

          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }

    else if (strcmp (buf, "newfaces") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        face.SetNP (3);
        ist >> face.PNum(1);
        ist >> ch;    // ','
        ist >> face.PNum(2);
        ist >> ch;    // ','
        ist >> face.PNum(3);
        ist >> ch;    // ')'

        faces.Append (face);

        ist >> ch;
        while (ch != ';')
        {
          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }

    else if (strcmp (buf, "freezone") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> p.X();
        ist >> ch;    // ','
        ist >> p.Y();
        ist >> ch;    // ','
        ist >> p.Z();
        ist >> ch;    // ')'

        freezone.Append (p);

        ist >> ch;
        while (ch != ';')
        {
          if (ch == '{')
          {
            LoadVMatrixLine (ist, tempoldutofreezone,
                             3 * freezone.Size() - 2);

            ist >> ch; // '{'
            LoadVMatrixLine (ist, tempoldutofreezone,
                             3 * freezone.Size() - 1);

            ist >> ch; // '{'
            LoadVMatrixLine (ist, tempoldutofreezone,
                             3 * freezone.Size()    );
          }

          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }

    else if (strcmp (buf, "freeset") == 0)
    {
      freesets.Append (new ARRAY<int>);

      ist >> ch;

      while (ch != ';')
      {
        ist.putback (ch);
        ist >> i;
        freesets.Last()->Append(i);
        ist >> ch;
      }
    }

    else if (strcmp (buf, "elements") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        elements.Append (Element());

        elements.Last().SetNP(1);
        ist >> elements.Last().PNum(1);
        ist >> ch;    // ','

        if (ch == ',')
        {
          elements.Last().SetNP(2);
          ist >> elements.Last().PNum(2);
          ist >> ch;    // ','
        }
        if (ch == ',')
        {
          elements.Last().SetNP(3);
          ist >> elements.Last().PNum(3);
          ist >> ch;    // ','
        }
        if (ch == ',')
        {
          elements.Last().SetNP(4);
          ist >> elements.Last().PNum(4);
          ist >> ch;    // ','
        }

        orientations.Append (fourint());
        orientations.Last().i1 = elements.Last().PNum(1);
        orientations.Last().i2 = elements.Last().PNum(2);
        orientations.Last().i3 = elements.Last().PNum(3);
        orientations.Last().i4 = elements.Last().PNum(4);

        ist >> ch;
        while (ch != ';')
        {
          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }

    else if (strcmp (buf, "orientations") == 0)

    {
      ist >> ch;

      while (ch == '(')
      {
        //        fourint a = fourint();
        orientations.Append (fourint());

        ist >> orientations.Last().i1;
        ist >> ch;    // ','
        ist >> orientations.Last().i2;
        ist >> ch;    // ','
        ist >> orientations.Last().i3;
        ist >> ch;    // ','
        ist >> orientations.Last().i4;
        ist >> ch;    // ','


        ist >> ch;
        while (ch != ';')
        {
          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }


    else if (strcmp (buf, "endrule") != 0)
    {
      cout << "Unknown: " << buf << endl;
    }
  }
  while (!ist.eof() && strcmp (buf, "endrule") != 0);


  oldutonewu.SetSize (3 * (points.Size() - noldp), 3 * noldp);
  oldutofreezone.SetSize (3 * freezone.Size(), 3 * noldp);

  for (i = 1; i <= oldutonewu.Height(); i++)
    for (j = 1; j <= oldutonewu.Width(); j++)
      if (tempoldutonewu(i, j) &&
          j >= 4 && j != 5 && j != 6 && j != 9)
        oldutonewu(i, j) = tempoldutonewu(i, j);

  for (i = 1; i <= oldutofreezone.Height(); i++)
    for (j = 1; j <= oldutofreezone.Width(); j++)
      if (tempoldutofreezone(i, j) &&
          j >= 4 && j != 5 && j != 6 && j != 9)
        oldutofreezone(i, j) = tempoldutofreezone(i, j);


  for (i = 1; i <= elements.Size(); i++)
  {
    orientations.Append (fourint());
    orientations.Last().i1 = elements[i].PNum(1);
    orientations.Last().i2 = elements[i].PNum(2);
    orientations.Last().i3 = elements[i].PNum(3);
    orientations.Last().i4 = elements[i].PNum(4);
  }


  if (freesets.Size() == 0)
  {
    freesets.Append (new ARRAY<int>);
    for (i = 1; i <= freezone.Size(); i++)
      freesets[1]->Append(i);
  }


  //  testout << "Freezone: " << endl;
  for (fs = 1; fs <= freesets.Size(); fs++)
  {
    freefaces.Append (new ARRAY<threeint>);

    //    testout << "Freeset: " << endl;

    ARRAY<int> & freeset = *freesets[fs];
    ARRAY<threeint> & freesetfaces = *freefaces.Last();

    for (ii1 = 1; ii1 <= freeset.Size(); ii1++)
      for (ii2 = 1; ii2 <= freeset.Size(); ii2++)
        for (ii3 = 1; ii3 <= freeset.Size(); ii3++)
          if (ii1 < ii2 && ii1 < ii3 && ii2 != ii3)
          {
            i1 = freeset[ii1];
            i2 = freeset[ii2];
            i3 = freeset[ii3];

            Vec3d v1, v2, n;

            v1 = freezone[i3] - freezone[i1];
            v2 = freezone[i2] - freezone[i1];
            n = Cross (v1, v2);
            n /= n.Length();

            ok = 1;
            for (ii = 1; ii <= freeset.Size(); ii++)
            {
              i = freeset[ii];
              if (i != i1 && i != i2 && i != i3)
                if ( (freezone[i] - freezone[i1]) * n < 0 ) ok = 0;
            }

            if (ok)
            {
              freesetfaces.Append (threeint());
              freesetfaces.Last().i1 = i1;
              freesetfaces.Last().i2 = i2;
              freesetfaces.Last().i3 = i3;

              //              testout << "Freesetface: "
              //                      << i1 << "  " << i2 << "  " << i3 << endl;
            }
          }
  }

  for (fs = 1; fs <= freesets.Size(); fs++)
  {
    freefaceinequ.Append (new DenseMatrix (freefaces[fs]->Size(), 4));
  }


  {
    char ok;
    int minn;
    ARRAY<int> pnearness (noldp);

    for (i = 1; i <= pnearness.Size(); i++)
      pnearness[i] = 1000;

    for (j = 1; j <= 3; j++)
      pnearness[GetPointNr (1, j)] = 0;

    do
    {
      ok = 1;

      for (i = 1; i <= noldf; i++)
      {
        minn = 1000;
        for (j = 1; j <= 3; j++)
          minn = min (minn, pnearness[GetPointNr (i, j)]);

        for (j = 1; j <= 3; j++)
          if (pnearness[GetPointNr (i, j)] > minn+1)
          {
            ok = 0;
            pnearness[GetPointNr (i, j)] = minn+1;
          }
      }
    }
    while (!ok);

    fnearness.SetSize (noldf);

    for (i = 1; i <= noldf; i++)
    {
      fnearness[i] = 0;
      for (j = 1; j <= 3; j++)
        fnearness[i] += pnearness[GetPointNr (i, j)];
    }
  }
}








void Meshing3 :: LoadRules (char * filename)
{
  char buf[256];
  ifstream ist (filename);

  if (!ist.good())
  {
    cerr << "Rule description file " << filename << " not found" << endl;
    exit (1);
  }

  while (!ist.eof())
  {
    buf[0] = 0;
    ist >> buf;

    if (strcmp (buf, "rule") == 0)
    {
      vnetrule * rule = new vnetrule;
      rule -> LoadRule(ist);
      rules.Append (rule);
      if (!rule->TestOk())
      {
        cout << "Rule " << rules.Size() << " not ok" << endl;
        exit (1);
      }
    }
    else if (strcmp (buf, "tolfak") == 0)
    {
      ist >> tolfak;
    }
  }
}
