// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
// #include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <string.h>

#include <template.hh>
#include <array.hh>

#include <linalg/linalg.hh>
#include <geom/geom2d.hh>
#include <geom/geom3d.hh>

#include <meshing/global.hh>

#include <meshing/ruler2.hh>
#include <meshing/meshing2.hh>



void LoadMatrixLine (istream & ist, DenseMatrix & m, int line)
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
      m(line, 2 * pnum - 1) = f;
    if (ch == 'y' || ch == 'Y')
      m(line, 2 * pnum) = f;

    ist >> ch;
    if (ch == ',')
      ist >> ch;
  }
}


void netrule :: LoadRule (istream & ist)
{
  char buf[256];
  char ch;
  Point2d p;
  ILINE lin;
  int i, j;
  DenseMatrix tempoldutonewu(20, 20), tempoldutofreearea(20, 20);

  tempoldutonewu = 0;
  tempoldutofreearea = 0;

  noldp = 0;
  noldl = 0;

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

    else if (strcmp (buf, "mappoints") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> p.X();
        ist >> ch;    // ','
        ist >> p.Y();
        ist >> ch;    // ')'

        points.Append (p);
        noldp++;

        tolerances.SetSize (noldp);
        tolerances[noldp].f1 = 0;
        tolerances[noldp].f2 = 0;
        tolerances[noldp].f3 = 0;

        ist >> ch;
        while (ch != ';')
        {
          if (ch == '{')
          {
            ist >> tolerances[noldp].f1;
            ist >> ch;  // ','
            ist >> tolerances[noldp].f2;
            ist >> ch;  // ','
            ist >> tolerances[noldp].f3;
            ist >> ch;  // '}'
          }
          else if (ch == 'd')
          {
            //            delpoints.Append (noldp);
            ist >> ch; // 'e'
            ist >> ch; // 'l'
          }

          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }


    else if (strcmp (buf, "maplines") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> lin.I1();
        ist >> ch;    // ','
        ist >> lin.I2();
        ist >> ch;    // ')'

        lines.Append (lin);
        linevecs.Append (points[lin.I2()] - points[lin.I1()]);
        noldl++;
        linetolerances.SetSize (noldl);
        linetolerances[noldl].f1 = 0;
        linetolerances[noldl].f2 = 0;
        linetolerances[noldl].f3 = 0;

        ist >> ch;
        while (ch != ';')
        {
          if (ch == '{')
          {
            ist >> linetolerances[noldl].f1;
            ist >> ch;  // ','
            ist >> linetolerances[noldl].f2;
            ist >> ch;  // ','
            ist >> linetolerances[noldl].f3;
            ist >> ch;  // '}'
          }
          else if (ch == 'd')
          {
            dellines.Append (noldl);
            ist >> ch; // 'e'
            ist >> ch; // 'l'
          }

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
        ist >> ch;    // ')'

        points.Append (p);

        ist >> ch;
        while (ch != ';')
        {
          if (ch == '{')
          {
            LoadMatrixLine (ist, tempoldutonewu,
                            2 * (points.Size()-noldp) - 1);

            ist >> ch; // '{'
            LoadMatrixLine (ist, tempoldutonewu,
                            2 * (points.Size()-noldp));
          }

          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }

    else if (strcmp (buf, "newlines") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> lin.I1();
        ist >> ch;    // ','
        ist >> lin.I2();
        ist >> ch;    // ')'

        lines.Append (lin);
        linevecs.Append (points[lin.I2()] - points[lin.I1()]);

        ist >> ch;
        while (ch != ';')
        {
          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
    }

    else if (strcmp (buf, "freearea") == 0)
    {
      ist >> ch;

      while (ch == '(')
      {
        ist >> p.X();
        ist >> ch;    // ','
        ist >> p.Y();
        ist >> ch;    // ')'

        freezone.Append (p);

        ist >> ch;
        while (ch != ';')
        {
          if (ch == '{')
          {
            LoadMatrixLine (ist, tempoldutofreearea,
                            2 * freezone.Size() - 1);

            ist >> ch; // '{'
            LoadMatrixLine (ist, tempoldutofreearea,
                            2 * freezone.Size());
          }

          ist >> ch;
        }

        ist >> ch;
      }

      ist.putback (ch);
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
        //        threeint a = threeint();
        orientations.Append (threeint());

        ist >> orientations.Last().i1;
        ist >> ch;    // ','
        ist >> orientations.Last().i2;
        ist >> ch;    // ','
        ist >> orientations.Last().i3;
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

  oldutonewu.SetSize (2 * (points.Size() - noldp), 2 * noldp);
  oldutofreearea.SetSize (2 * freezone.Size(), 2 * noldp);

  for (i = 1; i <= oldutonewu.Height(); i++)
    for (j = 1; j <= oldutonewu.Width(); j++)
      oldutonewu(i, j) = tempoldutonewu(i, j);

  for (i = 1; i <= oldutofreearea.Height(); i++)
    for (j = 1; j <= oldutofreearea.Width(); j++)
      oldutofreearea(i, j) = tempoldutofreearea(i, j);

  freesetinequ.SetSize (freezone.Size(), 3);
}



void Meshing2 :: LoadRules (char * filename)
{
  char buf[256];
  ifstream ist (filename);

  while (!ist.eof())
  {
    buf[0] = 0;
    ist >> buf;

    if (strcmp (buf, "rule") == 0)
    {
      netrule * rule = new netrule;
      rule -> LoadRule(ist);
      rules.Append (rule);
    }
  }
}
