// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <fstream.h>
#include <string.h>
#include <math.h>

#include <template.hh>
#include <myadt.hh>
#include <linalg/linalg.hh>


//extern ofstream myerr;




SparseMatrix :: SparseMatrix () : BaseMatrix (), lins()
{};

SparseMatrix :: SparseMatrix (INDEX h, INDEX w) : BaseMatrix (h, w), lins (h)
{
  if (!w) w = h;

  for (INDEX i = 1; i <= h; i++)
  {
    lins[i].size = 0;
    lins[i].maxsize = 0;
    lins[i].col = NULL;
  }
}

SparseMatrix :: SparseMatrix (const SparseMatrix & m2) : BaseMatrix(), lins()
{
  (*this) = m2;
}


SparseMatrix :: ~SparseMatrix ()
{
  DeleteElements ();
}



void SparseMatrix :: SetSize (INDEX h, INDEX w)
{
  if (!w) w = h;
  DeleteElements ();

  if (height == h && width == w) return;

  height = h;
  width = w;
  lins.SetSize (height);

  if (lins.Size () == height)
  {
    for (INDEX i = 1; i <= h; i++)
    {
      lins[i].size = 0;
      lins[i].maxsize = 0;
      lins[i].col = NULL;
    }
  }
  else
  {
    height = width = 0;
    //    myerr << "SparseMatrix::SetSize: Out of memory" << endl;
  }
}


void SparseMatrix :: ChangeSize (INDEX h, INDEX w)
{
  INDEX i;

  if (!w) w = h;
  if (height == h && width == w) return;

  lins.SetSize (h);

  if (lins.Size () == h)
  {
    for (i = height+1; i <= h; i++)
    {
      lins[i].size = 0;
      lins[i].maxsize = 0;
      lins[i].col = NULL;
    }
    for (i = h+1; i <= height; i++)
    {
      if (lins[i].col)
      {
        DeleteColStruct (lins[i].col, lins[i].maxsize);

        lins[i].col = NULL;
        lins[i].size = 0;
        lins[i].maxsize = 0;
      }
    }

    height = h;
    width = w;
  }
  else
  {
    height = width = 0;
    //    myerr << "SparseMatrix::SetSize: Out of memory" << endl;
  }
}


void SparseMatrix :: SetSymmetric (int sym)
{
  INDEX i, j;
  int nr;

  if (sym != Symmetric())
  {
    BaseMatrix :: SetSymmetric (sym);

    if (!sym)
    {
      for (i = 1; i <= Height(); i++)
        for (nr = 1; nr <= ElementsInLine (i); nr++)
        {
          j = GetIndex (i, nr);
          if (j < i)
            Elem (j, i) = GetData (i, nr);
        }
    }
    else
    {
      DeleteElements();
    }
  }
}


void SparseMatrix :: DeleteElements ()
{
  for (INDEX i = 1; i <= lins.Size(); i++)
  {
    if (lins[i].col)
    {
      DeleteColStruct (lins[i].col, lins[i].maxsize);

      lins[i].col = NULL;
      lins[i].size = 0;
      lins[i].maxsize = 0;
    }
  }
}



double & SparseMatrix::operator() (INDEX i, INDEX j)
{
  if (i >= 1 && j >= 1 && i <= height && j <= width)
  {
    return Elem(i, j);
  }
  //  else myerr << "\nindex (" << i << "," << j << ") out of range (1.."
  //            << height << ",1.." << width << ")\n";
  return shit;
}

double SparseMatrix::operator() (INDEX i, INDEX j) const
{
  if (i >= 1 && j >= 1 && i <= height && j <= width)
  {
    return Get(i, j);
  }
  //  else myerr << "\nindex (" << i << "," << j << ") out of range (1.."
  //            << height << ",1.." << width << ")\n";
  return 0;
}

SparseMatrix & SparseMatrix :: operator= (const SparseMatrix & m2)
{
  INDEX i, j;

  SetSize (m2.Height(), m2.Width());
  SetSymmetric (m2.Symmetric());

  for (i = 1; i <= m2.Height(); i++)
    for (j = 1; j <= m2.ElementsInLine(i); j++)
      (*this).Elem(i, m2.GetIndex(i, j)) = m2.GetData(i, j);
  return *this;
}


SparseMatrix & SparseMatrix :: operator*= (double v)
{
  INDEX i, j;

  for (i = 1; i <= Height(); i++)
    for (j = 1; j <= ElementsInLine(i); j++)
      GetDataRef(i, j) *= v;
  return *this;
}



char SparseMatrix :: Used (INDEX i, INDEX j) const
{
  if (Symmetric() && j > i) swap (i, j);

  if (i < 1 || i > height) return 0;

  const colstruct * col = lins.Get(i).col;
  INDEX max = lins.Get(i).size;

  for (int k = 0; k < max; k++, col++)
    if (col->colnr == j) return 1;

  return 0;
}



double SparseMatrix :: Get(INDEX i, INDEX j) const
{
  if (Symmetric() && j > i) swap (i, j);

  const colstruct * col = lins.Get(i).col;
  INDEX max = lins.Get(i).size;

  for (INDEX k = 0; k < max; k++, col++)
    if (col->colnr == j) return col->data;

  return 0;
}


void SparseMatrix :: Set(INDEX i, INDEX j, double v)
{
  Elem (i, j) = v;
}



double & SparseMatrix :: Elem(INDEX i, INDEX j)
{
  if (Symmetric() && j > i) swap (i, j);

  linestruct * lin = &lins.Elem(i);
  colstruct * col, *ncol;

  int size = lin->size;
  int pos;

  if (size)
  {

    // bereits Elemente in der Liste

    if (j > lin->col[size-1].colnr)
    {

      // neues Element an letzter Position einfuegen

      if (lin->maxsize > size)
      {
        lin->size++;
        lin->col[size].colnr = j;
        return lin->col[size].data = 0;
      }

      if ( (ncol = NewColStruct(lin->maxsize+4) /* new colstruct[lin->maxsize+4] */) != NULL)
      {
        memcpy (ncol, lin->col, sizeof(colstruct) * size);

        DeleteColStruct (lin->col, lin->maxsize);

        lin->maxsize += 4;
        lin->col = ncol;
        lin->size++;
        ncol[size].colnr = j;
        return ncol[size].data = 0;
      }
      else
      {
        //      myerr << "SparseMatrix::Elem: Out of memory 1" << endl;
        return shit;
      }

    }
    else
    {

      for (col = lin->col; col->colnr < j; col++) ;
      // Zeilenliste durchsuchen

      if (col->colnr == j) return col->data;
      // Element exisitiert bereits

      if (lin->maxsize > size)
      {
        memmove (col+1, col, size_t((char*)&lin->col[size] - (char*)col));

        lin->size++;
        col->colnr = j;
        return col->data = 0;
      }

      pos = size_t (col - lin->col);

      if ( (ncol = NewColStruct(lin->maxsize+4) ) != NULL)
      {

        if (pos) memcpy (ncol, lin->col, sizeof(colstruct) * pos);
        memcpy (ncol+(pos+1), col, sizeof(colstruct) * (size-pos));

        DeleteColStruct (lin->col, lin->maxsize);
        //        delete lin->col;

        lin->maxsize += 4;
        lin->col = ncol;
        lin->size++;
        ncol[pos].colnr = j;
        return ncol[pos].data = 0;
      }
      else
      {
        //       myerr << "SparseMatrix::Elem: Out of memory 2" << endl;
        return shit;
      }
    }
  }
  else
  {
    // kein Element in der Liste

    if (lin->maxsize)
    {
      // Liste bereits angelegt

      lin->size = 1;
      lin->col->colnr = j;
      return lin->col->data = 0;
    }
    else
    {
      if ( (lin->col = NewColStruct(6) ) != NULL )
      {
        lin->maxsize = 6;
        lin->size = 1;
        lin->col->colnr = j;
        return lin->col->data = 0;
      }
      else
      {
        //        myerr << "SparseMatrix::Elem: Out of memory 3" << endl;
        return shit;
      }
    }
  }
}



void SparseMatrix :: Delete (INDEX i, int nr)
{
  colstruct * col = lins[i].col;

  nr--;
  while (nr < lins[i].size-1)
  {
    col[nr].data = col[nr+1].data;
    col[nr].colnr = col[nr+1].colnr;
    nr++;
  }
  lins[i].size--;
}

void SparseMatrix :: DeleteElem (INDEX i, INDEX j)
{
  int nr;
  int dec = 0;

  if (Symmetric() && j > i) swap (i, j);

  colstruct * col = lins[i].col;

  for (nr = 0; nr < lins[i].size; nr++)
  {
    if (dec)
    {
      col[nr-1].data = col[nr].data;
      col[nr-1].colnr = col[nr].colnr;
    }
    if (col[nr].colnr == j) dec = 1;
  }
  if (dec)
    lins[i].size--;
}


void SparseMatrix :: SetLineAllocSize (INDEX i, int j)
{
  colstruct * ncol;


  if (lins[i].maxsize < j)
  {
    if ( (ncol = NewColStruct(j)) != NULL)
    {
      memcpy (ncol, lins[i].col, sizeof(colstruct) * lins[i].size);

      if (lins[i].maxsize)
        DeleteColStruct (lins[i].col, lins[i].maxsize);

      lins[i].maxsize = j;
      lins[i].col = ncol;
    }
    //    else
    //      myerr << "SPARSE_MATIRX :: SetLineAllocSize: Out of Memory" << endl;
  }
}


ostream & SparseMatrix :: Print (ostream & s) const
{
  INDEX i, j;

  if (Symmetric()) s << "Lower half of matrix:" << endl;

  for (i = 1; i <= Height(); i++)
  {
    s << "Line " << i << " ";
    for (j = 1 ; j <= ElementsInLine (i); j++)
      s << "(" << GetIndex (i, j) << ": " << GetData (i, j) << ") ";
    s << endl;
  }
  return s;
}












void SparseMatrix :: Mult (const BaseVector & bv, BaseVector & bprod) const
{
  double sum, vi;
  INDEX i, j, n;
  colstruct * col;
  const linestruct * lin;

  const Vector & v = bv.CastToVector();
  Vector & prod = bprod.CastToVector();

  prod.SetLength (Height());

  if (prod.Length() != Height() || v.Length() != Width())
  {
    //   myerr << "SparseMatrix::Mult: Dimensions don't fit" << endl;
    return;
  }

  if (!Symmetric())
  {
    n = Height();
    for (i = 1; i <= n; i++)
    {
      lin = &lins.Get(i);
      sum = 0;
      col = lin->col;

      for (j = lin->size; j > 0; j--, col++)
        sum += col->data * v.Get(col->colnr);

      prod.Set (i, sum);
    }
  }
  else
  {
    prod = 0;
    n = Height();
    for (i = 1; i <= n; i++)
    {
      lin = &lins.Get(i);
      sum = 0;
      col = lin->col;
      vi = v.Get(i);

      for (j = lin->size; j > 0; j--, col++)
      {
        sum += col->data * v.Get(col->colnr);
        if (col->colnr != i)
          prod.Elem(col->colnr) += col->data * vi;
      }

      prod.Elem(i) += sum;
    }
  }
}

void SparseMatrix :: MultTrans (const BaseVector & bv, BaseVector & bprod) const
{
  INDEX i, j, n;
  colstruct * col;
  const linestruct * lin;
  const Vector & v = bv.CastToVector();
  Vector & prod = bprod.CastToVector();


  prod.SetLength (Width());

  if (prod.Length() != Width() || v.Length() != Height())
  {
    //   myerr << "SparseMatrix::Mult: Dimensions don't fit" << endl;
    return;
  }

  if (!Symmetric())
  {
    n = Height();

    for (i = 1; i <= Width(); i++)
      prod.Set (i, 0);

    for (i = 1; i <= n; i++)
    {
      lin = &lins[i];
      col = lin->col;

      for (j = lin->size; j > 0; j--, col++)
        prod.Elem(col->colnr) += col->data * v.Get(i);
    }
  }
  else
  {
    Mult (v, prod);
  }
}




void SparseMatrix :: Residuum (const BaseVector & bx, const BaseVector & bb,
                               BaseVector & bres) const
{
  double sum, xi;
  INDEX i, j, n;
  colstruct * col;
  const linestruct * lin;
  const Vector & x = bx.CastToVector();
  const Vector & b = bb.CastToVector();
  Vector & res = bres.CastToVector();

  res.SetLength (b.Length());

  if (res.Length() != b.Length() || b.Length() != Height() ||
      x.Length() != Width())
  {
    //    myerr << "SparseMatrix::Residuum: Dimensions don't fit" << endl;
    return;
  }

  n = Height();
  if (!Symmetric())
  {
    for (i = 1; i <= Height(); i++)
    {
      lin = &lins.Get(i);
      sum = b.Get(i);
      col = lin->col;

      for (j = lin->size; j > 0; j--, col++)
        sum -= col->data * x.Get(col->colnr);

      res.Set (i, sum);
    }
  }
  else
  {
    res = b;
    for (i = 1; i <= n; i++)
    {
      lin = &lins.Get(i);
      sum = 0;
      col = lin->col;
      xi = x.Get(i);

      for (j = lin->size; j > 0; j--, col++)
      {
        sum -= col->data * x.Get(col->colnr);
        if (col->colnr != i)
          res.Elem(col->colnr) -= col->data * xi;
      }

      res.Elem(i) += sum;
    }
  }
}


void SparseMatrix :: ResiduumTrans (const BaseVector & bx, const BaseVector & bb,
                                    BaseVector & bres) const
{
  INDEX i, j, n;
  colstruct * col;
  const linestruct * lin;
  const Vector & x = bx.CastToVector();
  const Vector & b = bb.CastToVector();
  Vector & res = bres.CastToVector();


  res.SetLength (Width());

  if (res.Length() != b.Length() || b.Length() != Width() ||
      x.Length() != Height())
  {
    //    myerr << "SparseMatrix::ResiduumTrans: Dimensions don't fit" << endl;
    return;
  }

  if (!Symmetric())
  {
    n = Height();

    res = b;

    for (i = 1; i <= n; i++)
    {
      lin = &lins.Get(i);
      col = lin->col;

      for (j = lin->size; j > 0; j--, col++)
        res.Elem(col->colnr) -= col->data * x.Get(i);
    }
  }
  else
  {
    Residuum (bx, bb, bres);
  }
}








SparseMatrix operator* (const SparseMatrix & m1,
                        const SparseMatrix & m2)
{
  SparseMatrix m(m1.Height(), m2.Width());
  INDEX i, j, k, nr;

  if (!m1.Symmetric() && !m2.Symmetric())
  {
    for (i = 1; i <= m1.Height(); i++)
    {
      for (j = 1; j <= m1.ElementsInLine(i); j++)
      {
        nr = m1.GetIndex (i, j);
        for (k = 1; k <= m2.ElementsInLine(nr); k++)
        {
          m(i, m2.GetIndex(nr, k)) += m1.GetData(i, j) * m2.GetData(nr, k);
        }
      }
    }
  }
  else
  {
    //   myerr << "SparseMatrix :: operator* not implemented for symmetric case" << endl;
  }
  return m;
}


SparseMatrix & SparseMatrix :: operator+= (const SparseMatrix & m2)
{
  INDEX i, k;
  int j;

  if (Symmetric() == m2.Symmetric())
  {
    for (i = 1; i <= Height(); i++)
      for (j = 1; j <= m2.ElementsInLine (i); j++)
      {
        k = m2.GetIndex (i, j);
        Elem(i, k) += m2.GetData (i, j);
      }
  }
  else
  {
    //    myerr << "SparseMatrix :: operator+= not implemented for different cases" << endl;
  }
  return *this;
}



SparseMatrix & SparseMatrix :: operator*= (const SparseMatrix & m2)
{
  INDEX i, k;
  int j, l;
  colstruct * cs;
  int ms, s;

  if (!Symmetric() && !m2.Symmetric())
  {
    for (i = 1; i <= Height(); i++)
    {
      cs = lins[i].col;
      s = lins[i].size;
      ms = lins[i].maxsize;

      lins[i].col = NULL;
      lins[i].size = 0;
      lins[i].maxsize = 0;


      for (j = 0; j < s; j++)
      {
        k = cs[j].colnr;

        for (l = 1; l <= m2.ElementsInLine (k); l++)
          Elem(i, m2.GetIndex(k, l)) += cs[j].data * m2.GetData(k, l);
      }

      DeleteColStruct (cs, ms);
    }
  }
  else
  {
    //   myerr << "SparseMatrix :: operator*= not implemented for Symmetric matrices" << endl;
  }

  return *this;
}


void SparseMatrix :: Solve (const Vector & v, Vector & sol) const
{
  SparseMatrix temp = *this;
  INDEX i, j, nr, k;
  double q;

  sol = v;


  if (!Symmetric())
  {
    for (i = 1; i <= Height(); i++)
    {
      if (temp.ElementsInLine(i) < 1 || temp.GetIndex(i, 1) != i)
      {
        //        myerr << "i = " << i << endl;
        //        myerr << "Solve: Matrix singular" << endl;
        char ch;
        cin >> ch;
      }
      for (j = 2; j <= temp.ElementsInLine(i); j++)
      {
        nr = temp.GetIndex(i, j);
        if (temp.GetIndex(nr, 1) != i)
        {
          //          myerr << temp << endl;
          //          myerr << "i = " << i << "j = " << j << "nr = " << nr << endl;
          //          myerr << "Solve: Graph not symmetrix" << endl;
          char ch;
          cin >> ch;
        }

        q = temp.GetData (nr, 1) / temp.GetData(i, 1);
        temp.Delete (nr, 1);

        for (k = 2; k <= temp.ElementsInLine (i); k++)
          temp.Elem(nr, temp.GetIndex(i, k)) -= q * temp.GetData(i, k);

        sol(nr) -= q * sol(i);
      }
    }

    for (i = temp.Height(); i >= 1; i--)
    {
      for (j = 2; j <= temp.ElementsInLine(i); j++)
      {
        sol(i) -= temp.GetData(i, j) * sol(temp.GetIndex(i, j));
      }
      sol(i) /= temp.GetData(i, 1);
    }
  }
  //  else
  //    myerr << "SparseMatrix :: Solve not implemented for symmetic case" << endl;
}





void Transpose (const SparseMatrix & m1, SparseMatrix & m2)
{
  INDEX i, j;

  m2.SetSize(m1.Width(), m1.Height());
  m2.SetSymmetric(m1.Symmetric());

  if (!m1.Symmetric())
  {
    for (i = 1; i <= m1.Height(); i++)
      for (j = 1; j <= m1.ElementsInLine(i); j++)
        m2(m1.GetIndex(i, j), i) = m1.GetData(i, j);
  }
  else
    m2 = m1;
}



BaseMatrix * SparseMatrix :: Copy () const
{
  return new SparseMatrix (*this);
}





void SparseMatrix :: AddElementMatrix (const ARRAY<INDEX> & pnum,
                                       const BaseMatrix & elemmat)
{
  int i1, i2;

  if (Symmetric())
  {
    for (i1 = 1; i1 <= pnum.Size(); i1++)
      for (i2 = 1; i2 <= i1; i2++)
        Elem(pnum[i1], pnum[i2]) += elemmat(i1, i2);
  }
  else
  {
    for (i1 = 1; i1 <= pnum.Size(); i1++)
      for (i2 = 1; i2 <= pnum.Size(); i2++)
        Elem(pnum[i1], pnum[i2]) += elemmat(i1, i2);
  }

}



void SparseMatrix :: AddRowMatrix (INDEX row, const SparseMatrix & m2)
{
  int i1, i2, i;
  colstruct * cs1, * cs2, * ncs;
  int s1, s2, s;

  s1 = lins[row].size;
  s2 = m2.lins[1].size;
  cs1 = lins[row].col;
  cs2 = m2.lins[1].col;

  i1 = 0;
  i2 = 0;
  i = 0;


  if (Symmetric())
  {
    while (s2 && cs2[s2-1].colnr > row) s2--;
  }


  while (i1 < s1 && i2 < s2)
  {
    if (cs1[i1].colnr < cs2[i2].colnr) i1++;
    else if (cs1[i1].colnr > cs2[i2].colnr) i2++;
    else { i1++; i2++; }
    i++;
  }

  i += (s1 - i1);
  i += (s2 - i2);

  s = i;

  if (s > s1)
  {
    ncs = NewColStruct (s);

    i1 = 0;
    i2 = 0;
    i = 0;

    while (i1 < s1 && i2 < s2)
    {
      if (cs1[i1].colnr < cs2[i2].colnr)
      {
        ncs[i].colnr = cs1[i1].colnr;
        ncs[i].data = cs1[i1].data;
        i1++;
      }
      else if (cs1[i1].colnr > cs2[i2].colnr)
      {
        ncs[i].colnr = cs2[i2].colnr;
        ncs[i].data = cs2[i2].data;
        i2++;
      }
      else
      {
        ncs[i].colnr = cs1[i1].colnr;
        ncs[i].data = cs1[i1].data + cs2[i2].data;
        i1++;
        i2++;
      }
      i++;
    }

    while (i1 < s1)
    {
      ncs[i].colnr = cs1[i1].colnr;
      ncs[i].data = cs1[i1].data;
      i1++;
      i++;
    }

    while (i2 < s2)
    {
      ncs[i].colnr = cs2[i2].colnr;
      ncs[i].data = cs2[i2].data;
      i2++;
      i++;
    }

    if (lins[row].maxsize)
      DeleteColStruct (cs1, lins[row].maxsize);

    lins[row].col = ncs;
    lins[row].size = s;
    lins[row].maxsize = s;
  }


  else
  {
    i1 = 0;
    i2 = 0;

    while (i2 < s2)
    {
      if (cs1[i1].colnr == cs2[i2].colnr)
      {
        cs1[i1].data += cs2[i2].data;
        i2++;
      }
      i1++;
    }
  }
}


double SparseMatrix :: RowTimesVector (INDEX i, const Vector & v) const
{
  const linestruct * lin;
  const colstruct * col;
  double sum;
  int j;

  if (Width() != v.Length())
  {
    cerr << "SparseMatrix::RowTimesVector: Size doesn't fit" << endl;
    return 0;
  }

  lin = &lins.Get(i);
  sum = 0;
  col = lin->col;

  for (j = lin->size; j > 0; j--, col++)
    sum += col->data * v.Get(col->colnr);

  return sum;
}

void SparseMatrix :: AddRowToVector (INDEX i, double s, Vector & v) const
{
  const linestruct * lin;
  const colstruct * col;
  int j;

  if (Width() != v.Length())
  {
    cerr << "SparseMatrix::RowTimesVector: Size doesn't fit" << endl;
    return;
  }

  lin = &lins.Get(i);
  col = lin->col;

  for (j = lin->size; j > 0; j--, col++)
    v.Elem(col->colnr) += s * col->data;
}


static ARRAY<void*> poolarray;
static ARRAY<int> poolsizes;


SparseMatrix::colstruct * SparseMatrix :: NewColStruct (int s)
{
  //  return new colstruct[s];


  int i, j;
  void * p;
  colstruct * cs;

  if (s > 20) return new colstruct[s];

  for (i = 1; i <= poolsizes.Size(); i++)
    if (poolsizes.Get(i) == s) break;

  if (i > poolsizes.Size())
  {
    poolsizes.Append(s);
    poolarray.Append((void*)NULL);
    i = poolsizes.Size();
  }

  p = poolarray.Get(i);
  if (!p)
  {
    poolarray.Elem(i) = p = cs = new colstruct[10 * s];
    for (j = 0; j < 9; j++)
      *(void**)(void*)(&cs[s * j]) = &cs[s * (j+1)];
    *(void**)(void*)(&cs[s * 9]) = NULL;
  }

  poolarray.Elem(i) = *(void**)p;
  return (colstruct*)p;
}

void SparseMatrix :: DeleteColStruct (colstruct * cs, int s)
{
  //  delete cs;
  //  return;

  int i;

  if (s > 20)
  {
    delete cs;
    return;
  }

  for (i = 1; i <= poolsizes.Size(); i++)
    if (poolsizes.Get(i) == s) break;


  if (i > poolsizes.Size())
  {
    //    myerr << "SparseMatrix :: DeleteColStruct: Size not found" << endl;
    return;
  }


  *(void**)(void*)cs = poolarray.Get(i);
  poolarray.Elem(i) = cs;
}
