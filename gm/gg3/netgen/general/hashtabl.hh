// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_HASHTABL
#define FILE_HASHTABL

/**************************************************************************/
/* File:   hashtabl.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   Abstract data type HASHTABLE
 */



class BASE_INDEX_HASHTABLE
{
protected:
  TABLE<INDEX> hash;

public:
  BASE_INDEX_HASHTABLE (int size)
    : hash (size) { };

protected:
  int HashValue (const INDEX & ind) const;
  int Position (int bnr, const INDEX & ind) const;
};

template <class T>
class INDEX_HASHTABLE : private BASE_INDEX_HASHTABLE
{
  TABLE<T> cont;

public:
  INDEX_HASHTABLE (int size)
    : BASE_INDEX_HASHTABLE (size), cont(size)
  {
    ;
  }

  void Set (const INDEX & hash, const T & acont)
  {
    int bnr = HashValue (hash);
    int pos = Position (bnr, hash);
    if (pos)
      cont.Set (bnr, pos, acont);
    else
    {
      data.Add (bnr, hash);
      cont.Add (bnr, acont);
    }
  }

  const T & Get (const INDEX & ahash) const
  {
    int bnr = HashValue (hash);
    int pos = Position (bnr, hash);
    return cont.Get (bnr, pos);
  }

  int Used (const INDEX & ahash) const
  {
    return (Position (Hashvalue (hash), hash)) ? 1 : 0;
  }

  int GetNBags () const
  {
    return data.Size();
  }

  int GetBagSize (int bnr) const
  {
    return data.Get(bnr).size;
  }

  void GetData (int bnr, int colnr, INDEX & ahash, T & acont)
  {
    ahash = data.Get(bnr).col[colnr-1];
    acont = cont.Get(bnr).col[colnr-1];
  }
};













class INDEX_2
{
  INDEX i1, i2;

public:
  INDEX_2 () { }
  INDEX_2 (INDEX ai1, INDEX ai2)
  { i1 = ai1; i2 = ai2; }

  INDEX_2 (const INDEX_2 & in2)
  { i1 = in2.i1; i2 = in2.i2; }

  int operator== (const INDEX_2 & in2) const
  { return i1 == in2.i1 && i2 == in2.i2; }

  void Sort ();
  INDEX & I1 () { return i1; }
  INDEX & I2 () { return i2; }
  const INDEX & I1 () const { return i1; }
  const INDEX & I2 () const { return i2; }
};


class BASE_INDEX_2_HASHTABLE
{
protected:
  TABLE<INDEX_2> hash;

public:
  BASE_INDEX_2_HASHTABLE (int size)
    : hash (size) { };

protected:
  int HashValue (const INDEX_2 & ind) const
  {
    return (ind.I1() + ind.I2()) % hash.Size() + 1;
  }
  int Position (int bnr, const INDEX_2 & ind) const;
};


template <class T>
class INDEX_2_HASHTABLE : private BASE_INDEX_2_HASHTABLE
{
  TABLE<T> cont;

public:
  INDEX_2_HASHTABLE (int size)
    : BASE_INDEX_2_HASHTABLE (size), cont(size)
  {
    ;
  }

  void Set (const INDEX_2 & ahash, const T & acont)
  {
    int bnr = HashValue (ahash);
    int pos = Position (bnr, ahash);
    if (pos)
      cont.Set (bnr, pos, acont);
    else
    {
      hash.Add (bnr, ahash);
      cont.Add (bnr, acont);
    }
  }

  const T & Get (const INDEX_2 & ahash) const
  {
    int bnr = HashValue (ahash);
    int pos = Position (bnr, ahash);
    return cont.Get (bnr, pos);
  }

  int Used (const INDEX_2 & ahash) const
  {
    return (Position (HashValue (ahash), ahash)) ? 1 : 0;
  }

  int GetNBags () const
  {
    return cont.Size();
  }

  int GetBagSize (int bnr) const
  {
    return cont.EntrySize (bnr);
  }

  void GetData (int bnr, int colnr, INDEX_2 & ahash, T & acont)
  {
    ahash = hash.Get(bnr, colnr);
    acont = cont.Get(bnr, colnr);
  }
};


class INDEX_3
{
  INDEX i1, i2, i3;

public:
  INDEX_3 () { }
  INDEX_3 (INDEX ai1, INDEX ai2, INDEX ai3)
  { i1 = ai1; i2 = ai2; i3 = ai3; }

  INDEX_3 (const INDEX_3 & in2)
  { i1 = in2.i1; i2 = in2.i2; i3 = in2.i3; }

  void Sort ();

  int operator== (const INDEX_3 & in2) const
  { return i1 == in2.i1 && i2 == in2.i2 && i3 == in2.i3;}

  INDEX & I1 () { return i1; }
  INDEX & I2 () { return i2; }
  INDEX & I3 () { return i3; }
  const INDEX & I1 () const { return i1; }
  const INDEX & I2 () const { return i2; }
  const INDEX & I3 () const { return i3; }
};


class BASE_INDEX_3_HASHTABLE
{
protected:
  TABLE<INDEX_3> hash;

public:
  BASE_INDEX_3_HASHTABLE (int size)
    : hash (size) { };

protected:
  int HashValue (const INDEX_3 & ind) const
  {
    return (ind.I1() + ind.I2() + ind.I3()) % hash.Size() + 1;
  }
  int Position (int bnr, const INDEX_3 & ind) const;
};


template <class T>
class INDEX_3_HASHTABLE : private BASE_INDEX_3_HASHTABLE
{
  TABLE<T> cont;

public:
  INDEX_3_HASHTABLE (int size)
    : BASE_INDEX_3_HASHTABLE (size), cont(size)
  {
    ;
  }

  void Set (const INDEX_3 & ahash, const T & acont)
  {
    int bnr = HashValue (ahash);
    int pos = Position (bnr, ahash);
    if (pos)
      cont.Set (bnr, pos, acont);
    else
    {
      hash.Add (bnr, ahash);
      cont.Add (bnr, acont);
    }
  }

  const T & Get (const INDEX_3 & ahash) const
  {
    int bnr = HashValue (ahash);
    int pos = Position (bnr, ahash);
    return cont.Get (bnr, pos);
  }

  int Used (const INDEX_3 & ahash) const
  {
    return (Position (HashValue (ahash), ahash)) ? 1 : 0;
  }

  int GetNBags () const
  {
    return cont.Size();
  }

  int GetBagSize (int bnr) const
  {
    return cont.EntrySize (bnr);
  }

  void GetData (int bnr, int colnr, INDEX_3 & ahash, T & acont)
  {
    ahash = hash.Get(bnr, colnr);
    acont = cont.Get(bnr, colnr);
  }
};

#endif
