// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// abstract data type HASHTABLE

// by Joachim Sch”berl







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



























/*



   template <class HT, class T>

   class HASHTABLE

   {

   class colstruct

    {

    public:



    HT hash;

    T cont;



    colstruct () { }

    colstruct (const HT & ahash, const T & acont)

        : hash(ahash), cont(acont) { }

    };



   class linestruct

    {

    public:



    int size, maxsize;

    colstruct * col;

    };





   ARRAY<linestruct> data;





   public:



   HASHTABLE (int size);

   ~HASHTABLE ();



   void Set (const HT & ahash, const T & acont);

   const T & Get (const HT & ahash) const;

   int Used (const HT & ahash) const;



   int GetNBags () const { return data.Size(); }

   int GetBagSize (int bnr) const { return data.Get(bnr).size; }

   void GetData (int bnr, int colnr, HT & ahash, T & acont)

    {

    ahash = data.Get(bnr).col[colnr-1].hash;

    acont = data.Get(bnr).col[colnr-1].cont;

    }



   void Print (ostream & ost) const;

   };





   template<class HT, class T>

   HASHTABLE<HT, T> :: HASHTABLE (int size)

    : data(size)

   {

   int i;



   for (i = 1; i <= size; i++)

    {

    data.Elem(i).maxsize = 0;

    data.Elem(i).size = 0;

    data.Elem(i).col = NULL;

    }

   }



   template<class HT, class T>

   HASHTABLE<HT, T> :: ~HASHTABLE ()

   {

   int i;



   for (i = 1; i <= data.Size(); i++)

    {

    delete [] data[i].col;

    }

   }





   template<class HT, class T>

   void HASHTABLE<HT, T> :: Set (const HT & ahash, const T & acont)

   {

   int i;

   int hval = ahash.HashValue (data.Size());



   if (hval < 1 || hval > data.Size())

    {

    MyError ("Hashvalue out of Range");

    return;

    }



   linestruct & line = data.Elem (hval);



   for (i = 0; i < line.size; i++)

    {

    if (line.col[i].hash == ahash)

      {

      line.col[i].cont = acont;

      return;

      }

    }



   if (line.size == line.maxsize)

    {

    colstruct * ncol = new colstruct [line.maxsize+5];



    if (!ncol)

      {

      MyError ("Hashtable::Add: Out of memory");

      return;

      }



    memcpy (ncol, line.col, line.maxsize * sizeof (colstruct));

    delete line.col;

    line.col = ncol;

    line.maxsize += 5;

    }



   line.col[line.size].hash = ahash;

   line.col[line.size].cont = acont;

   line.size++;

   }



   template<class HT, class T>

   const T & HASHTABLE<HT, T> :: Get (const HT & ahash) const

   {

   int i;

   int hval = ahash.HashValue (data.Size());



   if (hval < 1 || hval > data.Size())

    {

    MyError ("Hashvalue out of Range");

    return * new T;

    }



   const linestruct & line = data.Get (hval);



   for (i = 0; i < line.size; i++)

    {

    if (line.col[i].hash == ahash)

      {

      return line.col[i].cont;

      }

    }



   MyError ("Hashtable :: Get: Element not found");





   return * new T;

   }





   template<class HT, class T>

   int HASHTABLE<HT, T> :: Used (const HT & ahash) const

   {

   int i;

   int hval = ahash.HashValue (data.Size());



   if (hval < 1 || hval > data.Size())

    {

    MyError ("Hashvalue out of Range");

    return 0;

    }



   const linestruct & line = data.Get (hval);



   for (i = 0; i < line.size; i++)

    {

    if (line.col[i].hash == ahash)

      {

      return 1;

      }

    }



   return 0;

   }





   template<class HT, class T>

   void HASHTABLE<HT, T> :: Print (ostream & ost) const

   {

   int i;



   ost << data.Size() << " bags" << endl;



   for (i = 1; i <= data.Size(); i++)

    ost << data[i].size << endl;

   }







   class HASHINDEX

   {

   INDEX i;



   public:

   HASHINDEX () { };

   HASHINDEX (INDEX ai) { i = ai; }



   operator INDEX () const { return i; }

   int HashValue (int bags) const { return i % bags + 1; }

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



   int HashValue (int bags) const { return (i1+i2) % bags + 1; }



   INDEX & I1 () { return i1; }

   INDEX & I2 () { return i2; }

   };



   class SORT_INDEX_2

   {

   INDEX i1, i2;



   public:

   SORT_INDEX_2 () { }

   SORT_INDEX_2 (INDEX ai1, INDEX ai2)

    {

    if (ai1 < ai2) { i1 = ai1; i2 = ai2; }

           else    { i1 = ai2; i2 = ai1; }

    }



   SORT_INDEX_2 (const SORT_INDEX_2 & si2)

    { i1 = si2.i1; i2 = si2.i2; }



   int operator== (const SORT_INDEX_2 & si2) const

    { return i1 == si2.i1 && i2 == si2.i2; }



   int HashValue (int bags) const { return i1 % bags + 1; }



   INDEX & I1 () { return i1; }

   INDEX & I2 () { return i2; }

   };



 */
