// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
class BITARRAY

{

  INDEX size;

  unsigned char * data;



public:

  BITARRAY ();

  BITARRAY (INDEX asize);

  ~BITARRAY ();



  void SetSize (INDEX asize);

  inline INDEX Size () const;



  void Set ();

  inline void Set (INDEX i);

  void Clear ();

  inline void Clear (INDEX i);

  inline int Test (INDEX i) const;



  //  private:

  inline unsigned char Mask (INDEX i) const;

  inline INDEX Addr (INDEX i) const;

};





inline

INDEX BITARRAY :: Size () const

{

  return size;

}



inline

unsigned char BITARRAY :: Mask (INDEX i) const

{

  return char(1) << (i % CHAR_BIT);

}



inline

INDEX BITARRAY :: Addr (INDEX i) const

{

  return (i / CHAR_BIT);

}



inline

void BITARRAY :: Set (INDEX i)

{

  data[Addr(i)] |= Mask(i);

}



inline

void BITARRAY :: Clear (INDEX i)

{

  data[Addr(i)] &= ~Mask(i);

}



inline

int BITARRAY :: Test (INDEX i) const

{

  return (data[Addr(i)] & Mask(i)) ? 1 : 0;

}
