// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* Geometric Algorithms */



#define EPSGEOM 1E-5



extern void MyError (char * ch);

class ostream;



class POINT2D;

class VEC2D;

class LINE2D;







class POINT2D

{

  friend VEC2D;



protected:

  double px, py;



public:

  POINT2D() { /* px = py = 0; */ }

  POINT2D(double ax, double ay) { px = ax; py = ay; }

  POINT2D(const POINT2D & p2) { px = p2.px; py = p2.py; }



  POINT2D & operator= (const POINT2D & p2)

  { px = p2.px; py = p2.py; return *this; }



  double & X() { return px; }

  double & Y() { return py; }

  double X() const { return px; }

  double Y() const { return py; }



  friend inline VEC2D operator- (const POINT2D & p1, const POINT2D & p2);

  friend inline POINT2D operator- (const POINT2D & p1, const VEC2D & v);

  friend inline POINT2D operator+ (const POINT2D & p1, const VEC2D & v);



  friend inline POINT2D Center (const POINT2D & p1, const POINT2D & p2);



  friend double Dist (const POINT2D & p1, const POINT2D & p2)

  { return sqrt ( sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ); }

  friend double Dist2 (const POINT2D & p1, const POINT2D & p2)

  { return sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ; }



  friend int CW (const POINT2D & p1, const POINT2D & p2, const POINT2D & p3);

  friend int CCW (const POINT2D & p1, const POINT2D & p2, const POINT2D & p3);



  friend inline void PpSmV (const POINT2D & p1, double s, const VEC2D & v, POINT2D & p2);

  friend inline void PmP (const POINT2D & p1, const POINT2D & p2, VEC2D & v);



  friend ostream & operator<<(ostream  & s, const POINT2D & p);

};







class VEC2D

{

protected:

  double vx, vy;



public:

  VEC2D() { /* vx = vy = 0; */ }

  VEC2D(double ax, double ay)

  { vx = ax; vy = ay; }

  VEC2D(const VEC2D & v2)

  { vx = v2.vx; vy = v2.vy; }

  VEC2D(const POINT2D & p1, const POINT2D & p2)

  { vx = p2.px - p1.px; vy = p2.py - p1.py; }



  VEC2D & operator= (const VEC2D & p2)

  { vx = p2.vx; vy = p2.vy; return *this; }



  double & X() { return vx; }

  double & Y() { return vy; }

  double X() const { return vx; }

  double Y() const { return vy; }



  double Length() const { return sqrt (vx * vx + vy * vy); }

  double Length2() const { return vx * vx + vy * vy; }



  inline VEC2D & operator+= (const VEC2D & v2);

  inline VEC2D & operator-= (const VEC2D & v2);

  inline VEC2D & operator*= (double s);

  inline VEC2D & operator/= (double s);



  friend inline VEC2D operator- (const POINT2D & p1, const POINT2D & p2);

  friend inline POINT2D operator- (const POINT2D & p1, const VEC2D & v);

  friend inline POINT2D operator+ (const POINT2D & p1, const VEC2D & v);

  friend inline VEC2D operator- (const VEC2D & p1, const VEC2D & v);

  friend inline VEC2D operator+ (const VEC2D & p1, const VEC2D & v);

  friend inline VEC2D operator* (double scal, const VEC2D & v);



  friend double operator* (const VEC2D & v1, const VEC2D & v2)

  { return v1.X() * v2.X() + v1.Y() * v2.Y(); }





  friend double Cross (const VEC2D & v1, const VEC2D & v2)

  {
    return double(v1.X()) * double(v2.Y()) -

           double(v1.Y()) * double(v2.X());
  }



  friend inline void PpSmV (const POINT2D & p1, double s, const VEC2D & v, POINT2D & p2);

  friend inline void PmP (const POINT2D & p1, const POINT2D & p2, VEC2D & v);





  friend double Angle (const VEC2D & v);

  friend double FastAngle (const VEC2D & v);

  friend double Angle (const VEC2D & v1, const VEC2D & v2);

  friend double FastAngle (const VEC2D & v1, const VEC2D & v2);



  friend ostream & operator<<(ostream  & s, const VEC2D & v);

};







class LINE2D

{

protected:

  POINT2D p1, p2;



public:

  LINE2D() : p1(), p2() { };

  LINE2D(const POINT2D & ap1, const POINT2D & ap2)

  { p1 = ap1; p2 = ap2; }



  LINE2D & operator= (const LINE2D & l2)

  { p1 = l2.p1; p2 = l2.p2; return *this;}



  POINT2D & P1() { return p1; }

  POINT2D & P2() { return p2; }

  const POINT2D & P1() const { return p1; }

  const POINT2D & P2() const { return p2; }



  double XMax() const { return max (p1.X(), p2.X()); }

  double YMax() const { return max (p1.Y(), p2.Y()); }

  double XMin() const { return min (p1.X(), p2.X()); }

  double YMin() const { return min (p1.Y(), p2.Y()); }





  VEC2D Delta () const { return VEC2D (p2.X()-p1.X(), p2.Y()-p1.Y()); }

  double Length () const { return Delta().Length(); }

  double Length2 () const

  {
    return sqr (p1.X() - p2.X()) +

           sqr (p1.Y() - p2.Y());
  }



  friend POINT2D CrossPoint (const LINE2D & l1, const LINE2D & l2);

  friend int Parallel (const LINE2D & l1, const LINE2D & l2, double peps = EPSGEOM);

  friend int IsOnLine (const LINE2D & l, const POINT2D & p, double heps = EPSGEOM);

  friend int IsOnLongLine (const LINE2D & l, const POINT2D & p);

  friend int Hit (const LINE2D & l1, const LINE2D & l2, double heps = EPSGEOM);



  friend ostream & operator<<(ostream  & s, const LINE2D & l);

};



class PLINE2D

{

protected:

  POINT2D const * p1, *p2;



public:

  PLINE2D() { };

  PLINE2D(POINT2D const * ap1, POINT2D const * ap2)

  { p1 = ap1; p2 = ap2; }



  PLINE2D & operator= (const PLINE2D & l2)

  { p1 = l2.p1; p2 = l2.p2; return *this;}



  const POINT2D *& P1() { return p1; }

  const POINT2D *& P2() { return p2; }

  const POINT2D & P1() const { return *p1; }

  const POINT2D & P2() const { return *p2; }



  double XMax() const { return max (p1->X(), p2->X()); }

  double YMax() const { return max (p1->Y(), p2->Y()); }

  double XMin() const { return min (p1->X(), p2->X()); }

  double YMin() const { return min (p1->Y(), p2->Y()); }





  VEC2D Delta () const { return VEC2D (p2->X()-p1->X(), p2->Y()-p1->Y()); }

  double Length () const { return Delta().Length(); }

  double Length2 () const

  {
    return sqr (p1->X() - p2->X()) +

           sqr (p1->Y() - p2->Y());
  }



  friend POINT2D CrossPoint (const PLINE2D & l1, const PLINE2D & l2);

  friend int Parallel (const PLINE2D & l1, const PLINE2D & l2, double peps = EPSGEOM);

  friend int IsOnLine (const PLINE2D & l, const POINT2D & p, double heps = EPSGEOM);

  friend int IsOnLongLine (const PLINE2D & l, const POINT2D & p);

  friend int Hit (const PLINE2D & l1, const LINE2D & l2, double heps = EPSGEOM);



  friend ostream & operator<<(ostream  & s, const LINE2D & l);

};







class ILINE

{

  INDEX i[2];



public:

  ILINE() {};

  ILINE(INDEX i1, INDEX i2) { i[0] = i1; i[1] = i2; }

  ILINE(const ILINE & l) { i[0] = l.i[0]; i[1] = l.i[1]; }



  ILINE & operator= (const ILINE & l)

  { i[0] = l.i[0]; i[1] = l.i[1]; return *this; }



  const INDEX & I(int ai) const { return i[ai-1]; }

  const INDEX & X() const { return i[0]; }

  const INDEX & Y() const { return i[1]; }

  const INDEX & I1() const { return i[0]; }

  const INDEX & I2() const { return i[1]; }



  INDEX & I(int ai) { return i[ai-1]; }

  INDEX & X() { return i[0]; }

  INDEX & Y() { return i[1]; }

  INDEX & I1() { return i[0]; }

  INDEX & I2() { return i[1]; }

};









class TRIANGLE2D

{

private:

  POINT2D p1, p2, p3;



public:

  TRIANGLE2D() { };

  TRIANGLE2D (const POINT2D & ap1, const POINT2D & ap2,

              const POINT2D & ap3)

  { p1 = ap1; p2 = ap2; p3 = ap3;}



  TRIANGLE2D & operator= (const TRIANGLE2D & t2)

  { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }



  POINT2D & P1() { return p1; }

  POINT2D & P2() { return p2; }

  POINT2D & P3() { return p3; }

  const POINT2D & P1() const { return p1; }

  const POINT2D & P2() const { return p2; }

  const POINT2D & P3() const { return p3; }



  double XMax() const { return max (p1.X(), p2.X(), p3.X()); }

  double YMax() const { return max (p1.Y(), p2.Y(), p3.Y()); }

  double XMin() const { return min (p1.X(), p2.X(), p3.X()); }

  double YMin() const { return min (p1.Y(), p2.Y(), p3.Y()); }



  inline POINT2D Center () const

  { return POINT2D( (p1.X()+p2.X()+p3.X())/3, (p1.Y()+p2.Y()+p3.Y())/3); }



  int Regular() const;

  int CW () const;

  int CCW () const;



  int IsOn (const POINT2D & p) const;

  int IsIn (const POINT2D & p) const;

  friend ostream & operator<<(ostream  & s, const TRIANGLE2D & t);

};





class PTRIANGLE2D

{

private:

  POINT2D const *p1, *p2, *p3;



public:

  PTRIANGLE2D() { };

  PTRIANGLE2D (const POINT2D * ap1, const POINT2D * ap2,

               const POINT2D * ap3)

  { p1 = ap1; p2 = ap2; p3 = ap3;}



  PTRIANGLE2D & operator= (const PTRIANGLE2D & t2)

  { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }



  const POINT2D *& P1() { return p1; }

  const POINT2D *& P2() { return p2; }

  const POINT2D *& P3() { return p3; }

  const POINT2D * P1() const { return p1; }

  const POINT2D * P2() const { return p2; }

  const POINT2D * P3() const { return p3; }



  double XMax() const { return max (p1->X(), p2->X(), p3->X()); }

  double YMax() const { return max (p1->Y(), p2->Y(), p3->Y()); }

  double XMin() const { return min (p1->X(), p2->X(), p3->X()); }

  double YMin() const { return min (p1->Y(), p2->Y(), p3->Y()); }



  POINT2D Center () const

  { return POINT2D( (p1->X()+p2->X()+p3->X())/3, (p1->Y()+p2->Y()+p3->Y())/3); }





  int Regular() const;

  int CW () const;

  int CCW () const;



  int IsOn (const POINT2D & p) const;

  int IsIn (const POINT2D & p) const;

  friend ostream & operator<<(ostream & s, const PTRIANGLE2D & t);

};







/*

   class POLYGON2D

   {

   protected:

   PLIST<POINT2D> points;



   public:

   POLYGON2D () : points() { };



   void Insert (POINT2D * p) { points.Insert(p); }



   double Area () const { return fabs (HArea()); }

   char CW () const { return HArea() > 0; }

   char CCW () const { return HArea() < 0; }



   char IsOn (const POINT2D & p) const;

   char IsIn (const POINT2D & p) const;

   char Convex () const;

   char IsStarPoint (const POINT2D & p) const;

   POINT2D Center() const;

   POINT2D EqualAreaPoint () const;

   private:

   double HArea () const;

   };

 */





extern double Fastatan2 (double x, double y);











// Inline - Functions:



inline VEC2D & VEC2D :: operator+= (const VEC2D & v2)

{

  vx += v2.vx;

  vy += v2.vy;

  return *this;

}



inline VEC2D & VEC2D :: operator-= (const VEC2D & v2)

{

  vx -= v2.vx;

  vy -= v2.vy;

  return *this;

}



inline VEC2D & VEC2D :: operator*= (double s)

{

  vx *= s;

  vy *= s;

  return *this;

}



inline VEC2D & VEC2D :: operator/= (double s)

{

  if (s != 0)

  {

    vx /= s;

    vy /= s;

  }

  else

  {

    MyError ("VEC3D::operator /=: Divisioin by zero");

  }

  return *this;

}





inline VEC2D operator- (const POINT2D & p1, const POINT2D & p2)

{

  return VEC2D (p1.X() - p2.X(), p1.Y() - p2.Y());

}



inline POINT2D operator- (const POINT2D & p1, const VEC2D & v)

{

  return POINT2D (p1.X() - v.X(), p1.Y() - v.Y());

}



inline POINT2D operator+ (const POINT2D & p1, const VEC2D & v)

{

  return POINT2D (p1.X() + v.X(), p1.Y() + v.Y());

}



inline POINT2D Center (const POINT2D & p1, const POINT2D & p2)

{

  return POINT2D ((p1.X() + p2.X()) / 2, (p1.Y() + p2.Y()) / 2);

}



inline VEC2D operator- (const VEC2D & v1, const VEC2D & v2)

{

  return VEC2D (v1.X() - v2.X(), v1.Y() - v2.Y());

}



inline VEC2D operator+ (const VEC2D & v1, const VEC2D & v2)

{

  return VEC2D (v1.X() + v2.X(), v1.Y() + v2.Y());

}



inline VEC2D operator* (double scal, const VEC2D & v)

{

  return VEC2D (scal * v.X(), scal * v.Y());

}



inline void PpSmV (const POINT2D & p1, double s,

                   const VEC2D & v, POINT2D & p2)

{

  p2.X() = p1.X() + s * v.X();

  p2.Y() = p1.Y() + s * v.Y();

}



inline void PmP (const POINT2D & p1, const POINT2D & p2, VEC2D & v)

{

  v.X() = p1.X() - p2.X();

  v.Y() = p1.Y() - p2.Y();

}











inline int TRIANGLE2D :: Regular() const

{

  return fabs(Cross ( p2 - p1, p3 - p2)) > EPSGEOM;

}



inline int TRIANGLE2D :: CW () const

{

  return Cross ( p2 - p1, p3 - p2) < 0;

}



inline int TRIANGLE2D :: CCW () const

{

  return Cross ( p2 - p1, p3 - p2) > 0;

}







inline int PTRIANGLE2D :: Regular() const

{

  return fabs(Cross ( *p2 - *p1, *p3 - *p2)) > EPSGEOM;

}



inline int PTRIANGLE2D :: CW () const

{

  return Cross ( *p2 - *p1, *p3 - *p2) < 0;

}



inline int PTRIANGLE2D :: CCW () const

{

  return Cross ( *p2 - *p1, *p3 - *p2) > 0;

}
