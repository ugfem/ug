// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* Geometric Algorithms */

#define EPSGEOM 1E-5

extern void MyError (char * ch);
class ostream;

class Point2d;
class Vec2d;
class LINE2D;



class Point2d
{
  friend Vec2d;

protected:
  double px, py;

public:
  Point2d() { /* px = py = 0; */ }
  Point2d(double ax, double ay) { px = ax; py = ay; }
  Point2d(const Point2d & p2) { px = p2.px; py = p2.py; }

  Point2d & operator= (const Point2d & p2)
  { px = p2.px; py = p2.py; return *this; }

  double & X() { return px; }
  double & Y() { return py; }
  double X() const { return px; }
  double Y() const { return py; }

  friend inline Vec2d operator- (const Point2d & p1, const Point2d & p2);
  friend inline Point2d operator- (const Point2d & p1, const Vec2d & v);
  friend inline Point2d operator+ (const Point2d & p1, const Vec2d & v);

  friend inline Point2d Center (const Point2d & p1, const Point2d & p2);

  friend double Dist (const Point2d & p1, const Point2d & p2)
  { return sqrt ( sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ); }
  friend double Dist2 (const Point2d & p1, const Point2d & p2)
  { return sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ; }

  friend int CW (const Point2d & p1, const Point2d & p2, const Point2d & p3);
  friend int CCW (const Point2d & p1, const Point2d & p2, const Point2d & p3);

  friend inline void PpSmV (const Point2d & p1, double s, const Vec2d & v, Point2d & p2);
  friend inline void PmP (const Point2d & p1, const Point2d & p2, Vec2d & v);

  friend ostream & operator<<(ostream  & s, const Point2d & p);
};



class Vec2d
{
protected:
  double vx, vy;

public:
  Vec2d() { /* vx = vy = 0; */ }
  Vec2d(double ax, double ay)
  { vx = ax; vy = ay; }
  Vec2d(const Vec2d & v2)
  { vx = v2.vx; vy = v2.vy; }
  Vec2d(const Point2d & p1, const Point2d & p2)
  { vx = p2.px - p1.px; vy = p2.py - p1.py; }

  Vec2d & operator= (const Vec2d & p2)
  { vx = p2.vx; vy = p2.vy; return *this; }

  double & X() { return vx; }
  double & Y() { return vy; }
  double X() const { return vx; }
  double Y() const { return vy; }

  double Length() const { return sqrt (vx * vx + vy * vy); }
  double Length2() const { return vx * vx + vy * vy; }

  inline Vec2d & operator+= (const Vec2d & v2);
  inline Vec2d & operator-= (const Vec2d & v2);
  inline Vec2d & operator*= (double s);
  inline Vec2d & operator/= (double s);

  friend inline Vec2d operator- (const Point2d & p1, const Point2d & p2);
  friend inline Point2d operator- (const Point2d & p1, const Vec2d & v);
  friend inline Point2d operator+ (const Point2d & p1, const Vec2d & v);
  friend inline Vec2d operator- (const Vec2d & p1, const Vec2d & v);
  friend inline Vec2d operator+ (const Vec2d & p1, const Vec2d & v);
  friend inline Vec2d operator* (double scal, const Vec2d & v);

  friend double operator* (const Vec2d & v1, const Vec2d & v2)
  { return v1.X() * v2.X() + v1.Y() * v2.Y(); }


  friend double Cross (const Vec2d & v1, const Vec2d & v2)
  {
    return double(v1.X()) * double(v2.Y()) -
           double(v1.Y()) * double(v2.X());
  }

  friend inline void PpSmV (const Point2d & p1, double s, const Vec2d & v, Point2d & p2);
  friend inline void PmP (const Point2d & p1, const Point2d & p2, Vec2d & v);


  friend double Angle (const Vec2d & v);
  friend double FastAngle (const Vec2d & v);
  friend double Angle (const Vec2d & v1, const Vec2d & v2);
  friend double FastAngle (const Vec2d & v1, const Vec2d & v2);

  friend ostream & operator<<(ostream  & s, const Vec2d & v);
};



class LINE2D
{
protected:
  Point2d p1, p2;

public:
  LINE2D() : p1(), p2() { };
  LINE2D(const Point2d & ap1, const Point2d & ap2)
  { p1 = ap1; p2 = ap2; }

  LINE2D & operator= (const LINE2D & l2)
  { p1 = l2.p1; p2 = l2.p2; return *this;}

  Point2d & P1() { return p1; }
  Point2d & P2() { return p2; }
  const Point2d & P1() const { return p1; }
  const Point2d & P2() const { return p2; }

  double XMax() const { return max (p1.X(), p2.X()); }
  double YMax() const { return max (p1.Y(), p2.Y()); }
  double XMin() const { return min (p1.X(), p2.X()); }
  double YMin() const { return min (p1.Y(), p2.Y()); }


  Vec2d Delta () const { return Vec2d (p2.X()-p1.X(), p2.Y()-p1.Y()); }
  double Length () const { return Delta().Length(); }
  double Length2 () const
  {
    return sqr (p1.X() - p2.X()) +
           sqr (p1.Y() - p2.Y());
  }

  friend Point2d CrossPoint (const LINE2D & l1, const LINE2D & l2);
  friend int Parallel (const LINE2D & l1, const LINE2D & l2, double peps = EPSGEOM);
  friend int IsOnLine (const LINE2D & l, const Point2d & p, double heps = EPSGEOM);
  friend int IsOnLongLine (const LINE2D & l, const Point2d & p);
  friend int Hit (const LINE2D & l1, const LINE2D & l2, double heps = EPSGEOM);

  friend ostream & operator<<(ostream  & s, const LINE2D & l);
};

class PLINE2D
{
protected:
  Point2d const * p1, *p2;

public:
  PLINE2D() { };
  PLINE2D(Point2d const * ap1, Point2d const * ap2)
  { p1 = ap1; p2 = ap2; }

  PLINE2D & operator= (const PLINE2D & l2)
  { p1 = l2.p1; p2 = l2.p2; return *this;}

  const Point2d *& P1() { return p1; }
  const Point2d *& P2() { return p2; }
  const Point2d & P1() const { return *p1; }
  const Point2d & P2() const { return *p2; }

  double XMax() const { return max (p1->X(), p2->X()); }
  double YMax() const { return max (p1->Y(), p2->Y()); }
  double XMin() const { return min (p1->X(), p2->X()); }
  double YMin() const { return min (p1->Y(), p2->Y()); }


  Vec2d Delta () const { return Vec2d (p2->X()-p1->X(), p2->Y()-p1->Y()); }
  double Length () const { return Delta().Length(); }
  double Length2 () const
  {
    return sqr (p1->X() - p2->X()) +
           sqr (p1->Y() - p2->Y());
  }

  friend Point2d CrossPoint (const PLINE2D & l1, const PLINE2D & l2);
  friend int Parallel (const PLINE2D & l1, const PLINE2D & l2, double peps = EPSGEOM);
  friend int IsOnLine (const PLINE2D & l, const Point2d & p, double heps = EPSGEOM);
  friend int IsOnLongLine (const PLINE2D & l, const Point2d & p);
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
  Point2d p1, p2, p3;

public:
  TRIANGLE2D() { };
  TRIANGLE2D (const Point2d & ap1, const Point2d & ap2,
              const Point2d & ap3)
  { p1 = ap1; p2 = ap2; p3 = ap3;}

  TRIANGLE2D & operator= (const TRIANGLE2D & t2)
  { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }

  Point2d & P1() { return p1; }
  Point2d & P2() { return p2; }
  Point2d & P3() { return p3; }
  const Point2d & P1() const { return p1; }
  const Point2d & P2() const { return p2; }
  const Point2d & P3() const { return p3; }

  double XMax() const { return max (p1.X(), p2.X(), p3.X()); }
  double YMax() const { return max (p1.Y(), p2.Y(), p3.Y()); }
  double XMin() const { return min (p1.X(), p2.X(), p3.X()); }
  double YMin() const { return min (p1.Y(), p2.Y(), p3.Y()); }

  inline Point2d Center () const
  { return Point2d( (p1.X()+p2.X()+p3.X())/3, (p1.Y()+p2.Y()+p3.Y())/3); }

  int Regular() const;
  int CW () const;
  int CCW () const;

  int IsOn (const Point2d & p) const;
  int IsIn (const Point2d & p) const;
  friend ostream & operator<<(ostream  & s, const TRIANGLE2D & t);
};


class PTRIANGLE2D
{
private:
  Point2d const *p1, *p2, *p3;

public:
  PTRIANGLE2D() { };
  PTRIANGLE2D (const Point2d * ap1, const Point2d * ap2,
               const Point2d * ap3)
  { p1 = ap1; p2 = ap2; p3 = ap3;}

  PTRIANGLE2D & operator= (const PTRIANGLE2D & t2)
  { p1 = t2.p1; p2 = t2.p2; p3 = t2.p3; return *this; }

  const Point2d *& P1() { return p1; }
  const Point2d *& P2() { return p2; }
  const Point2d *& P3() { return p3; }
  const Point2d * P1() const { return p1; }
  const Point2d * P2() const { return p2; }
  const Point2d * P3() const { return p3; }

  double XMax() const { return max (p1->X(), p2->X(), p3->X()); }
  double YMax() const { return max (p1->Y(), p2->Y(), p3->Y()); }
  double XMin() const { return min (p1->X(), p2->X(), p3->X()); }
  double YMin() const { return min (p1->Y(), p2->Y(), p3->Y()); }

  Point2d Center () const
  { return Point2d( (p1->X()+p2->X()+p3->X())/3, (p1->Y()+p2->Y()+p3->Y())/3); }


  int Regular() const;
  int CW () const;
  int CCW () const;

  int IsOn (const Point2d & p) const;
  int IsIn (const Point2d & p) const;
  friend ostream & operator<<(ostream & s, const PTRIANGLE2D & t);
};



/*
   class POLYGON2D
   {
   protected:
   PLIST<Point2d> points;

   public:
   POLYGON2D () : points() { };

   void Insert (Point2d * p) { points.Insert(p); }

   double Area () const { return fabs (HArea()); }
   char CW () const { return HArea() > 0; }
   char CCW () const { return HArea() < 0; }

   char IsOn (const Point2d & p) const;
   char IsIn (const Point2d & p) const;
   char Convex () const;
   char IsStarPoint (const Point2d & p) const;
   Point2d Center() const;
   Point2d EqualAreaPoint () const;
   private:
   double HArea () const;
   };
 */


extern double Fastatan2 (double x, double y);





// Inline - Functions:

inline Vec2d & Vec2d :: operator+= (const Vec2d & v2)
{
  vx += v2.vx;
  vy += v2.vy;
  return *this;
}

inline Vec2d & Vec2d :: operator-= (const Vec2d & v2)
{
  vx -= v2.vx;
  vy -= v2.vy;
  return *this;
}

inline Vec2d & Vec2d :: operator*= (double s)
{
  vx *= s;
  vy *= s;
  return *this;
}

inline Vec2d & Vec2d :: operator/= (double s)
{
  if (s != 0)
  {
    vx /= s;
    vy /= s;
  }
  else
  {
    MyError ("Vec3d::operator /=: Divisioin by zero");
  }
  return *this;
}


inline Vec2d operator- (const Point2d & p1, const Point2d & p2)
{
  return Vec2d (p1.X() - p2.X(), p1.Y() - p2.Y());
}

inline Point2d operator- (const Point2d & p1, const Vec2d & v)
{
  return Point2d (p1.X() - v.X(), p1.Y() - v.Y());
}

inline Point2d operator+ (const Point2d & p1, const Vec2d & v)
{
  return Point2d (p1.X() + v.X(), p1.Y() + v.Y());
}

inline Point2d Center (const Point2d & p1, const Point2d & p2)
{
  return Point2d ((p1.X() + p2.X()) / 2, (p1.Y() + p2.Y()) / 2);
}

inline Vec2d operator- (const Vec2d & v1, const Vec2d & v2)
{
  return Vec2d (v1.X() - v2.X(), v1.Y() - v2.Y());
}

inline Vec2d operator+ (const Vec2d & v1, const Vec2d & v2)
{
  return Vec2d (v1.X() + v2.X(), v1.Y() + v2.Y());
}

inline Vec2d operator* (double scal, const Vec2d & v)
{
  return Vec2d (scal * v.X(), scal * v.Y());
}

inline void PpSmV (const Point2d & p1, double s,
                   const Vec2d & v, Point2d & p2)
{
  p2.X() = p1.X() + s * v.X();
  p2.Y() = p1.Y() + s * v.Y();
}

inline void PmP (const Point2d & p1, const Point2d & p2, Vec2d & v)
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
