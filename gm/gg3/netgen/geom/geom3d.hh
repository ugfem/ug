// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
extern void MyError (char * ch);
class ostream;

class Point3d;
class Vec3d;

class Point3d
{
protected:
  double px, py, pz;

public:
  Point3d() { px = py = pz = 0; }
  Point3d(double ax, double ay, double az) { px = ax; py = ay; pz = az; }
  Point3d(const Point3d & p2) { px = p2.px; py = p2.py; pz = p2.pz; }

  Point3d & operator= (const Point3d & p2)
  { px = p2.px; py = p2.py; pz = p2.pz; return *this; }

  double & X() { return px; }
  double & Y() { return py; }
  double & Z() { return pz; }
  double X() const { return px; }
  double Y() const { return py; }
  double Z() const { return pz; }

  friend inline Vec3d operator- (const Point3d & p1, const Point3d & p2);
  friend inline Point3d operator- (const Point3d & p1, const Vec3d & v);
  friend inline Point3d operator+ (const Point3d & p1, const Vec3d & v);
  friend double Dist (const Point3d & p1, const Point3d & p2)
  {
    return sqrt ( sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) +
                  sqr (p1.Z()-p2.Z()) );
  }
  friend double Dist2 (const Point3d & p1, const Point3d & p2)
  {
    return sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) +
           sqr (p1.Z()-p2.Z()) ;
  }
  friend inline Point3d Center (const Point3d & p1, const Point3d & p2);
  friend inline Point3d Center (const Point3d & p1, const Point3d & p2, const Point3d & p3);
  friend ostream & operator<<(ostream  & s, const Point3d & p);
};


class Vec3d
{
protected:
  double vx, vy, vz;

public:
  Vec3d() { vx = vy = vz = 0; }
  Vec3d(double ax, double ay, double az) { vx = ax; vy = ay; vz = az; }
  Vec3d(Vec3d & v2) { vx = v2.vx; vy = v2.vy; vz = v2.vz; }

  Vec3d & operator= (const Vec3d & v2)
  { vx = v2.vx; vy = v2.vy; vz = v2.vz; return *this; }

  double & X() { return vx; }
  double & Y() { return vy; }
  double & Z() { return vz; }
  double X() const { return vx; }
  double Y() const { return vy; }
  double Z() const { return vz; }

  double Length() const { return sqrt (vx * vx + vy * vy + vz * vz); }
  double Length2() const { return vx * vx + vy * vy + vz * vz; }

  Vec3d & operator+= (const Vec3d & v2);
  Vec3d & operator-= (const Vec3d & v2);
  Vec3d & operator*= (double s);
  Vec3d & operator/= (double s);

  friend inline Vec3d operator- (const Point3d & p1, const Point3d & p2);
  friend inline Point3d operator- (const Point3d & p1, const Vec3d & v);
  friend inline Point3d operator+ (const Point3d & p1, const Vec3d & v);
  friend inline Vec3d operator- (const Vec3d & p1, const Vec3d & v);
  friend inline Vec3d operator+ (const Vec3d & p1, const Vec3d & v);
  friend inline Vec3d operator* (double scal, const Vec3d & v);

  friend inline double operator* (const Vec3d & v1, const Vec3d & v2);
  friend inline Vec3d Cross (const Vec3d & v1, const Vec3d & v2);

  void GetNormal (Vec3d & n);
  friend double Angle (const Vec3d & v);
  friend double FastAngle (const Vec3d & v);
  friend double Angle (const Vec3d & v1, const Vec3d & v2);
  friend double FastAngle (const Vec3d & v1, const Vec3d & v2);

  friend ostream & operator<<(ostream  & s, const Vec3d & v);
};




inline Point3d Center (const Point3d & p1, const Point3d & p2)
{
  return Point3d (0.5 * (p1.X() + p2.X()),
                  0.5 * (p1.Y() + p2.Y()),
                  0.5 * (p1.Z() + p2.Z()));
}

inline Point3d Center (const Point3d & p1, const Point3d & p2,
                       const Point3d & p3)
{
  return Point3d (1.0/3.0 * (p1.X() + p2.X() + p3.X()),
                  1.0/3.0 * (p1.Y() + p2.Y() + p3.Y()),
                  1.0/3.0 * (p1.Z() + p2.Z() + p3.Z()));
}


inline Vec3d & Vec3d :: operator+= (const Vec3d & v2)
{
  vx += v2.vx;
  vy += v2.vy;
  vz += v2.vz;
  return *this;
}

inline Vec3d & Vec3d :: operator-= (const Vec3d & v2)
{
  vx -= v2.vx;
  vy -= v2.vy;
  vz -= v2.vz;
  return *this;
}

inline Vec3d & Vec3d :: operator*= (double s)
{
  vx *= s;
  vy *= s;
  vz *= s;
  return *this;
}

inline Vec3d & Vec3d :: operator/= (double s)
{
  if (s != 0)
  {
    vx /= s;
    vy /= s;
    vz /= s;
  }
  else
  {
    MyError ("Vec3d::operator /=: Divisioin by zero");
  }
  return *this;
}





inline Vec3d operator- (const Point3d & p1, const Point3d & p2)
{
  return Vec3d (p1.X() - p2.X(), p1.Y() - p2.Y(),p1.Z() - p2.Z());
}

inline Point3d operator- (const Point3d & p1, const Vec3d & v)
{
  return Point3d (p1.X() - v.X(), p1.Y() - v.Y(),p1.Z() - v.Z());
}

inline Point3d operator+ (const Point3d & p1, const Vec3d & v)
{
  return Point3d (p1.X() + v.X(), p1.Y() + v.Y(),p1.Z() + v.Z());
}

inline Vec3d operator- (const Vec3d & v1, const Vec3d & v2)
{
  return Vec3d (v1.X() - v2.X(), v1.Y() - v2.Y(),v1.Z() - v2.Z());
}

inline Vec3d operator+ (const Vec3d & v1, const Vec3d & v2)
{
  return Vec3d (v1.X() + v2.X(), v1.Y() + v2.Y(),v1.Z() + v2.Z());
}

inline Vec3d operator* (double scal, const Vec3d & v)
{
  return Vec3d (scal * v.X(), scal * v.Y(), scal * v.Z());
}


inline double operator* (const Vec3d & v1, const Vec3d & v2)
{
  return v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z();
}


inline Vec3d Cross (const Vec3d & v1, const Vec3d & v2)
{
  return Vec3d
           ( v1.Y() * v2.Z() - v1.Z() * v2.Y(),
           v1.Z() * v2.X() - v1.X() * v2.Z(),
           v1.X() * v2.Y() - v1.Y() * v2.X());
}


class box3d
{
  double minx, maxx, miny, maxy, minz, maxz;
  double diam, inner;
  Point3d c;


public:

  box3d () { };

  box3d ( double aminx, double amaxx,
          double aminy, double amaxy,
          double aminz, double amaxz );


  const Point3d & Center () const { return c; }
  double Diam () const { return diam; }
  double Inner () const { return inner; }
  double MinX () const { return minx; }
  double MaxX () const { return maxx; }
  double MinY () const { return miny; }
  double MaxY () const { return maxy; }
  double MinZ () const { return minz; }
  double MaxZ () const { return maxz; }

  void CalcSubBox (int i, box3d & sbox) const;
  void GetPointNr (int i, Point3d & point) const;

private:
  void CalcDiamCenter ();
};
