// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
extern void MyError (char * ch);

class ostream;



class POINT3D;

class VEC3D;



class POINT3D

{

protected:

  double px, py, pz;



public:

  POINT3D() { px = py = pz = 0; }

  POINT3D(double ax, double ay, double az) { px = ax; py = ay; pz = az; }

  POINT3D(const POINT3D & p2) { px = p2.px; py = p2.py; pz = p2.pz; }



  POINT3D & operator= (const POINT3D & p2)

  { px = p2.px; py = p2.py; pz = p2.pz; return *this; }



  double & X() { return px; }

  double & Y() { return py; }

  double & Z() { return pz; }

  double X() const { return px; }

  double Y() const { return py; }

  double Z() const { return pz; }



  friend inline VEC3D operator- (const POINT3D & p1, const POINT3D & p2);

  friend inline POINT3D operator- (const POINT3D & p1, const VEC3D & v);

  friend inline POINT3D operator+ (const POINT3D & p1, const VEC3D & v);

  friend double Dist (const POINT3D & p1, const POINT3D & p2)

  {
    return sqrt ( sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) +

                  sqr (p1.Z()-p2.Z()) );
  }

  friend double Dist2 (const POINT3D & p1, const POINT3D & p2)

  {
    return sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) +

           sqr (p1.Z()-p2.Z()) ;
  }

  friend inline POINT3D Center (const POINT3D & p1, const POINT3D & p2);

  friend inline POINT3D Center (const POINT3D & p1, const POINT3D & p2, const POINT3D & p3);

  friend ostream & operator<<(ostream  & s, const POINT3D & p);

};





class VEC3D

{

protected:

  double vx, vy, vz;



public:

  VEC3D() { vx = vy = vz = 0; }

  VEC3D(double ax, double ay, double az) { vx = ax; vy = ay; vz = az; }

  VEC3D(VEC3D & v2) { vx = v2.vx; vy = v2.vy; vz = v2.vz; }



  VEC3D & operator= (const VEC3D & v2)

  { vx = v2.vx; vy = v2.vy; vz = v2.vz; return *this; }



  double & X() { return vx; }

  double & Y() { return vy; }

  double & Z() { return vz; }

  double X() const { return vx; }

  double Y() const { return vy; }

  double Z() const { return vz; }



  double Length() const { return sqrt (vx * vx + vy * vy + vz * vz); }

  double Length2() const { return vx * vx + vy * vy + vz * vz; }



  VEC3D & operator+= (const VEC3D & v2);

  VEC3D & operator-= (const VEC3D & v2);

  VEC3D & operator*= (double s);

  VEC3D & operator/= (double s);



  friend inline VEC3D operator- (const POINT3D & p1, const POINT3D & p2);

  friend inline POINT3D operator- (const POINT3D & p1, const VEC3D & v);

  friend inline POINT3D operator+ (const POINT3D & p1, const VEC3D & v);

  friend inline VEC3D operator- (const VEC3D & p1, const VEC3D & v);

  friend inline VEC3D operator+ (const VEC3D & p1, const VEC3D & v);

  friend inline VEC3D operator* (double scal, const VEC3D & v);



  friend inline double operator* (const VEC3D & v1, const VEC3D & v2);

  friend inline VEC3D Cross (const VEC3D & v1, const VEC3D & v2);



  friend double Angle (const VEC3D & v);

  friend double FastAngle (const VEC3D & v);

  friend double Angle (const VEC3D & v1, const VEC3D & v2);

  friend double FastAngle (const VEC3D & v1, const VEC3D & v2);



  friend ostream & operator<<(ostream  & s, const VEC3D & v);

};









inline POINT3D Center (const POINT3D & p1, const POINT3D & p2)

{

  return POINT3D (0.5 * (p1.X() + p2.X()),

                  0.5 * (p1.Y() + p2.Y()),

                  0.5 * (p1.Z() + p2.Z()));

}



inline POINT3D Center (const POINT3D & p1, const POINT3D & p2,

                       const POINT3D & p3)

{

  return POINT3D (1.0/3.0 * (p1.X() + p2.X() + p3.X()),

                  1.0/3.0 * (p1.Y() + p2.Y() + p3.Y()),

                  1.0/3.0 * (p1.Z() + p2.Z() + p3.Z()));

}





inline VEC3D & VEC3D :: operator+= (const VEC3D & v2)

{

  vx += v2.vx;

  vy += v2.vy;

  vz += v2.vz;

  return *this;

}



inline VEC3D & VEC3D :: operator-= (const VEC3D & v2)

{

  vx -= v2.vx;

  vy -= v2.vy;

  vz -= v2.vz;

  return *this;

}



inline VEC3D & VEC3D :: operator*= (double s)

{

  vx *= s;

  vy *= s;

  vz *= s;

  return *this;

}



inline VEC3D & VEC3D :: operator/= (double s)

{

  if (s != 0)

  {

    vx /= s;

    vy /= s;

    vz /= s;

  }

  else

  {

    MyError ("VEC3D::operator /=: Divisioin by zero");

  }

  return *this;

}











inline VEC3D operator- (const POINT3D & p1, const POINT3D & p2)

{

  return VEC3D (p1.X() - p2.X(), p1.Y() - p2.Y(),p1.Z() - p2.Z());

}



inline POINT3D operator- (const POINT3D & p1, const VEC3D & v)

{

  return POINT3D (p1.X() - v.X(), p1.Y() - v.Y(),p1.Z() - v.Z());

}



inline POINT3D operator+ (const POINT3D & p1, const VEC3D & v)

{

  return POINT3D (p1.X() + v.X(), p1.Y() + v.Y(),p1.Z() + v.Z());

}



inline VEC3D operator- (const VEC3D & v1, const VEC3D & v2)

{

  return VEC3D (v1.X() - v2.X(), v1.Y() - v2.Y(),v1.Z() - v2.Z());

}



inline VEC3D operator+ (const VEC3D & v1, const VEC3D & v2)

{

  return VEC3D (v1.X() + v2.X(), v1.Y() + v2.Y(),v1.Z() + v2.Z());

}



inline VEC3D operator* (double scal, const VEC3D & v)

{

  return VEC3D (scal * v.X(), scal * v.Y(), scal * v.Z());

}





inline double operator* (const VEC3D & v1, const VEC3D & v2)

{

  return v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z();

}





inline VEC3D Cross (const VEC3D & v1, const VEC3D & v2)

{

  return VEC3D

           ( v1.Y() * v2.Z() - v1.Z() * v2.Y(),

           v1.Z() * v2.X() - v1.X() * v2.Z(),

           v1.X() * v2.Y() - v1.Y() * v2.X());

}





class box3d

{

  double minx, maxx, miny, maxy, minz, maxz;

  double diam, inner;

  POINT3D c;





public:



  box3d () { };



  box3d ( double aminx, double amaxx,

          double aminy, double amaxy,

          double aminz, double amaxz );





  const POINT3D & Center () const { return c; }

  double Diam () const { return diam; }

  double Inner () const { return inner; }

  double MinX () const { return minx; }

  double MaxX () const { return maxx; }

  double MinY () const { return miny; }

  double MaxY () const { return maxy; }

  double MinZ () const { return minz; }

  double MaxZ () const { return maxz; }



  void CalcSubBox (int i, box3d & sbox) const;

  void GetPointNr (int i, POINT3D & point) const;



private:

  void CalcDiamCenter ();

};
