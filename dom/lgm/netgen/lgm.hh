// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
class geompoint3d
{
public:
  Point3d p;
  double refatpoint;
};

class splinesegment3d
{
public:
  virtual void Draw (const ROT3D & r) const;
  virtual double Length () const;
  virtual void GetPoint (double t, Point3d & p) const = 0;
  void Partition (double h, double elto0,
                  ARRAY<Point3d> & points, ARRAY<int> & lp1, ARRAY<int> & lp2) const;
  /*  void Partition (int intervalls,
        ARRAY<Point3d> & points, ARRAY<int> & lp1, ARRAY<int> & lp2) const;
     virtual geompoint3d * StartPI () const = 0;
     virtual geompoint3d * EndPI () const = 0;*/
};

class linesegment3d : public splinesegment3d
{
  geompoint3d *p1, *p2;
public:
  linesegment3d (geompoint3d * ap1, geompoint3d * ap2);
  virtual void Draw (const ROT3D & r) const;
  virtual double Length () const;
  virtual void GetPoint (double t, Point3d & p) const;
  /*  virtual geompoint3d * StartPI () const { return p1; };
     virtual geompoint3d * EndPI () const { return p2; }*/
};

/*class INPUTElement : public Element
   {
   int neighbour[4];
   public:
   int & Neighbour (int i) { return neighbour[i-1]; }

   }*/

extern void LoadGeo (char * filename, ARRAY<geompoint3d> & geompoints,
                     ARRAY<splinesegment3d*> & splines, double & elto0);
