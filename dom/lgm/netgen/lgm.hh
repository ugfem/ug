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

class InputElement : public Element
{
  int neighbour[4];
public:
  int & Neighbour (int i) { return neighbour[i-1]; }
};

class surfacemeshing
{
  ADFRONT2 * adfront;
  ARRAY<netrule*> rules;
  ARRAY<int> ruleused;
  double cxx, cyy, czz, cxy, cxz, cyz, cx, cy, cz, c1;
  Vec3d ex, ey, ez;
  Point3d globp1;

public:
  surfacemeshing (char * rulefilename);
  virtual ~surfacemeshing ();

  void LoadRules (char * filename);
  void Mesh (double gh);

  void ImproveMesh (ARRAY<Point3d> & points, const ARRAY<Element> & elements,
                    int improveedges, int numboundarypoints, double h, int steps, int err2);

  void AddPoint (const Point3d & p, INDEX globind);
  void AddBoundaryElement (INDEX i1, INDEX i2, int surfind);
  virtual void TestPoint (const Point3d & /* p */,int flag) { };

  virtual int Write_Surface_Grid ();
  virtual int SavePoint (const Point3d & p);
  virtual void SaveElement (const Element & elem);
  virtual void Write2Shell (int n);

  //  friend int StartNetgen (double h, int smooth, int display);

protected:
  virtual void StartMesh ();
  virtual void EndMesh ();
  virtual int DefineTransformation (INDEX surfind, Point3d & p1, Point3d & p2);
  virtual void TransformToPlain (INDEX ind, const Point3d & locpoint,
                                 Point2d & plainpoint, double h);
  virtual void TransformFromPlain (INDEX surfind, Point2d & plainpoint,
                                   Point3d & locpoint, double h);
public:
  virtual void ProjectPoint (INDEX surfind, Point3d & p) const;
  virtual void ProjectPoint2 (INDEX surfind, INDEX surfind2, Point3d & p) const;
  virtual void GetNormalVector(INDEX surfind, const Point3d & p, Vec3d & n) const;

  virtual double CalcLocalH (const Point3d & p, int surfind, double gh) const;
};


void LoadGeo (char * filename, ARRAY<geompoint3d> & geompoints,ARRAY<splinesegment3d*> & splines, double & elto0);
int GetEdgeId(const Point3d & ep1, Point3d & ep2);
int Calc_Coord_Vectors(const Point3d p1,const Point3d p2,const int mi,Vec3d & nx,Vec3d & ny,Vec3d & nz);
void Project_Point2Surface(Point3d &inpoint, Point3d &outpoint);
int GenerateTriangle(ARRAY<Point2d> & lpoints, ARRAY<ILINE> & llines,ARRAY<Element> & elements, ARRAY<INDEX> & dellines,double h);
void BFGS (Vector & x, double (*f)(const Vector & x, Vector & g));
int Project_Point2Surface_2(Point3d &inpoint, Point3d &outpoint, Vec3d n);
int GetEdgeId(const Point3d & ep1, Point3d & ep2);
void Smooth_SurfaceMesh (ARRAY<Point3d> & points, const ARRAY<Element> & elements, int numboundarypoints, int steps);
Vec3d NormalVector(InputElement e);
double Project2Plane(Point3d & p,Vec3d & np,Point3d & p0,Vec3d & n0,Vec3d & n1,Vec3d & n2);
int Calc_Vectors(const Point3d p0,const Point3d p1,const Point3d p2,Vec3d & n0,Vec3d & n1,Vec3d & n2);
double Calc_Angle(InputElement e1, InputElement e2);
double Calc_Local_Coordinates(const Point3d p0,const Point3d p1,const Point3d p2,const Point3d p,double & lam0,double & lam1,double & lam2);
