// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_MESHING2
#define FILE_MESHING2

#include <meshing/meshtype.hh>
#include <meshing/ruler2.hh>
#include <meshing/adfront2.hh>

class Meshing2
{
  ADFRONT2 * adfront;
  ARRAY<netrule*> rules;
  ARRAY<int> ruleused;

public:
  Meshing2 (char * rulefilename);
  virtual ~Meshing2 ();

  void LoadRules (char * filename);
  void Mesh (double gh);

  void ImproveMesh (ARRAY<Point3d> & points, const ARRAY<Element> & elements,
                    int improveedges, int numboundarypoints, double h, int steps, int err2);

  void AddPoint (const Point3d & p, INDEX globind);
  void AddBoundaryElement (INDEX i1, INDEX i2, int surfind);

  virtual int SavePoint (const Point3d & /* p */) { return 0; }
  virtual void SaveElement (const Element & /* elem */) { };
  virtual void TestPoint (const Point3d & /* p */,int flag) { };

  virtual void DrawTransformation (const ROT3D & /* rot */) const { };

protected:
  virtual void StartMesh ();
  virtual void EndMesh ();
  virtual double CalcLocalH (const Point3d & p, int surfind, double gh) const;

  virtual void DefineTransformation (INDEX surfind, Point3d & p1, Point3d & p2);
  virtual void TransformToPlain (INDEX ind, const Point3d & locpoint,
                                 Point2d & plainpoint, double h);
  virtual void TransformFromPlain (INDEX surfind, Point2d & plainpoint,
                                   Point3d & locpoint, double h);
  virtual void ProjectPoint (INDEX /* surfind */, Point3d & /* p */) const { };
  virtual void ProjectPoint2 (INDEX /* surfind */, INDEX /* surfind2 */, Point3d & /* p */) const { };
  virtual void GetNormalVector(INDEX surfind, const Point3d & p, Vec3d & n) const;

  friend double Opti2FunctionValueGrad (const Vector & x, Vector & grad);
  friend double Opti2EdgeFunctionValueGrad (const Vector & x, Vector & grad);
};

#endif
