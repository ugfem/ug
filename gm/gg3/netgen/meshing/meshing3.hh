// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_MESHING3
#define FILE_MESHING3

#include <meshing/meshtype.hh>
#include <meshing/adfront3.hh>
#include <meshing/ruler3.hh>

class Meshing3
{
protected:
  ADFRONT3 * adfront;
  ARRAY<vnetrule*> rules;
  ARRAY<int> ruleused;
  ARRAY<char*> problems;
  double tolfak;

public:
  Meshing3 (char * rulefilename);
  virtual ~Meshing3 ();

  void LoadRules (char * filename);
  void Mesh (double gh,int prism);

  void ImproveMesh (ARRAY<Point3d> & points,
                    const ARRAY<Element> & surfelements,
                    const ARRAY<Element> & elements,
                    double h);

  void ImproveMesh (ARRAY<Point3d> & points,
                    const ARRAY<Element> & elements,
                    int nboundnodes,
                    double h);

  void CombineImprove (ARRAY<Point3d> & points,
                       ARRAY<Element> & elements,
                       int nboundnodes,
                       double h);

  void SplitImprove (ARRAY<Point3d> & points,
                     ARRAY<Element> & elements, ARRAY<Element> & surfelements,
                     double h);

  void SwapImprove (ARRAY<Point3d> & points,
                    ARRAY<Element> & elements, ARRAY<Element> & surfelements,
                    double h);


  void AddPoint (const Point3d & p, INDEX globind);
  void AddBoundaryElement (const Element & elem, int inverse, int prism_flag);

  virtual int SavePoint (const Point3d & /* p */) { return 0; }
  virtual void SaveElement (const Element & /* elem */) { };
  virtual void Get_Local_h_3d(double *in,double *out) { };

  friend void PlotVolMesh (const ROT3D & r, char key);
  friend void TestRules ();

};
#endif
