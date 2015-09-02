#ifndef FILE_STLTOOL
#define FILE_STLTOOL


#include "../general/array.hpp"
#include "../gprim/adtree.hpp"
#include "../gprim/geom3d.hpp"

#include "stltopology.hpp"

/**************************************************************************/
/* File:   stlgeom.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

namespace netgen {

// use one normal vector for whole chart
extern int usechartnormal;
extern int chartdebug;

extern int geomsearchtreeon;
int AddPointIfNotExists(Array<Point3d>& ap, const Point3d& p, double eps = 1e-8);
//get distance from line lp1-lp2 to point p
double GetDistFromLine(const Point<3>& lp1, const Point<3>& lp2, Point<3>& p);
double GetDistFromInfiniteLine(const Point<3>& lp1, const Point<3>& lp2, const Point<3>& p);


void FIOReadInt(std::istream& ios, int& i);
void FIOWriteInt(std::ostream& ios, const int& i);
void FIOReadDouble(std::istream& ios, double& i);
void FIOWriteDouble(std::ostream& ios, const double& i);
void FIOReadFloat(std::istream& ios, float& i);
void FIOWriteFloat(std::ostream& ios, const float& i);
void FIOReadString(std::istream& ios, char* str, int len);
void FIOReadStringE(std::istream& ios, char* str, int len);
void FIOWriteString(std::ostream& ios, char* str, int len);


typedef Array <int> * ArrayINTPTR;

class STLGeometry;

class STLChart
{
private:
  STLGeometry * geometry;
  Array<int>* charttrigs; // trigs which only belong to this chart
  Array<int>* outertrigs; // trigs which belong to other charts
  Box3dTree * searchtree; // ADT containing outer trigs

  Array<twoint>* olimit; //outer limit of outer chart
  Array<twoint>* ilimit; //outer limit of inner chart


public:
  
  STLChart(STLGeometry * ageometry);
  void AddChartTrig(int i);
  void AddOuterTrig(int i);
  
  int IsInWholeChart(int nr) const;

  int GetChartTrig(int i) const {return charttrigs->Get(i);}
  int GetOuterTrig(int i) const {return outertrigs->Get(i);}
  //get all trigs:
  int GetTrig(int i) const
    {
      if (i <= charttrigs->Size()) {return charttrigs->Get(i);}
      else {return outertrigs->Get(i-charttrigs->Size());}
    }
  
  int GetNChartT() const {return charttrigs->Size();}
  int GetNOuterT() const {return outertrigs->Size();}
  int GetNT() const {return charttrigs->Size()+outertrigs->Size(); }

  void GetTrianglesInBox (const Point3d & pmin,
			  const Point3d & pmax,
			  Array<int> & trias) const;
  void AddOLimit(twoint l) {olimit->Append(l);}
  void AddILimit(twoint l) {ilimit->Append(l);}

  void ClearOLimit() {olimit->SetSize(0);}
  void ClearILimit() {ilimit->SetSize(0);}

  int GetNOLimit() const {return olimit->Size();}
  int GetNILimit() const {return ilimit->Size();}

  twoint GetOLimit(int i) const {return olimit->Get(i);}
  twoint GetILimit(int i) const {return ilimit->Get(i);}

  //move triangles trigs (local chart-trig numbers) to outer chart
  void MoveToOuterChart(const Array<int>& trigs);
  void DelChartTrigs(const Array<int>& trigs);


  // define local coordinate system, JS:
private:
  Vec<3> normal;
  Point<3> pref;
  Vec<3> t1, t2;
public:
  void SetNormal (const Point<3> & apref, const Vec<3> & anormal);
  const Vec<3> & GetNormal () const { return normal; }
  Point<2> Project2d (const Point<3> & p3d) const
  {
    Vec<3> v = p3d-pref;
    return Point<2> (t1 * v, t2 * v);
  }
};

class STLBoundarySeg
{
  Point<3> p1, p2, center;
  Point<2> p2d1, p2d2;
  Box<2> boundingbox;
  //  Point<2> p2dmin, p2dmax;

  double rad;
  int i1, i2;
  int smoothedge;
public:
  STLBoundarySeg () { ; }
  STLBoundarySeg (int ai1, int ai2, const Array<Point<3> > & points,
		  const STLChart * chart)
    : p1(points.Get(ai1)), p2(points.Get(ai2)),
      i1(ai1), i2(ai2)
  {
    center = ::netgen::Center (p1, p2);
    rad = Dist (p1, center);
    
    p2d1 = chart->Project2d (p1);
    p2d2 = chart->Project2d (p2);
    
    boundingbox.Set (p2d1);
    boundingbox.Add (p2d2);
  }

  int operator== (const STLBoundarySeg & s2) const
    { return i1 == s2.i1 && i2 == s2.i2; }
  void Swap ();
  int I1() const { return i1; }
  int I2() const { return i2; }
  const Point<3> & P1() const { return p1; }
  const Point<3> & P2() const { return p2; }
  const Point<2> & P2D1() const { return p2d1; }
  const Point<2> & P2D2() const { return p2d2; }
  const Point<2> & P2DMin() const { return boundingbox.PMin(); }
  const Point<2> & P2DMax() const { return boundingbox.PMax(); }
  const Point<3> & Center() const { return center; }
  const Box<2> & BoundingBox() const { return boundingbox; }
  double Radius () const { return rad; }

  void SetSmoothEdge (int se) { smoothedge = se; }
  int IsSmoothEdge () const { return smoothedge; }
  friend class STLBoundary;
};

class STLBoundary
{
private:
  STLGeometry * geometry;
  const STLChart * chart;
  Array<STLBoundarySeg> boundary;
public:
  STLBoundary(STLGeometry * ageometry);
  // : boundary() {};

  void Clear() {boundary.SetSize(0);};
  void SetChart (const STLChart * achart) { chart = achart; }
  //don't check, if already exists!
  void AddNewSegment(const STLBoundarySeg & seg) {boundary.Append(seg);};
  //check if segment exists
  void AddOrDelSegment(const STLBoundarySeg & seg);
  //addordelsegment for all 3 triangle segments!
  void AddTriangle(const STLTriangle & t);
  int NOSegments() {return boundary.Size();};
  const STLBoundarySeg & GetSegment(int i) {return boundary.Get(i);}

  int TestSeg(const Point<3> & p1, const Point<3> & p2, const Vec<3> & sn, 
	      double sinchartangle, int divisions, Array<Point<3> >& points,
	      double eps);

  int TestSegChartNV(const Point3d& p1, const Point3d& p2, const Vec3d& sn);
};


class STLDoctorParams
{
public:
  int drawmeshededges;
  double geom_tol_fact;

  double longlinefact;
  int showexcluded;

  int selectmode; //0==trig, 1==edge, 2==point, 3==multiedge, 4==line cluster
  int edgeselectmode;

  int useexternaledges;
  int showfaces;
  int showedgecornerpoints;
  int showtouchedtrigchart;
  int conecheck;
  int spiralcheck;
  int selecttrig;
  int nodeofseltrig;
  int selectwithmouse;
  int showmarkedtrigs;
  double dirtytrigfact;
  double smoothangle;

  double smoothnormalsweight;

  int showvicinity;
  int vicinity;
  ///
  STLDoctorParams();
  ///
  void Print (std::ostream & ost) const;
};

extern STLDoctorParams stldoctor;



class STLParameters
{
public:
  /// angle for edge detection
  double yangle;
  double contyangle; //edges continued with contyangle
  /// angle of geometry edge at which the mesher should set a point
  double edgecornerangle;
  /// angle inside on chart
  double chartangle;
  /// angle for overlapping parts of char
  double outerchartangle;
  /// 0 .. no, 1 .. local, (2 .. global)
  int usesearchtree;
  ///
  double resthatlasfac; 
  int resthatlasenable;
  double atlasminh;

  double resthsurfcurvfac; 
  int resthsurfcurvenable;

  double resthchartdistfac;
  int resthchartdistenable;

  double resthcloseedgefac;
  int resthcloseedgeenable;
  
  double resthedgeanglefac;
  int resthedgeangleenable;
  
  double resthsurfmeshcurvfac;
  int resthsurfmeshcurvenable;
  
  double resthlinelengthfac;
  int resthlinelengthenable;

  ///
  int recalc_h_opt;
  ///
  STLParameters();
  ///
  void Print (std::ostream & ost) const;
};

extern STLParameters stlparam;


void STLMeshing (STLGeometry & geom,
		 class Mesh & mesh);


int STLSurfaceMeshing (STLGeometry & geom,
			class Mesh & mesh);

void STLSurfaceOptimization (STLGeometry & geom,
			     class Mesh & mesh,
			     class MeshingParameters & mparam);


}

#endif
