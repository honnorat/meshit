#ifndef FILE_TRIAPPROX
#define FILE_TRIAPPROX

/**************************************************************************/
/* File:   triapprox.hh                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   2. Mar. 98                                                    */
/**************************************************************************/


namespace netgen
{

  /**
     Triangulated approxiamtion to true surface
  */
 

  class TATriangle
  {
    int pi[3];
    int surfind;
  public:
    TATriangle () { ; }

    TATriangle (int si, int pi1, int pi2, int pi3)
    { surfind = si; pi[0] = pi1; pi[1] = pi2; pi[2] = pi3; }

    int SurfaceIndex() const { return surfind; }
    int & SurfaceIndex() { return surfind; }

    int & operator[] (int i) { return pi[i]; }
    const int & operator[] (int i) const { return pi[i]; }
  };


  class TriangleApproximation
  {
    Array<Point<3> > points;
    Array<Vec<3> > normals;
    Array<TATriangle> trigs;

  public:
    TriangleApproximation();
    int GetNP () const { return points.size(); }
    int GetNT () const { return trigs.size(); }

    int AddPoint (const Point<3> & p) { points.push_back (p); return points.size()-1; }
    int AddNormal (const Vec<3> & n) { normals.push_back (n); return normals.size()-1; }
    int AddTriangle (const TATriangle & tri, bool invert = 0);

    const Point<3> & GetPoint (int i) const { return points[i]; }
    const TATriangle & GetTriangle (int i) const { return trigs[i]; }
    const Vec<3> & GetNormal (int i) const { return normals[i]; }

    void RemoveUnusedPoints ();

    friend class CSGeometry;
  };

}

#endif
