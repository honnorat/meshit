#ifndef FILE_GEOMTEST3D
#define FILE_GEOMTEST3D

/* *************************************************************************/
/* File:   geomtest3d.hh                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   13. Feb. 98                                                     */
/* *************************************************************************/

#include "geomobjects.hpp"
#include "geom3d.hpp"

namespace meshit {


    int IntersectTriangleLine(const Point3d ** tri, const Point3d ** line);

    // 1, iff not regular triangulation
    int IntersectTriangleTriangle(const Point3d ** tri1, const Point3d ** tri2);

    void LocalCoordinates(const Vec3d & e1, const Vec3d & e2, const Vec3d & v, double & lam1, double & lam2);

    /// return 1 = degenerated sphere
    int CalcSphereCenter(const Point3d ** pts, Point3d & c);

    /// return 1 = degenerated triangle
    int CalcTriangleCenter(const Point3d ** pts, Point3d & c);

    /*
      Compute radius of cylinder fitting 4 points.
      cylinder axis is in the direction of p1-p2
     */
    double ComputeCylinderRadius(const Point3d & p1, const Point3d & p2, const Point3d & p3, const Point3d & p4);

    /*
      Two triangles T1 and T2 have normals n1 and n2.
      The height over the common edge is h1, and h2.
      Radius of cylinder fitting both triangles
     */
    double ComputeCylinderRadius(const Vec3d & n1, const Vec3d & n2, double h1, double h2);

    /// Minimal distance of point p to the line segment [lp1,lp2]
    double MinDistLP2(const Point2d & lp1, const Point2d & lp2, const Point2d & p);

    /// Minimal distance of point p to the line segment [lp1,lp2]
    double MinDistLP2(const Point3d & lp1, const Point3d & lp2, const Point3d & p);

    /// Minimal distance of point p to the triangle segment [tp1,tp2,pt3]
    double MinDistTP2(const Point3d & tp1, const Point3d & tp2, const Point3d & tp3, const Point3d & p);

    /// Minimal distance of the 2 lines [l1p1,l1p2] and [l2p1,l2p2]
    double MinDistLL2(const Point3d & l1p1, const Point3d & l1p2, const Point3d & l2p1, const Point3d & l2p2);

}

#endif
