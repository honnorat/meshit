#ifndef FILE_GEOMTEST3D_HPP
#define FILE_GEOMTEST3D_HPP

/* *************************************************************************/
/* File:   geomtest3d.hh                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   13. Feb. 98                                                     */
/* *************************************************************************/

#include "geomobjects.hpp"
#include "geom3d.hpp"

namespace meshit
{
    int IntersectTriangleLine(const Point2d** tri, const Point2d** line);

    // 1, iff not regular triangulation
    int IntersectTriangleTriangle(const Point2d** tri1, const Point2d** tri2);

    void LocalCoordinates(const Vec3d& e1, const Vec3d& e2, const Vec3d& v, double& lam1, double& lam2);

    /// return 1 = degenerated sphere
    int CalcSphereCenter(const Point3d** pts, Point3d& c);

    /// Minimal distance of point p to the line segment [lp1,lp2]
    double MinDistLP2(const Point2d& lp1, const Point2d& lp2, const Point2d& p);

    /// Minimal distance of point p to the line segment [lp1,lp2]
    double MinDistLP2(const Point3d& lp1, const Point3d& lp2, const Point3d& p);

    /// Minimal distance of point p to the triangle segment [tp1,tp2,pt3]
    double MinDistTP2(const Point3d& tp1, const Point3d& tp2, const Point3d& tp3, const Point3d& p);

}  // namespace meshit

#endif
