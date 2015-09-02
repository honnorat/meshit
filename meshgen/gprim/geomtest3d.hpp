#ifndef FILE_GEOMTEST3D
#define FILE_GEOMTEST3D

/* *************************************************************************/
/* File:   geomtest3d.hh                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   13. Feb. 98                                                     */
/* *************************************************************************/

#include "geomobjects.hpp"
#include "geom3d.hpp"

namespace netgen
{


int IntersectTriangleLine (const Point<3> ** tri, const Point<3> ** line);



/**
  Returns 0, iff
  closure (tet)  cup  closure (tri)  is empty, one corner point of tet,
  one edge of tet or one face of tet
 */
int 
IntersectTetTriangle (const Point<3> ** tet, const Point<3> ** tri,
		      const int * tetpi = NULL, const int * tripi = NULL);

/**
  Same test as above, but tet int reference position (0, ex, ey, ez),
  tetpi = 1, 2, 4, 5
 */
int 
IntersectTetTriangleRef (const Point3d ** tri, const int * tripi = NULL);


// 1, iff not regular triangulation
int 
IntersectTriangleTriangle (const Point<3> ** tri1, const Point<3> ** tri2);


void
LocalCoordinates (const Vec3d & e1, const Vec3d & e2,
		  const Vec3d & v, double & lam1, double & lam2);

/// return 1 = degenerated sphere
int
CalcSphereCenter (const Point<3> ** pts, Point<3> & c);

/// return 1 = degenerated triangle
int
CalcTriangleCenter (const Point3d ** pts, Point3d & c);



/*
  Compute radius of cylinder fitting 4 points.
  cylinder axis is in the direction of p1-p2
*/
double ComputeCylinderRadius (const Point3d & p1, const Point3d & p2,
				     const Point3d & p3, const Point3d & p4);

/*
  Two triangles T1 and T2 have normals n1 and n2.
  The height over the common edge is h1, and h2.
  Radius of cylinder fitting both triangles
*/
double ComputeCylinderRadius (const Vec3d & n1, const Vec3d & n2,
				     double h1, double h2);

/// Minimal distance of point p to the line segment [lp1,lp2]
double MinDistLP2 (const Point2d & lp1, const Point2d & lp2, const Point2d & p);

/// Minimal distance of point p to the line segment [lp1,lp2]
double MinDistLP2 (const Point3d & lp1, const Point3d & lp2, const Point3d & p);

/// Minimal distance of point p to the triangle segment [tp1,tp2,pt3]
double MinDistTP2 (const Point3d & tp1, const Point3d & tp2, 
			  const Point3d & tp3, const Point3d & p);

/// Minimal distance of the 2 lines [l1p1,l1p2] and [l2p1,l2p2]
double MinDistLL2 (const Point3d & l1p1, const Point3d & l1p2,
			  const Point3d & l2p1, const Point3d & l2p2);

}

#endif
