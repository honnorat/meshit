#ifndef FILE_GEOMTEST3D_HPP
#define FILE_GEOMTEST3D_HPP

/* *************************************************************************/
/* File:   geomtest3d.hh                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   13. Feb. 98                                                     */
/* *************************************************************************/

#include "geom3d.hpp"
#include "geomobjects.hpp"

namespace meshit {

bool IntersectTriangleLine(const Point2d** tri, const Point2d** line);
bool IntersectTriangleTriangle(const Point2d** tri1, const Point2d** tri2);

}  // namespace meshit

#endif
