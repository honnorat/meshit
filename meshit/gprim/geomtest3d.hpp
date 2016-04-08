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
    bool IntersectTriangleLine(const Point2d** tri, const Point2d** line);

    // 1, iff not regular triangulation
    bool IntersectTriangleTriangle(const Point2d** tri1, const Point2d** tri2);

}  // namespace meshit

#endif
