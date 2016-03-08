#ifndef FILE_GEOM2DMESH_H
#define FILE_GEOM2DMESH_H

/**************************************************************************/
/* File:   geom2dmesh.hh                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   22. Jan. 01                                                    */
/**************************************************************************/

#include "../gprim/geomobjects.hpp"
#include "../meshing/refine.hpp"

namespace meshit {

    class Refinement2d : public Refinement
    {
        const class SplineGeometry2d& geometry;

     public:
        explicit Refinement2d(const class SplineGeometry2d& ageometry)
                : Refinement(), geometry(ageometry) { }

        virtual ~Refinement2d() { }

        virtual void PointBetween(
                const Point<3>& p1, const Point<3>& p2,
                double secpoint,
                const PointGeomInfo& gi1,
                const PointGeomInfo& gi2,
                Point<3>& newp, PointGeomInfo& newgi) const;

        virtual void PointBetween(
                const Point<3>& p1, const Point<3>& p2,
                double secpoint,
                const EdgePointGeomInfo& ap1,
                const EdgePointGeomInfo& ap2,
                Point<3>& newp, EdgePointGeomInfo& newgi) const;
    };

}  // namespace meshit

#endif
