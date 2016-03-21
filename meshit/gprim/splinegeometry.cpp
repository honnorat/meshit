/*
 * 2d Spline curve for Mesh generator
 */

#include "splinegeometry.hpp"

namespace meshit {

    SplineGeometry::~SplineGeometry()
    {
        for (size_t i = 0; i < splines.size(); i++) {
            delete splines[i];
        }
    }

    void SplineGeometry::GetBoundingBox(Box<2>& box) const
    {
        if (!splines.size()) {
            Point<2> auxp = 0.;
            box.Set(auxp);
            return;
        }

        std::vector<Point<2> > points;
        for (size_t i = 0; i < splines.size(); i++) {
            splines[i]->GetPoints(20, points);

            if (i == 0) box.Set(points[0]);
            for (size_t j = 0; j < points.size(); j++) {
                box.Add(points[j]);
            }
        }
    }
}


