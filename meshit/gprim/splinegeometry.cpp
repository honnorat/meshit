/*
 * 2d Spline curve for Mesh generator
 */

#include "splinegeometry.hpp"

namespace meshit {

    template<int D>
    SplineGeometry<D>::~SplineGeometry()
    {
        for (size_t i = 0; i < splines.size(); i++) {
            delete splines[i];
        }
    }

    template<int D>
    void SplineGeometry<D>::GetBoundingBox(Box<D>& box) const
    {
        if (!splines.size()) {
            Point<D> auxp = 0.;
            box.Set(auxp);
            return;
        }

        Array<Point<D> > points;
        for (size_t i = 0; i < splines.size(); i++) {
            splines[i]->GetPoints(20, points);

            if (i == 0) box.Set(points[0]);
            for (size_t j = 0; j < points.size(); j++) {
                box.Add(points[j]);
            }
        }
    }

    template
    class SplineGeometry<2>;

    template
    class SplineGeometry<3>;
}


