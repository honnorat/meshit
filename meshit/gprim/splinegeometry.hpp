#ifndef _FILE_SPLINEGEOMETRY
#define _FILE_SPLINEGEOMETRY

/*

JS, Nov 2007

The 2D/3D template-base classes should go into the libsrc/gprim directory

in geom2d only 2D - Geometry classes (with material properties etc.)
 */

#include <vector>
#include <string>
#include "spline.hpp"

namespace meshit {

    class SplineGeometry
    {
     public:
        std::vector<GeomPoint<2> > geompoints;
        std::vector<SplineSeg*> splines;

        virtual ~SplineGeometry();

        const std::vector<SplineSeg*>& GetSplines() const
        {
            return splines;
        }
        void GetBoundingBox(Box<2>& box) const;

        Box<2> GetBoundingBox() const
        {
            Box<2> box;
            GetBoundingBox(box);
            return box;
        }

        size_t GetNP() const
        {
            return geompoints.size();
        }

        const GeomPoint<2>& GetPoint(size_t i) const
        {
            return geompoints[i];
        }
    };
}

#endif // _FILE_SPLINEGEOMETRY
