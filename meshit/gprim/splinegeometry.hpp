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

    template<int D>
    class SplineGeometry
    {
     public:
        std::vector<GeomPoint<D> > geompoints;
        std::vector<SplineSeg<D>*> splines;

        virtual ~SplineGeometry();

        const std::vector<SplineSeg<D>*>& GetSplines() const
        {
            return splines;
        }
        void GetBoundingBox(Box<D>& box) const;

        Box<D> GetBoundingBox() const
        {
            Box<D> box;
            GetBoundingBox(box);
            return box;
        }

        size_t GetNP() const
        {
            return geompoints.size();
        }

        const GeomPoint<D>& GetPoint(size_t i) const
        {
            return geompoints[i];
        }
    };
}

#endif // _FILE_SPLINEGEOMETRY
