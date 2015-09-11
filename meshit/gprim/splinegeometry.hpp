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

    template < int D >
    class SplineGeometry
    {
        // protected:
      public:
        std::vector<GeomPoint<D> > geompoints;
        std::vector<SplineSeg<D>*> splines;

        ~SplineGeometry();

        int Load(const Array<double> & raw_data, const int startpos = 0);
        void GetRawData(Array<double> & raw_data) const;

        const std::vector<SplineSeg<D>*> & GetSplines() const
        {
            return splines;
        }

        int GetNSplines(void) const
        {
            return splines.size();
        }

        std::string GetSplineType(const int i) const
        {
            return splines[i]->GetType();
        }

        SplineSeg<D> & GetSpline(const int i)
        {
            return *splines[i];
        }

        const SplineSeg<D> & GetSpline(const int i) const
        {
            return *splines[i];
        }

        void GetBoundingBox(Box<D> & box) const;

        Box<D> GetBoundingBox() const
        {
            Box<D> box;
            GetBoundingBox(box);
            return box;
        }

        int GetNP() const
        {
            return geompoints.size();
        }

        const GeomPoint<D> & GetPoint(int i) const
        {
            return geompoints[i];
        }
    };
}

#endif // _FILE_SPLINEGEOMETRY
