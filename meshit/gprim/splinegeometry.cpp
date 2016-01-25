/*

2d Spline curve for Mesh generator

 */

#include <stdexcept>
#include <meshit/meshit.hpp>
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
    void SplineGeometry<D>::GetRawData(Array<double> & raw_data) const
    {
        raw_data.push_back(D);
        raw_data.push_back(splines.size());
        for (size_t i = 0; i < splines.size(); i++) {
            splines[i]->GetRawData(raw_data);
        }
    }

    template<int D>
    int SplineGeometry<D>::Load(const Array<double> & raw_data, const int startpos)
    {
        int pos = startpos;
        if (raw_data[pos] != D)
            throw std::runtime_error("wrong dimension of spline raw_data");

        pos++;

        splines.resize(int(raw_data[pos]));
        pos++;

        Array< Point<D> > pts(3);

        for (size_t i = 0; i < splines.size(); i++) {
            int type = int(raw_data[pos]);
            pos++;

            for (int j = 0; j < type; j++) {
                for (int k = 0; k < D; k++) {
                    pts[j](k) = raw_data[pos];
                    pos++;
                }
            }
            if (type == 2) {
                splines[i] = new LineSeg<D>(
                        GeomPoint<D>(pts[0], 1),
                        GeomPoint<D>(pts[1], 1));
            }
            else if (type == 3) {
                splines[i] = new SplineSeg3<D>(
                        GeomPoint<D>(pts[0], 1),
                        GeomPoint<D>(pts[1], 1),
                        GeomPoint<D>(pts[2], 1));
            }
            else
                throw std::runtime_error("something wrong with spline raw data");
        }
        return pos;
    }

    template<int D>
    void SplineGeometry<D>::GetBoundingBox(Box<D> & box) const
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
            for (int j = 0; j < points.size(); j++) {
                box.Add(points[j]);
            }
        }
    }

    template class SplineGeometry<2>;
    template class SplineGeometry<3>;
}


