#include "../meshing/bisect.hpp"
#include "geometry2d.hpp"

namespace meshit {

    void Refinement2d::PointBetween(
            const Point<3>& p1, const Point<3>& p2,
            double secpoint,
            const PointGeomInfo& gi1,
            const PointGeomInfo& gi2,
            Point<3>& newp, PointGeomInfo& newgi) const
    {
        newp = p1 + secpoint * (p2 - p1);
        newgi.trignum = 1;
    }

    void Refinement2d::PointBetween(
            const Point<3>& p1, const Point<3>& p2,
            double secpoint,
            const EdgePointGeomInfo& ap1,
            const EdgePointGeomInfo& ap2,
            Point<3>& newp, EdgePointGeomInfo& newgi) const
    {
        Point<2> p2d;

        p2d = geometry.GetSplines()[ap1.edgenr + 1]->GetPoint(((1 - secpoint) * ap1.dist + secpoint * ap2.dist));

        newp = Point3d(p2d(0), p2d(1), 0);
        newgi.edgenr = ap1.edgenr;
        newgi.dist = ((1 - secpoint) * ap1.dist + secpoint * ap2.dist);
    };

}  // namespace meshit
