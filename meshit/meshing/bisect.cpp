#include <sstream>

#include "bisect.hpp"

namespace meshit {

    Refinement::Refinement() { }

    void Refinement::PointBetween(
            const Point<3>& p1, const Point<3>& p2, double secpoint,
            const PointGeomInfo& gi1,
            const PointGeomInfo& gi2,
            Point<3>& newp, PointGeomInfo& newgi) const
    {
        newp = p1 + secpoint * (p2 - p1);
    }

    void Refinement::PointBetween(const Point<3>& p1, const Point<3>& p2, double secpoint,
                                  const EdgePointGeomInfo& ap1,
                                  const EdgePointGeomInfo& ap2,
                                  Point<3>& newp, EdgePointGeomInfo& newgi) const
    {
        newp = p1 + secpoint * (p2 - p1);
    }
}
