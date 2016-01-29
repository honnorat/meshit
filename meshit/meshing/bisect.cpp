#include <sstream>
#include <stdexcept>

#include "bisect.hpp"
#include "../meshit.hpp"
#include "../gprim/geomfuncs.hpp"

namespace meshit {

    Refinement::Refinement()
    {
        optimizer2d = NULL;
    }

    void Refinement::PointBetween(
            const Point<3> & p1, const Point<3> & p2, double secpoint,
            int surfi,
            const PointGeomInfo & gi1,
            const PointGeomInfo & gi2,
            Point<3> & newp, PointGeomInfo & newgi) const
    {
        newp = p1 + secpoint * (p2 - p1);
    }

    void Refinement::PointBetween(const Point<3> & p1, const Point<3> & p2, double secpoint,
            int surfi1, int surfi2,
            const EdgePointGeomInfo & ap1,
            const EdgePointGeomInfo & ap2,
            Point<3> & newp, EdgePointGeomInfo & newgi) const
    {
        std::cout << "base class edge point between" << std::endl;
        newp = p1 + secpoint * (p2 - p1);
    }

    Vec<3> Refinement::GetTangent(const Point<3> & p, int surfi1, int surfi2,
            const EdgePointGeomInfo & ap1) const
    {
        std::cerr << "Refinement::GetTangent not overloaded" << std::endl;
        return Vec<3> (0, 0, 0);
    }

    Vec<3> Refinement::GetNormal(const Point<3> & p, int surfi1,
            const PointGeomInfo & gi) const
    {
        std::cerr << "Refinement::GetNormal not overloaded" << std::endl;
        return Vec<3> (0, 0, 0);
    }

    void Refinement::ProjectToSurface(Point<3> & p, int surfi) const
    {
        LOG_ERROR("Refinement::ProjectToSurface: no geometry set");
    };

    void Refinement::ProjectToEdge(Point<3> & p, int surfi1, int surfi2, const EdgePointGeomInfo & egi) const
    {
        LOG_ERROR("Refinement::ProjectToEdge not overloaded");
    }
}
