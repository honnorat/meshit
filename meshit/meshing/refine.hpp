#ifndef REFINE_HPP
#define REFINE_HPP

#include "improve2.hpp"

namespace meshit {

    class Refinement
    {
     public:
        Refinement() { }

        virtual ~Refinement() { }

        void Refine(Mesh& mesh);

        virtual void PointBetween(
                const Point<3>& p1, const Point<3>& p2,
                double secpoint, Point<3>& newp) const;

        virtual void PointBetween(
                const Point<3>& p1, const Point<3>& p2,
                double secpoint,
                const EdgePointGeomInfo& ap1,
                const EdgePointGeomInfo& ap2,
                Point<3>& newp, EdgePointGeomInfo& newgi) const;
    };

}  // namespace meshit

#endif
