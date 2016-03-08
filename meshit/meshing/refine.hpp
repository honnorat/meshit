#ifndef BISECT_H
#define BISECT_H

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
                double secpoint,
                const PointGeomInfo& gi1,
                const PointGeomInfo& gi2,
                Point<3>& newp, PointGeomInfo& newgi) const;

        virtual void PointBetween(
                const Point<3>& p1, const Point<3>& p2,
                double secpoint,
                const EdgePointGeomInfo& ap1,
                const EdgePointGeomInfo& ap2,
                Point<3>& newp, EdgePointGeomInfo& newgi) const;
    };

}  // namespace meshit

#endif
