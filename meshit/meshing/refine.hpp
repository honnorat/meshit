#ifndef REFINE_HPP
#define REFINE_HPP

#include "improve2.hpp"

namespace meshit
{
    class Refinement
    {
     public:
        Refinement() { }

        ~Refinement() { }

        void Refine(Mesh& mesh);

        void PointBetween(const MeshPoint& p1, const MeshPoint& p2, double secpoint, Point2d& newp) const;
    };

}  // namespace meshit

#endif
