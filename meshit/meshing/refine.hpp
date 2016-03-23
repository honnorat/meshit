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

        void PointBetween(const Point3d& p1, const Point3d& p2, double secpoint, Point3d& newp) const;
    };

}  // namespace meshit

#endif
