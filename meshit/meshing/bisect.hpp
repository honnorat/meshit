#ifndef BISECT
#define BISECT

#include "improve2.hpp"

namespace meshit {

    class Refinement
    {

     public:

        Refinement();

        virtual ~Refinement() { }

        void Refine(Mesh& mesh) const;
        void Refine(Mesh& mesh);

        virtual void PointBetween(
                const Point<3>& p1, const Point<3>& p2,
                double secpoint,
                int surfi,
                const PointGeomInfo& gi1,
                const PointGeomInfo& gi2,
                Point<3>& newp, PointGeomInfo& newgi) const;

        virtual void PointBetween(
                const Point<3>& p1, const Point<3>& p2,
                double secpoint,
                int surfi1, int surfi2,
                const EdgePointGeomInfo& ap1,
                const EdgePointGeomInfo& ap2,
                Point<3>& newp, EdgePointGeomInfo& newgi) const;

        virtual Vec<3> GetTangent(
                const Point<3>& p, int surfi1, int surfi2,
                const EdgePointGeomInfo& egi) const;

        virtual Vec<3> GetNormal(
                const Point<3>& p, int surfi1,
                const PointGeomInfo& gi) const;

        virtual void ProjectToSurface(Point<3>& p, int surfi) const;

        virtual void ProjectToEdge(Point<3>& p, int surfi1, int surfi2, const EdgePointGeomInfo& egi) const;

    };
}
#endif
