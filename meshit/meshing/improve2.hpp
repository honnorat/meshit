#ifndef IMPROVE2_HPP
#define IMPROVE2_HPP

#include "meshclass.hpp"

namespace meshit
{
    class MeshOptimize2d
    {
     public:
        MeshOptimize2d()
            : faceindex{0}, metricweight{0.0} { }

        ~MeshOptimize2d() { }

        void ImproveMesh(Mesh& mesh2d, const MeshingParameters& mp);
        void EdgeSwapping(Mesh& mesh, int usemetric);
        void CombineImprove(Mesh& mesh);

        void SetMetricWeight(double mw)
        {
            metricweight = mw;
        }

        friend class Opti2SurfaceMinFunction;

     protected:
        size_t faceindex;
        double metricweight;
    };

    double CalcTriangleBadness(
        const Point3d& p1,
        const Point3d& p2,
        const Point3d& p3,
        double metricweight,
        double h);

    double CalcTriangleBadness_2(
        const Point3d& p1,
        const Point3d& p2,
        const Point3d& p3,
        double n_z);

}  // namespace meshit

#endif

