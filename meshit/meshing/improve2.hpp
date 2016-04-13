#ifndef IMPROVE2_HPP
#define IMPROVE2_HPP

#include "meshclass.hpp"

namespace meshit
{
    class MeshOptimize2d
    {
     public:
        MeshOptimize2d(Mesh& mesh) : mesh_{mesh} { }

        ~MeshOptimize2d() { }

        void EdgeSwapping(bool use_metric, double metric_weight = 0.0);
        void ImproveMesh(double metric_weight);
        void CombineImprove(double metric_weight);

        friend class Opti2SurfaceMinFunction;

     protected:
        void EdgeSwapping(size_t face_index, bool use_metric, double metric_weight);
        void ImproveMesh(size_t face_index, double metric_weight);
        void CombineImprove(size_t face_index, double metric_weight);

     protected:
        Mesh& mesh_;
    };

    double CalcTriangleBadness(
        const Point2d& p1,
        const Point2d& p2,
        const Point2d& p3,
        double metric_weight,
        double h);

    double CalcTriangleBadness_2(
        const Point2d& p1,
        const Point2d& p2,
        const Point2d& p3,
        double n_z);

}  // namespace meshit

#endif

