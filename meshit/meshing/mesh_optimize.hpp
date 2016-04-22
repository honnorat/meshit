#ifndef MESH_OPTIMIZE_HPP
#define MESH_OPTIMIZE_HPP

#include "meshclass.hpp"

namespace meshit {

class MeshOptimize
{
 public:
    explicit MeshOptimize(Mesh& mesh)
        : mesh_{mesh} { }

    ~MeshOptimize() { }

    void EdgeSwapping(bool use_metric, double metric_weight = 0.0);
    void ImproveMesh(double metric_weight);
    void CombineImprove();

    friend class Opti2SurfaceMinFunction;

 protected:
    void EdgeSwapping(size_t face_index, bool use_metric, double metric_weight);
    void ImproveMesh(size_t face_index, double metric_weight);
    void CombineImprove(size_t face_index);

 protected:
    Mesh& mesh_;
};

double CalcTriangleBadness(const Point2d& p1, const Point2d& p2, const Point2d& p3,
                           double metric_weight, double h);

double CalcTriangleBadness_2(const Point2d& p1, const Point2d& p2, const Point2d& p3, double n_z);

}  // namespace meshit

#endif
