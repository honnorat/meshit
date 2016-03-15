#ifndef IMPROVE2_HPP
#define IMPROVE2_HPP

#include "meshclass.hpp"

namespace meshit {

    class MeshOptimize2d
    {
        int faceindex;
        double metricweight;
        int writestatus;

     public:
        MeshOptimize2d();

        virtual ~MeshOptimize2d() { }

        void ImproveMesh(Mesh& mesh2d, const MeshingParameters& mp);

        void EdgeSwapping(Mesh& mesh, int usemetric);
        void CombineImprove(Mesh& mesh);

        void GenericImprove(Mesh& mesh);

        void SetFaceIndex(int fi)
        {
            faceindex = fi;
        }

        void SetImproveEdges(int ie) { }

        void SetMetricWeight(double mw)
        {
            metricweight = mw;
        }

        void SetWriteStatus(int ws)
        {
            writestatus = ws;
        }

        friend class Opti2SurfaceMinFunction;
    };

    double CalcTriangleBadness(
            const Point3d& p1,
            const Point3d& p2,
            const Point3d& p3,
            double metricweight,
            double h);

    double CalcTriangleBadness(
            const Point3d& p1,
            const Point3d& p2,
            const Point3d& p3,
            const Vec3d& n,
            double metricweight,
            double h);
}  // namespace meshit

#endif

