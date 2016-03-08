#ifndef FILE_IMPROVE2
#define FILE_IMPROVE2

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

        /// project point, use gi as initial value, and compute new gi
        int ProjectPointGI(INDEX surfind, Point3d& p, PointGeomInfo& gi) const
        {
            return CalcPointGeomInfo(surfind, gi, p);
        }

        /// liefert zu einem 3d-Punkt die geominfo (Dreieck) und liefert 1, wenn erfolgreich, 
        /// 0, wenn nicht (Punkt ausserhalb von chart)
        int CalcPointGeomInfo(PointGeomInfo& gi, const Point3d& /*p3*/) const
        {
            gi.trignum = 1;
            return 1;
        };

        int CalcPointGeomInfo(int /* surfind */, PointGeomInfo& gi, const Point3d& p3) const
        {
            return CalcPointGeomInfo(gi, p3);
        }

        void GetNormalVector(INDEX surfind, const Point3d& p, PointGeomInfo& gi, Vec3d& n) const;
        void GetNormalVector(INDEX surfind, const Point3d& p, Vec3d& n) const;

        friend class Opti2SurfaceMinFunction;
    };

    void CalcTriangleBadness(
            double x2, double x3, double y3,
            double metricweight,
            double h, double& badness,
            double& g1x, double& g1y);

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
}

#endif

