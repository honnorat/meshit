#include "improve2.hpp"
#include "global.hpp"

#include "../linalg/opti.hpp"

namespace meshit {

    static const double c_trig = 0.1443375672974;  // sqrt(3.0) / 12
    static const double c_trig4 = 0.577350269189;  // sqrt(3.0) / 3

    inline double CalcTriangleBadness(double x2, double x3, double y3, double metricweight, double h)
    {
        // badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1
        // p1 = (0, 0), p2 = (x2, 0), p3 = (x3, y3);

        double cir_2 = (x2 * x2 + x3 * x3 + y3 * y3 - x2 * x3);
        double area = x2 * y3;

        if (area <= 1e-24 * cir_2)
            return 1e10;

        double badness = c_trig4 * cir_2 / area - 1;

        if (metricweight > 0) {
            // add:  metricweight * (area / h^2 + h^2 / area - 2)

            double areahh = area / (h * h);
            badness += metricweight * (areahh + 1 / areahh - 2);
        }
        return badness;
    }

    double CalcTriangleBadness(
            const Point3d& p1,
            const Point3d& p2,
            const Point3d& p3,
            double metricweight,
            double h)
    {
        // badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1

        Vec3d e12 = p2 - p1;
        Vec3d e13 = p3 - p1;
        Vec3d e23 = p3 - p2;

        double cir_2 = e12.Length2() + e13.Length2() + e23.Length2();
        double area = 0.5 * Cross(e12, e13).Length();

        if (area <= 1e-24 * cir_2)
            return 1e10;

        double badness = c_trig * cir_2 / area - 1;

        if (metricweight > 0) {
            // add:  metricweight * (area / h^2 + h^2 / area - 2)
            area *= 2;   // optimum for (2 area) is h^2
            double areahh = area / (h * h);
            badness += metricweight * (areahh + 1 / areahh - 2);
        }

        return badness;
    }

    double CalcTriangleBadnessGrad(
            const Point3d& p1,
            const Point3d& p2,
            const Point3d& p3,
            Vec3d& gradp1,
            double metricweight,
            double h)
    {
        // badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1

        Vec3d e12 = p2 - p1;
        Vec3d e13 = p3 - p1;
        Vec3d e23 = p3 - p2;

        double cir_2 = e12.Length2() + e13.Length2() + e23.Length2();
        Vec3d varea = Cross(e12, e13);
        double area = 0.5 * varea.Length();

        Vec3d dcir_2 = (-2) * (e12 + e13);
        Vec3d darea = (0.25 / area) * Cross(p2 - p3, varea);

        if (area <= 1e-24 * cir_2) {
            gradp1 = 0;
            return 1e10;
        }

        double badness = c_trig * cir_2 / area - 1;
        gradp1 = c_trig * (1.0 / area * dcir_2 - cir_2 / (area * area) * darea);

        if (metricweight > 0) {
            // add:  metricweight * (area / h^2 + h^2 / area - 2)
            area *= 2;  // optimum for (2 area) is h^2

            double areahh = area / (h * h);
            badness += metricweight * (areahh + 1 / areahh - 2);

            gradp1 += (2 * metricweight * (1 / (h * h) - (h * h) / (area * area))) * darea;
        }

        return badness;
    }

    double CalcTriangleBadness(
            const Point3d& p1,
            const Point3d& p2,
            const Point3d& p3,
            const Vec3d& n,
            double metricweight,
            double h)
    {
        Vec3d v1 = p2 - p1;
        Vec3d v2 = p3 - p1;

        Vec3d e1 = v1;
        Vec3d e2 = v2;

        e1 -= (e1 * n) * n;
        e1 /= (e1.Length() + 1e-24);
        e2 = Cross(n, e1);

        return CalcTriangleBadness((e1 * v1), (e1 * v2), (e2 * v2), metricweight, h);
    }

    class Opti2dLocalData
    {
     public:
        MeshPoint sp1;
        Vec3d normal, t1, t2;
        Array<SurfaceElementIndex> locelements;
        Array<int> locrots;
        Array<double> lochs;
        Array<Point3d> loc_pnts2, loc_pnts3;
        double locmetricweight;
        double loch;
        int uselocalh;

     public:
        Opti2dLocalData()
        {
            locmetricweight = 0;
        }
    };

    class Opti2SurfaceMinFunction : public MinFunction
    {
        const Mesh& mesh;
        Opti2dLocalData& ld;

     public:
        Opti2SurfaceMinFunction(const Mesh& amesh, Opti2dLocalData& ald)
                : mesh(amesh), ld(ald) { }

        virtual double Func(const Vector& x) const
        {
            double badness = 0;

            Point3d pp1(ld.sp1.X() - x(1),
                        ld.sp1.Y() + x(0));

            for (size_t j = 0; j < ld.locelements.size(); j++) {
                const Point3d& pe2 = ld.loc_pnts2[j];
                const Point3d& pe3 = ld.loc_pnts3[j];
                Point3d e2(pe2.X() - pp1.X(), pe2.Y() - pp1.Y());
                Point3d e3(pe3.X() - pp1.X(), pe3.Y() - pp1.Y());

                if (ld.uselocalh) ld.loch = ld.lochs[j];

                if (e2.X() * e3.Y() - e2.Y() * e3.X() > 1e-8 * ld.loch * ld.loch) {
                    badness += CalcTriangleBadness(pp1, pe2, pe3, ld.locmetricweight, ld.loch);
                } else {
                    badness += 1e8;
                }
            }
            return badness;
        }

        virtual double FuncGrad(const Vector& x, Vector& g) const
        {

            double badness = 0;

            Vec3d vgrad;
            Point3d pp1(ld.sp1.X() - x(1),
                        ld.sp1.Y() + x(0));

            for (size_t j = 0; j < ld.locelements.size(); j++) {
                const Point3d& pe2 = ld.loc_pnts2[j];
                const Point3d& pe3 = ld.loc_pnts3[j];
                Point3d e2(pe2.X() - pp1.X(), pe2.Y() - pp1.Y());
                Point3d e3(pe3.X() - pp1.X(), pe3.Y() - pp1.Y());

                if (ld.uselocalh) ld.loch = ld.lochs[j];

                if (e2.X() * e3.Y() - e2.Y() * e3.X() > 1e-8 * ld.loch * ld.loch) {
                    Vec3d hgrad;
                    badness += CalcTriangleBadnessGrad(pp1, pe2, pe3, hgrad, ld.locmetricweight, ld.loch);
                    vgrad += hgrad;
                } else {
                    badness += 1e8;
                }
            }
            g[0] = +vgrad.Y();
            g[1] = -vgrad.X();
            return badness;
        }

        virtual double FuncDeriv(const Vector& x, const Vector& dir, double& deriv) const
        {
            deriv = 0;
            double badness = 0;

            Point3d pp1(ld.sp1.X() - x[1],
                        ld.sp1.Y() + x[0]);

            for (size_t j = 0; j < ld.locelements.size(); j++) {
                const Point3d& pe2 = ld.loc_pnts2[j];
                const Point3d& pe3 = ld.loc_pnts3[j];
                Point3d e2(pe2.X() - pp1.X(), pe2.Y() - pp1.Y());
                Point3d e3(pe3.X() - pp1.X(), pe3.Y() - pp1.Y());

                if (ld.uselocalh) ld.loch = ld.lochs[j];

                if (e2.X() * e3.Y() - e2.Y() * e3.X() > 1e-8 * ld.loch * ld.loch) {
                    Vec3d hgrad;
                    badness += CalcTriangleBadnessGrad(pp1, pe2, pe3, hgrad, ld.locmetricweight, ld.loch);
                    deriv += dir[0] * hgrad.Y() - dir[1] * hgrad.X();
                } else {
                    badness += 1e8;
                }
            }
            return badness;
        }
    };

    MeshOptimize2d::MeshOptimize2d()
    {
        SetFaceIndex(0);
        SetImproveEdges(0);
        SetMetricWeight(0);
        SetWriteStatus(1);
    }

    void MeshOptimize2d::ImproveMesh(Mesh& mesh, const MeshingParameters& mp)
    {
        if (!faceindex) {
            MESHIT_LOG_DEBUG("Smoothing");

            for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++) {
                ImproveMesh(mesh, mp);
            }
            faceindex = 0;
            return;
        }

        Opti2dLocalData ld;

        Array<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace(faceindex, seia);

        Vector x(2);

        Array<MeshPoint> savepoints(mesh.GetNP());

        ld.uselocalh = mp.uselocalh;

        Array<int> compress(mesh.GetNP());
        Array<PointIndex> icompress;
        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                compress[el[j]] = -1;
            }
        }
        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                if (compress[el[j]] == -1) {
                    compress[el[j]] = icompress.size();
                    icompress.push_back(el[j]);
                }
            }
        }
        Array<int> cnta(icompress.size());
        cnta = 0;
        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                cnta[compress[el[j]]]++;
            }
        }
        TABLE<SurfaceElementIndex> elementsonpoint(cnta);
        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                elementsonpoint.Add(compress[el[j]], seia[i]);
            }
        }

        ld.loch = mp.maxh;
        ld.locmetricweight = metricweight;

        Opti2SurfaceMinFunction surfminf(mesh, ld);
        OptiParameters par;
        par.maxit_linsearch = 8;
        par.maxit_bfgs = 5;

        for (size_t hi = 0; hi < icompress.size(); hi++) {
            PointIndex pi = icompress[hi];
            if (mesh[pi].Type() == SURFACEPOINT) {
                if (elementsonpoint[hi].size() == 0) continue;

                ld.sp1 = mesh[pi];
                ld.locelements.resize(0);
                ld.locrots.resize(0);
                ld.lochs.resize(0);
                ld.loc_pnts2.resize(0);
                ld.loc_pnts3.resize(0);

                for (size_t j = 0; j < elementsonpoint[hi].size(); j++) {
                    SurfaceElementIndex sei = elementsonpoint[hi][j];
                    const Element2d& bel = mesh.SurfaceElement(sei);
                    ld.locelements.push_back(sei);

                    for (size_t k = 1; k <= 3; k++) {
                        if (bel.PNum(k) == pi) {
                            ld.locrots.push_back(k);
                            ld.loc_pnts2.push_back(mesh[bel.PNumMod(k + 1)]);
                            ld.loc_pnts3.push_back(mesh[bel.PNumMod(k + 2)]);
                            break;
                        }
                    }
                    if (ld.uselocalh) {
                        Point3d pmid = Center(mesh[bel[0]], mesh[bel[1]], mesh[bel[2]]);
                        ld.lochs.push_back(mesh.GetH(pmid));
                    }
                }

                ld.normal = Vec3d(0, 0, 1);
                ld.t1 = Vec3d(0, 1, 0);
                ld.t2 = Vec3d(-1, 0, 0);

                // save points, and project to tangential plane
                for (size_t j = 0; j < ld.locelements.size(); j++) {
                    const Element2d& el = mesh.SurfaceElement(ld.locelements[j]);
                    for (size_t k = 0; k < 3; k++) {
                        savepoints[el[k]] = mesh[el[k]];
                    }
                }

                for (size_t j = 0; j < ld.locelements.size(); j++) {
                    const Element2d& el = mesh.SurfaceElement(ld.locelements[j]);
                    for (size_t k = 0; k < 3; k++) {
                        PointIndex hhpi = el[k];
                        double lam = ld.normal * (mesh[hhpi] - ld.sp1);
                        mesh.Point(hhpi) -= lam * ld.normal;
                    }
                }

                x = 0;
                par.typx = 0.3 * ld.lochs[0];
                BFGS(x, surfminf, par, 1e-6);

                // restore other points
                for (size_t j = 0; j < ld.locelements.size(); j++) {
                    const Element2d& el = mesh.SurfaceElement(ld.locelements[j]);
                    for (size_t k = 0; k < 3; k++) {
                        PointIndex hhpi = el[k];
                        if (hhpi != pi) mesh[hhpi] = savepoints[hhpi];
                    }
                }

                // optimizer pass (if whole distance is not possible, move only a bit!!!!)
                Vec3d hv(-x(1), x(0), 0.0);
                Point3d origp = mesh[pi];
                Point3d hnp = origp + hv;
                mesh.Point(pi).X() = hnp.X();
                mesh.Point(pi).Y() = hnp.Y();
                mesh.Point(pi).Z() = hnp.Z();
            }
        }
        mesh.SetNextTimeStamp();
    }

}  // namespace meshit
