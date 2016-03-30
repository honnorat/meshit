#include "smoothing2.hpp"
#include "improve2.hpp"
#include "global.hpp"

namespace meshit
{
    static const double c_trig0 = 0.144337567297406;  // sqrt(3.0) / 12
    static const double c_trig4 = 0.577350269189626;  // sqrt(3.0) / 3

    double CalcTriangleBadness(const Point3d& p1,
                               const Point3d& p2,
                               const Point3d& p3,
                               double metricweight,
                               double h)
    {
        // badness = B = sqrt(3.0)/12 * (\sum l_i^2) / area - 1
        double dx12 = p2.X() - p1.X();  // x component of e12 = p2 - p1
        double dx13 = p3.X() - p1.X();  // x component of e13 = p3 - p1
        double dx23 = p3.X() - p2.X();  // x component of e23 = p3 - p2
        double dy12 = p2.Y() - p1.Y();  // y component of e12 = p2 - p1
        double dy13 = p3.Y() - p1.Y();  // y component of e13 = p3 - p1
        double dy23 = p3.Y() - p2.Y();  // y component of e23 = p3 - p2

        // c^2 = Σ_i L_i^2 = e12.Length2() + e13.Length2() + e23.Length2()
        double c_2 = dx12 * dx12 + dx13 * dx13 + dx23 * dx23 +
                     dy12 * dy12 + dy13 * dy13 + dy23 * dy23;

        // A = (1/2) * || e_12 ^ e_13 || = 0.5 * Cross(e12, e13).Length()
        double cross_z = dx12 * dy13 - dy12 * dx13;  // z component of Cross(e12, e13)
        double area = 0.5 * fabs(cross_z);

        if (area <= 1e-24 * c_2)
            return 1e10;

        // B = sqrt(3.0)/12 * ( c^2 / A ) - 1
        double badness = c_trig0 * c_2 / area - 1.0;

        if (metricweight > 0) {
            // add:  metricweight * (A / h^2 + h^2 / A - 2)
            // optimum for (2A) is h^2
            double areahh = 2.0 * area / (h * h);

            // B += metricweight * (2A / h^2 + h^2 / 2A - 2)
            badness += metricweight * (areahh + 1.0 / areahh - 2.0);
        }

        return badness;
    }

    double CalcTriangleBadnessGrad(const Point3d& p1,
                                   const Point3d& p2,
                                   const Point3d& p3,
                                   Vec3d& d_bad,
                                   double metricweight,
                                   double h)
    {
        // badness = B = sqrt(3.0)/12 * (\sum l_i^2) / area - 1
        double dx12 = p2.X() - p1.X();  // x component of e12 = p2 - p1
        double dx13 = p3.X() - p1.X();  // x component of e13 = p3 - p1
        double dx23 = p3.X() - p2.X();  // x component of e23 = p3 - p2
        double dy12 = p2.Y() - p1.Y();  // y component of e12 = p2 - p1
        double dy13 = p3.Y() - p1.Y();  // y component of e13 = p3 - p1
        double dy23 = p3.Y() - p2.Y();  // y component of e23 = p3 - p2

        // c^2 = Σ_i L_i^2 = e12.Length2() + e13.Length2() + e23.Length2()
        double c_2 = dx12 * dx12 + dx13 * dx13 + dx23 * dx23 +
                     dy12 * dy12 + dy13 * dy13 + dy23 * dy23;

        // A = (1/2) * || e_12 x e_13 || = 0.5 * Cross(e12, e13).Length()
        double cross_z = dx12 * dy13 - dy12 * dx13;  // z component of Cross(e12, e13)
        double area = 0.5 * fabs(cross_z);

        if (area <= 1e-24 * c_2) {
            d_bad = 0;
            return 1e10;
        }

        // B = sqrt(3.0)/12 * ( c^2 / A ) - 1
        double badness = c_trig0 * c_2 / area - 1.0;

        // ∇B = sqrt(3.0)/12 * (1/A * ∇(c^2) - c^2 / (A^2) * ∇A);
        //  with ∇(c^2) = -2*(e_12 + e_13)
        //  and  ∇A = 1/(4*A) * (e_32 x (e_12 x e_13 ) ) =  (0.25 / area) * Cross(-e23, Cross(e12, e13))
        double beta = 0.125 * c_2 * cross_z / (area * area);
        d_bad.X() = -2.0 * (c_trig0 / area) * (dx12 + dx13 - beta * dy23);
        d_bad.Y() = -2.0 * (c_trig0 / area) * (dy12 + dy13 + beta * dx23);
        d_bad.Z() = 0.0;

        if (metricweight > 0) {
            // add:  metricweight * (A / h^2 + h^2 / A - 2)
            // optimum for (2A) is h^2
            double h_2 = h * h;
            double areahh = 2.0 * area / h_2;

            // B += metricweight * (2A / h^2 + h^2 / 2A - 2)
            badness += metricweight * (areahh + 1.0 / areahh - 2.0);

            // ∇B += metricweight * ( 2∇A/h^2 - 2 h^2/((2A)^2)∇A )
            double gamma = 0.5 * metricweight * (1.0 / h_2 - 0.25 * h_2 / (area * area)) / area;
            d_bad.X() -= gamma * (cross_z * dy23);
            d_bad.Y() += gamma * (cross_z * dx23);
        }
        return badness;
    }


    inline double CalcTriangleBadnessRect(double x2, double x3, double y3)
    {
        // badness = sqrt(3.0) / 12 * (\sum l_i^2) / area - 1
        // p1 = (0, 0), p2 = (x2, 0), p3 = (x3, y3);

        double c_2 = (x2 * x2 + x3 * x3 + y3 * y3 - x2 * x3);
        double area = x2 * y3;

        if (area <= 1e-24 * c_2)
            return 1e10;

        double badness = c_trig4 * c_2 / area - 1;
        return badness;
    }

    double CalcTriangleBadness_2(const Point3d& p1,
                                 const Point3d& p2,
                                 const Point3d& p3,
                                 double n_z)
    {
        double dp12x = p2.X() - p1.X();  // v1 = p2 - p1
        double dp12y = p2.Y() - p1.Y();
        double dp13x = p3.X() - p1.X();  // v2 = p3 - p1
        double dp13y = p3.Y() - p1.Y();

        // x2 = e1 . v1  with  e1 = v1 / ||v1||  ie.  x2 = ||v1||
        // x3 = e1 . v2  with  e2 = (0, 0, n_z) x e2
        // y3 = e2 . v2
        double x2 = sqrt(dp12x * dp12x + dp12y * dp12y);
        double x3 = (dp12x * dp13x + dp12y * dp13y);
        double y3 = n_z * (dp13y * dp12x - dp12y * dp13x) / x2;

        return CalcTriangleBadnessRect(x2, x3, y3);
    }

    double MinFunction_2d::Func(const double* x)
    {
        double badness = 0;

        Point3d pp1(ld.sp1.X() - x[1],
                    ld.sp1.Y() + x[0]);

        for (size_t j = 0; j < ld.loc_elements.size(); j++) {
            double loc_h = ld.lochs[j];
            const Point3d& pe2 = ld.loc_pnts2[j];
            const Point3d& pe3 = ld.loc_pnts3[j];
            Point3d e2(pe2.X() - pp1.X(), pe2.Y() - pp1.Y());
            Point3d e3(pe3.X() - pp1.X(), pe3.Y() - pp1.Y());

            if (e2.X() * e3.Y() - e2.Y() * e3.X() > 1e-8 * loc_h * loc_h) {
                badness += CalcTriangleBadness(pp1, pe2, pe3, ld.loc_metric_weight, loc_h);
            } else {
                badness += 1e8;
            }
        }
        return badness;
    }

    double MinFunction_2d::FuncGrad(const double* x, double* g)
    {
        double badness = 0.0;
        Vec3d vgrad;
        Point3d pp1(ld.sp1.X() - x[1],
                    ld.sp1.Y() + x[0]);

        for (size_t j = 0; j < ld.loc_elements.size(); j++) {
            double loc_h = ld.lochs[j];
            const Point3d& pe2 = ld.loc_pnts2[j];
            const Point3d& pe3 = ld.loc_pnts3[j];
            Point3d e2(pe2.X() - pp1.X(), pe2.Y() - pp1.Y());
            Point3d e3(pe3.X() - pp1.X(), pe3.Y() - pp1.Y());

            if (e2.X() * e3.Y() - e2.Y() * e3.X() > 1e-8 * loc_h * loc_h) {
                Vec3d hgrad;
                badness += CalcTriangleBadnessGrad(pp1, pe2, pe3, hgrad, ld.loc_metric_weight, loc_h);
                vgrad += hgrad;
            } else {
                badness += 1e8;
            }
        }
        g[0] = +vgrad.Y();
        g[1] = -vgrad.X();
        return badness;
    }

    double MinFunction_2d::FuncDeriv(const double* x, const double* dir, double& deriv)
    {
        deriv = 0;
        double badness = 0;

        Point3d pp1(ld.sp1.X() - x[1],
                    ld.sp1.Y() + x[0]);

        for (size_t j = 0; j < ld.loc_elements.size(); j++) {
            double loc_h = ld.lochs[j];
            const Point3d& pe2 = ld.loc_pnts2[j];
            const Point3d& pe3 = ld.loc_pnts3[j];
            Point3d e2(pe2.X() - pp1.X(), pe2.Y() - pp1.Y());
            Point3d e3(pe3.X() - pp1.X(), pe3.Y() - pp1.Y());

            if (e2.X() * e3.Y() - e2.Y() * e3.X() > 1e-8 * loc_h * loc_h) {
                Vec3d hgrad;
                badness += CalcTriangleBadnessGrad(pp1, pe2, pe3, hgrad, ld.loc_metric_weight, loc_h);
                deriv += dir[0] * hgrad.Y() - dir[1] * hgrad.X();
            } else {
                badness += 1e8;
            }
        }
        return badness;
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

        std::vector<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace(faceindex, seia);


        std::vector<MeshPoint> savepoints(mesh.GetNP());
        std::vector<int> compress(mesh.GetNP());
        std::vector<PointIndex> icompress;
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
        std::vector<int> cnta(icompress.size(), 0);
        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                cnta[compress[el[j]]]++;
            }
        }
        TABLE<SurfaceElementIndex> elements_on_point(cnta);
        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                elements_on_point.Add(compress[el[j]], seia[i]);
            }
        }

        ld.loc_metric_weight = metricweight;

        MinFunction_2d surfminf(mesh, ld);
        OptiParameters par;
        par.maxit_linsearch = 8;
        par.maxit_bfgs = 5;

        for (size_t hi = 0; hi < icompress.size(); hi++) {
            PointIndex pi = icompress[hi];
            MeshPoint& pp = mesh.Point(pi);

            if (pp.Type() == SURFACEPOINT) {
                std::vector<SurfaceElementIndex> elem_idx = elements_on_point[hi];
                size_t n_elems = elem_idx.size();

                if (n_elems == 0) continue;

                ld.sp1 = pp;
                ld.loc_elements.resize(n_elems);
                ld.lochs.resize(n_elems);
                ld.loc_pnts2.resize(0);
                ld.loc_pnts3.resize(0);

                for (size_t j = 0; j < n_elems; j++) {
                    SurfaceElementIndex sei = elem_idx[j];
                    ld.loc_elements[j] = sei;

                    const Element2d& bel = mesh.SurfaceElement(sei);
                    for (size_t k = 1; k <= 3; k++) {
                        if (bel.PNum(k) == pi) {
                            ld.loc_pnts2.push_back(mesh[bel.PNumMod(k + 1)]);
                            ld.loc_pnts3.push_back(mesh[bel.PNumMod(k + 2)]);
                            break;
                        }
                    }
                    Point3d pmid = Center(mesh[bel[0]], mesh[bel[1]], mesh[bel[2]]);
                    ld.lochs[j] = mesh.GetH(pmid);
                }

                // save points, and project to tangential plane
                for (size_t j = 0; j < n_elems; j++) {
                    const Element2d& el = mesh.SurfaceElement(ld.loc_elements[j]);
                    for (size_t k = 0; k < 3; k++) {
                        savepoints[el[k]] = mesh[el[k]];
                    }
                }

                double x[2] = {0.0, 0.0};
                par.typx = 0.3 * ld.lochs[0];
                BFGS_2d(x, surfminf, par, 1e-6);

                // restore other points
                for (size_t j = 0; j < n_elems; j++) {
                    const Element2d& el = mesh.SurfaceElement(ld.loc_elements[j]);
                    for (size_t k = 0; k < 3; k++) {
                        PointIndex hhpi = el[k];
                        if (hhpi != pi) mesh.Point(hhpi) = savepoints[hhpi];
                    }
                }

                // optimizer pass (if whole distance is not possible, move only a bit!!!!)
                pp.X() -= x[1];
                pp.Y() += x[0];
            }
        }
    }

}  // namespace meshit
