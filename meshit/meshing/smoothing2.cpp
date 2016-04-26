/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2015-2016 Marc Honnorat <marc.honnorat@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library in the file LICENSE.LGPL; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */
#include "smoothing2.hpp"
#include "mesh_optimize.hpp"

namespace meshit {

static const double c_trig0 = 0.144337567297406;  // sqrt(3.0) / 12
static const double c_trig4 = 0.577350269189626;  // sqrt(3.0) / 3

double CalcTriangleBadness(const Point2d& p1, const Point2d& p2, const Point2d& p3,
                           double metric_weight, double h)
{
    // badness = B = sqrt(3.0)/12 * (\sum l_i^2) / area - 1
    double dx12 = p2.X() - p1.X();  // x component of e12 = p2 - p1
    double dx13 = p3.X() - p1.X();  // x component of e13 = p3 - p1
    double dx23 = p3.X() - p2.X();  // x component of e23 = p3 - p2
    double dy12 = p2.Y() - p1.Y();  // y component of e12 = p2 - p1
    double dy13 = p3.Y() - p1.Y();  // y component of e13 = p3 - p1
    double dy23 = p3.Y() - p2.Y();  // y component of e23 = p3 - p2

    // c^2 = Σ_i L_i^2 = e12.Length2() + e13.Length2() + e23.Length2()
    double c_2 = dx12 * dx12 + dx13 * dx13 + dx23 * dx23 + dy12 * dy12 + dy13 * dy13 + dy23 * dy23;

    // A = (1/2) * || e_12 ^ e_13 || = 0.5 * Cross(e12, e13).Length()
    double cross_z = dx12 * dy13 - dy12 * dx13;  // z component of Cross(e12, e13)
    double area = 0.5 * fabs(cross_z);

    if (area <= 1e-24 * c_2) return 1e10;

    // B = sqrt(3.0)/12 * ( c^2 / A ) - 1
    double badness = c_trig0 * c_2 / area - 1.0;

    if (metric_weight > 0) {
        // add:  metric_weight * (A / h^2 + h^2 / A - 2)
        // optimum for (2A) is h^2
        double areahh = 2.0 * area / (h * h);

        // B += metric_weight * (2A / h^2 + h^2 / 2A - 2)
        badness += metric_weight * (areahh + 1.0 / areahh - 2.0);
    }

    return badness;
}

double CalcTriangleBadnessGrad(const Point2d& p1, const Point2d& p2, const Point2d& p3,
                               Vec2d& d_bad, double metric_weight, double h)
{
    // badness = B = sqrt(3.0)/12 * (\sum l_i^2) / area - 1
    double dx12 = p2.X() - p1.X();  // x component of e12 = p2 - p1
    double dx13 = p3.X() - p1.X();  // x component of e13 = p3 - p1
    double dx23 = p3.X() - p2.X();  // x component of e23 = p3 - p2
    double dy12 = p2.Y() - p1.Y();  // y component of e12 = p2 - p1
    double dy13 = p3.Y() - p1.Y();  // y component of e13 = p3 - p1
    double dy23 = p3.Y() - p2.Y();  // y component of e23 = p3 - p2

    // c^2 = Σ_i L_i^2 = e12.Length2() + e13.Length2() + e23.Length2()
    double c_2 = dx12 * dx12 + dx13 * dx13 + dx23 * dx23 + dy12 * dy12 + dy13 * dy13 + dy23 * dy23;

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

    if (metric_weight > 0) {
        // add:  metric_weight * (A / h^2 + h^2 / A - 2)
        // optimum for (2A) is h^2
        double h_2 = h * h;
        double areahh = 2.0 * area / h_2;

        // B += metric_weight * (2A / h^2 + h^2 / 2A - 2)
        badness += metric_weight * (areahh + 1.0 / areahh - 2.0);

        // ∇B += metric_weight * ( 2∇A/h^2 - 2 h^2/((2A)^2)∇A )
        double gamma = 0.5 * metric_weight * (1.0 / h_2 - 0.25 * h_2 / (area * area)) / area;
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

    if (area <= 1e-24 * c_2) return 1e10;

    double badness = c_trig4 * c_2 / area - 1;
    return badness;
}

double CalcTriangleBadness_2(const Point2d& p1, const Point2d& p2, const Point2d& p3, double n_z)
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

    Point2d pp1(ld.sp1.X() - x[1], ld.sp1.Y() + x[0]);

    for (size_t j = 0; j < ld.loc_elements.size(); j++) {
        double loc_h = ld.lochs[j];
        const MeshPoint& pe2 = ld.loc_pnts2[j];
        const MeshPoint& pe3 = ld.loc_pnts3[j];
        double e2x = pe2.X() - pp1.X();
        double e2y = pe2.Y() - pp1.Y();
        double e3x = pe3.X() - pp1.X();
        double e3y = pe3.Y() - pp1.Y();

        if (e2x * e3y - e2y * e3x > 1e-8 * loc_h * loc_h) {
            badness += CalcTriangleBadness(pp1, Point2d(pe2), Point2d(pe3), ld.loc_metric_weight, loc_h);
        } else {
            badness += 1e8;
        }
    }
    return badness;
}

double MinFunction_2d::FuncGrad(const double* x, double* g)
{
    double badness = 0.0;
    Vec2d vgrad;
    Point2d pp1(ld.sp1.X() - x[1], ld.sp1.Y() + x[0]);

    for (size_t j = 0; j < ld.loc_elements.size(); j++) {
        double loc_h = ld.lochs[j];
        const MeshPoint& pe2 = ld.loc_pnts2[j];
        const MeshPoint& pe3 = ld.loc_pnts3[j];
        double e2x = pe2.X() - pp1.X();
        double e2y = pe2.Y() - pp1.Y();
        double e3x = pe3.X() - pp1.X();
        double e3y = pe3.Y() - pp1.Y();

        if (e2x * e3y - e2y * e3x > 1e-8 * loc_h * loc_h) {
            Vec2d hgrad;
            badness +=
                CalcTriangleBadnessGrad(pp1, Point2d(pe2), Point2d(pe3), hgrad, ld.loc_metric_weight, loc_h);
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

    Point2d pp1(ld.sp1.X() - x[1], ld.sp1.Y() + x[0]);

    for (size_t j = 0; j < ld.loc_elements.size(); j++) {
        double loc_h = ld.lochs[j];
        const MeshPoint& pe2 = ld.loc_pnts2[j];
        const MeshPoint& pe3 = ld.loc_pnts3[j];
        double e2x = pe2.X() - pp1.X();
        double e2y = pe2.Y() - pp1.Y();
        double e3x = pe3.X() - pp1.X();
        double e3y = pe3.Y() - pp1.Y();

        if (e2x * e3y - e2y * e3x > 1e-8 * loc_h * loc_h) {
            Vec2d hgrad;
            badness += CalcTriangleBadnessGrad(pp1, Point2d(pe2), Point2d(pe3), hgrad, ld.loc_metric_weight, loc_h);
            deriv += dir[0] * hgrad.Y() - dir[1] * hgrad.X();
        } else {
            badness += 1e8;
        }
    }
    return badness;
}

void MeshOptimize::ImproveMesh(double metric_weight)
{
    MESHIT_LOG_DEBUG("Smoothing");

    for (size_t face_index = 1; face_index <= mesh_.GetNbFaces(); face_index++) {
        ImproveMesh(face_index, metric_weight);
    }
}

void MeshOptimize::ImproveMesh(DomainIndex faceindex, double metric_weight)
{
    std::vector<ElementIndex> seia;
    mesh_.GetElementsOfFace(faceindex, seia);

    std::vector<MeshPoint> savepoints(mesh_.GetNbPoints());
    std::vector<PointIndex> compress(mesh_.GetNbPoints());
    std::vector<PointIndex> icompress;
    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& el = mesh_.Element(seia[i]);
        compress[el.PointID(0)] = CONST<PointIndex>::undefined;
        compress[el.PointID(1)] = CONST<PointIndex>::undefined;
        compress[el.PointID(2)] = CONST<PointIndex>::undefined;
    }
    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& el = mesh_.Element(seia[i]);
        for (size_t j = 0; j < 3; j++) {
            if (compress[el.PointID(j)] == CONST<PointIndex>::undefined) {
                compress[el.PointID(j)] = icompress.size();
                icompress.push_back(el.PointID(j));
            }
        }
    }
    std::vector<int> cnta(icompress.size(), 0);
    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& el = mesh_.Element(seia[i]);
        cnta[compress[el.PointID(0)]]++;
        cnta[compress[el.PointID(1)]]++;
        cnta[compress[el.PointID(2)]]++;
    }
    TABLE<ElementIndex> elements_on_point(cnta);
    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& el = mesh_.Element(seia[i]);
        elements_on_point.Add(compress[el.PointID(0)], seia[i]);
        elements_on_point.Add(compress[el.PointID(1)], seia[i]);
        elements_on_point.Add(compress[el.PointID(2)], seia[i]);
    }

    Opti2dLocalData ld(metric_weight);
    MinFunction_2d surfminf(mesh_, ld);
    OptiParameters par;
    par.maxit_linsearch = 8;
    par.maxit_bfgs = 5;
    par.eps = 1e-2;

    for (size_t hi = 0; hi < icompress.size(); hi++) {
        PointIndex pi = icompress[hi];
        MeshPoint& pp = mesh_.Point(pi);

        if (pp.Type() == INNER_POINT) {
            std::vector<ElementIndex> elem_idx = elements_on_point[hi];
            size_t n_elems = elem_idx.size();

            if (n_elems == 0) continue;

            ld.sp1 = pp;
            ld.loc_elements.resize(n_elems);
            ld.lochs.resize(n_elems);
            ld.loc_pnts2.resize(0);
            ld.loc_pnts3.resize(0);

            for (size_t j = 0; j < n_elems; j++) {
                ElementIndex sei = elem_idx[j];
                ld.loc_elements[j] = sei;

                const Element2d& bel = mesh_.Element(sei);
                for (size_t k = 0; k < 3; k++) {
                    if (bel[k] == pi) {
                        ld.loc_pnts2.push_back(mesh_.Point(bel[(k + 1) % 3]));
                        ld.loc_pnts3.push_back(mesh_.Point(bel[(k + 2) % 3]));
                        break;
                    }
                }
                Point2d pmid = Center(Point2d(mesh_.Point(bel[0])), Point2d(mesh_.Point(bel[1])),
                                      Point2d(mesh_.Point(bel[2])));
                ld.lochs[j] = mesh_.GetH(pmid);
            }

            // save points, and project to tangential plane
            for (size_t j = 0; j < n_elems; j++) {
                const Element2d& el = mesh_.Element(ld.loc_elements[j]);
                for (size_t k = 0; k < 3; k++) {
                    savepoints[el[k]] = mesh_.Point(el[k]);
                }
            }

            double x[2] = {0.0, 0.0};
            par.typx = 0.3 * ld.lochs[0];
            BFGS_2d(x, surfminf, par);

            // restore other points
            for (size_t j = 0; j < n_elems; j++) {
                const Element2d& el = mesh_.Element(ld.loc_elements[j]);
                for (size_t k = 0; k < 3; k++) {
                    PointIndex hhpi = el[k];
                    if (hhpi != pi) mesh_.Point(hhpi) = savepoints[hhpi];
                }
            }

            // optimizer pass (if whole distance is not possible, move only a bit!!!!)
            pp.X() -= x[1];
            pp.Y() += x[0];
        }
    }
}

}  // namespace meshit
