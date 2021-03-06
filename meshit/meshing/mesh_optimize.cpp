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

#include "mesh_optimize.hpp"

namespace meshit {

class Neighbour
{
 public:
    Neighbour() { }

    void SetNeighborIndex(size_t side, ElementIndex anr) { nr[side] = anr; }
    ElementIndex GetNeighborIndex(size_t side) { return nr[side]; }

    void SetOrientation(size_t side, size_t aorient) { orient[side] = aorient; }
    size_t GetOrientation(size_t side) { return orient[side]; }

 protected:
    ElementIndex nr[3];
    size_t orient[3];
};

class trionedge
{
 public:
    GenericIndex tnr;
    size_t sidenr;

    trionedge()
        : tnr{0}, sidenr{0} { }

    trionedge(GenericIndex atnr, size_t asidenr)
        : tnr{atnr}, sidenr{asidenr} { }
};

void MeshOptimize::EdgeSwapping(bool use_metric, double metric_weight)
{
    if (use_metric) {
        MESHIT_LOG_DEBUG("Edgeswapping, metric");
    } else {
        MESHIT_LOG_DEBUG("Edgeswapping, topological");
    }

    for (DomainIndex face_index = 1; face_index <= mesh_.GetNbFaces(); face_index++) {
        EdgeSwapping(face_index, use_metric, metric_weight);
    }

    mesh_.IndexBoundaryEdges();
}

void MeshOptimize::EdgeSwapping(DomainIndex face_index, bool use_metric, double metric_weight)
{
    std::vector<ElementIndex> seia;
    mesh_.GetElementsOfFace(face_index, seia);

    std::vector<Neighbour> neighbors(mesh_.GetNbElements());
    IndexPair_map<trionedge> other(seia.size() + 2);

    std::vector<bool> swapped(mesh_.GetNbElements());
    std::vector<int> pdef(mesh_.GetNbPoints());
    std::vector<double> pangle(mesh_.GetNbPoints());

    static const double minangle[] = {0, 1.481, 2.565, 3.627, 4.683, 5.736, 7, 9};

    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& sel = mesh_.Element(seia[i]);
        for (int j = 0; j < 3; j++) {
            pangle[sel[j]] = 0.0;
        }
    }

    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& sel = mesh_.Element(seia[i]);
        for (size_t j = 0; j < 3; j++) {
            point_type_t typ = mesh_.Point(sel[j]).Type();
            if (typ == EDGE_POINT) {
                pangle[sel[j]] += Angle(mesh_.Point(sel[(j + 1) % 3]) - mesh_.Point(sel[j]),
                                        mesh_.Point(sel[(j + 2) % 3]) - mesh_.Point(sel[j]));
            }
        }
    }

    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& sel = mesh_.Element(seia[i]);
        for (size_t j = 0; j < 3; j++) {
            PointIndex pi = sel[j];
            if (mesh_.Point(pi).Type() == INNER_POINT) {
                pdef[pi] = -6;
            } else {
                for (int k = 0; k < 8; k++) {
                    if (pangle[pi] >= minangle[k]) {
                        pdef[pi] = -1 - k;
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& sel = mesh_.Element(seia[i]);
        for (size_t j = 0; j < 3; j++) {
            pdef[sel[j]]++;
        }
    }

    for (size_t i = 0; i < seia.size(); i++) {
        ElementIndex el_idx = seia[i];
        Neighbour& neighb = neighbors[el_idx];
        for (size_t j = 0; j < 3; j++) {
            neighb.SetNeighborIndex(j, CONST<ElementIndex>::undefined);
            neighb.SetOrientation(j, 0);
        }
    }

    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& sel = mesh_.Element(seia[i]);

        for (size_t j = 0; j < 3; j++) {
            PointIndex pi1 = sel.PointID((j + 1) % 3);
            PointIndex pi2 = sel.PointID((j + 2) % 3);

            IndexPair edge(pi1, pi2);
            edge.Sort();

            if (mesh_.IsSegment(edge)) continue;

            IndexPair ii2(pi1, pi2);
            if (other.count(ii2)) {
                GenericIndex i2 = other[ii2].tnr;
                size_t j2 = other[ii2].sidenr;
                neighbors[seia[i]].SetNeighborIndex(j, i2);
                neighbors[seia[i]].SetOrientation(j, j2);
                neighbors[i2].SetNeighborIndex(j2, seia[i]);
                neighbors[i2].SetOrientation(j2, j);
            } else {
                other[IndexPair(pi2, pi1)] = trionedge(seia[i], j);
            }
        }
    }
    for (size_t i = 0; i < seia.size(); i++) {
        swapped[seia[i]] = false;
    }
    int t = 4;
    bool done = false;
    while (!done && t >= 2) {
        for (size_t i = 0; i < seia.size(); i++) {
            ElementIndex t1 = seia[i];

            if (mesh_.Element(t1).IsDeleted()) continue;

            if (mesh_.Element(t1).FaceID() != face_index) continue;

            for (size_t o1 = 0; o1 < 3; o1++) {
                bool should;

                ElementIndex t2 = neighbors[t1].GetNeighborIndex(o1);
                size_t o2 = neighbors[t1].GetOrientation(o1);

                if (t2 == CONST<ElementIndex>::undefined) continue;
                if (swapped[t1] || swapped[t2]) continue;

                PointIndex pi1 = mesh_.Element(t1).PointID((o1 + 1) % 3);
                PointIndex pi2 = mesh_.Element(t1).PointID((o1 + 2) % 3);
                PointIndex pi3 = mesh_.Element(t1).PointID(o1 % 3);
                PointIndex pi4 = mesh_.Element(t2).PointID(o2 % 3);

                const MeshPoint& pm_1 = mesh_.Point(pi1);
                const MeshPoint& pm_2 = mesh_.Point(pi2);
                const MeshPoint& pm_3 = mesh_.Point(pi3);
                const MeshPoint& pm_4 = mesh_.Point(pi4);

                Point2d p_1 = Point2d(pm_1);
                Point2d p_2 = Point2d(pm_2);
                Point2d p_3 = Point2d(pm_3);
                Point2d p_4 = Point2d(pm_4);

                bool allowswap;

                Vec2d auxvec1 = p_3 - p_4;
                Vec2d auxvec2 = p_1 - p_4;

                allowswap = fabs(1.0 - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

                if (!allowswap) continue;

                // normal of new
                double bv1 = auxvec1.X() * auxvec2.Y() - auxvec1.Y() * auxvec2.X() > 0;

                Vec2d auxvec3 = p_4 - p_3;
                Vec2d auxvec4 = p_2 - p_3;

                allowswap = fabs(1. - (auxvec3 * auxvec4) / (auxvec3.Length() * auxvec4.Length())) > 1e-4;

                if (!allowswap) continue;

                bool bv2 = auxvec3.X() * auxvec4.Y() - auxvec3.Y() * auxvec4.X() > 0.0;

                // normals of original
                Vec2d auxvec5 = p_2 - p_4;
                Vec2d auxvec6 = p_1 - p_3;

                bool bv3 = auxvec2.X() * auxvec5.Y() - auxvec2.Y() * auxvec5.X() < 0.0;
                bool bv4 = auxvec4.X() * auxvec6.Y() - auxvec4.Y() * auxvec6.X() < 0.0;

                allowswap = bv1 && bv2 && bv3 && bv4;

                if (!allowswap) continue;

                if (use_metric) {
                    double loc_h = mesh_.GetH(p_1);
                    double bad1 = CalcTriangleBadness(p_4, p_3, p_1, metric_weight, loc_h);
                    double bad2 = CalcTriangleBadness(p_3, p_4, p_2, metric_weight, loc_h);
                    double bad3 = CalcTriangleBadness(p_1, p_2, p_3, metric_weight, loc_h);
                    double bad4 = CalcTriangleBadness(p_2, p_1, p_4, metric_weight, loc_h);
                    should = (bad1 + bad2) < (bad3 + bad4);
                } else {
                    int e = pdef[pi1] + pdef[pi2] - pdef[pi3] - pdef[pi4];
                    double d = Dist2(p_1, p_2) - Dist2(p_3, p_4);
                    should = e >= t && (e > 2 || d > 0);
                }

                if (should) {
                    // do swapping !
                    done = true;

                    mesh_.Element(t1).PointID(0) = pi1;
                    mesh_.Element(t1).PointID(1) = pi4;
                    mesh_.Element(t1).PointID(2) = pi3;

                    mesh_.Element(t2).PointID(0) = pi2;
                    mesh_.Element(t2).PointID(1) = pi3;
                    mesh_.Element(t2).PointID(2) = pi4;

                    pdef[pi1]--;
                    pdef[pi2]--;
                    pdef[pi3]++;
                    pdef[pi4]++;

                    swapped[t1] = true;
                    swapped[t2] = true;
                }
            }
        }
        t--;
    }
}

void MeshOptimize::CombineImprove()
{
    MESHIT_LOG_DEBUG("Combine improve");
    for (DomainIndex face_index = 1; face_index <= mesh_.GetNbFaces(); face_index++) {
        CombineImprove(face_index);
    }
}

void MeshOptimize::CombineImprove(DomainIndex face_index)
{
    std::vector<ElementIndex> seia;
    mesh_.GetElementsOfFace(face_index, seia);

    size_t np = mesh_.GetNbPoints();

    TABLE<ElementIndex> elements_on_node(np);
    std::vector<ElementIndex> has_one_pi, has_both_pi;

    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& el = mesh_.Element(seia[i]);
        elements_on_node.Add(el.PointID(0), seia[i]);
        elements_on_node.Add(el.PointID(1), seia[i]);
        elements_on_node.Add(el.PointID(2), seia[i]);
    }

    std::vector<bool> fixed(np, false);

    for (size_t i = 0; i < seia.size(); i++) {
        const Element2d& sel = mesh_.Element(seia[i]);
        for (size_t j = 0; j < 3; j++) {
            PointIndex pi1 = sel.PointID((j + 1) % 3);
            PointIndex pi2 = sel.PointID((j + 2) % 3);
            if (mesh_.IsSegment(IndexPair(pi1, pi2).Sort())) {
                fixed[pi1] = true;
                fixed[pi2] = true;
            }
        }
    }

    for (size_t i = 0; i < seia.size(); i++) {
        ElementIndex sei = seia[i];
        Element2d& elem = mesh_.Element(sei);
        if (elem.IsDeleted()) continue;

        for (size_t j = 0; j < 3; j++) {
            PointIndex pi1 = elem[j];
            PointIndex pi2 = elem[(j + 1) % 3];

            if (pi1 < 0 || pi2 < 0) continue;

            // more general
            if (fixed[pi2]) {
                std::swap(pi1, pi2);
            }
            if (fixed[pi2]) {
                continue;
            }

            IndexPair si2(pi1, pi2);
            si2.Sort();

            double nv_z = 0.0;
            has_one_pi.resize(0);
            has_both_pi.resize(0);

            std::vector<ElementIndex> elem_idx_1 = elements_on_node[pi1];
            std::vector<ElementIndex> elem_idx_2 = elements_on_node[pi2];

            for (size_t k = 0; k < elem_idx_1.size(); k++) {
                const Element2d& el2 = mesh_.Element(elem_idx_1[k]);

                if (el2.IsDeleted()) continue;

                if (el2[0] == pi2 || el2[1] == pi2 || el2[2] == pi2) {
                    has_both_pi.push_back(elem_idx_1[k]);

                    const MeshPoint& p1 = mesh_.Point(el2[0]);
                    const MeshPoint& p2 = mesh_.Point(el2[1]);
                    const MeshPoint& p3 = mesh_.Point(el2[2]);

                    // Vec3d nv = Cross(p2 -p1, p3 - p1);
                    nv_z = (p2.X() - p1.X()) * (p3.Y() - p1.Y()) - (p2.Y() - p1.Y()) * (p3.X() - p1.X());
                } else {
                    has_one_pi.push_back(elem_idx_1[k]);
                }
            }

            for (size_t k = 0; k < elem_idx_2.size(); k++) {
                const Element2d& el2 = mesh_.Element(elem_idx_2[k]);
                if (el2.IsDeleted()) continue;

                if (el2[0] != pi1 && el2[1] != pi1 && el2[2] != pi1) {
                    has_one_pi.push_back(elem_idx_2[k]);
                }
            }

            double bad1 = 0.0;
            for (size_t k = 0; k < has_one_pi.size(); k++) {
                const Element2d& el = mesh_.Element(has_one_pi[k]);
                bad1 += CalcTriangleBadness_2(Point2d(mesh_.Point(el[0])), Point2d(mesh_.Point(el[1])),
                                              Point2d(mesh_.Point(el[2])), nv_z);
            }

            for (size_t k = 0; k < has_both_pi.size(); k++) {
                const Element2d& el = mesh_.Element(has_both_pi[k]);
                bad1 += CalcTriangleBadness_2(Point2d(mesh_.Point(el[0])), Point2d(mesh_.Point(el[1])),
                                              Point2d(mesh_.Point(el[2])), nv_z);
            }
            bad1 /= (has_one_pi.size() + has_both_pi.size());

            MeshPoint p1_save = mesh_.Point(pi1);
            MeshPoint p2_save = mesh_.Point(pi2);

            MeshPoint pnew = p1_save;
            mesh_.Point(pi1) = pnew;
            mesh_.Point(pi2) = pnew;

            double bad2 = 0;
            for (size_t k = 0; k < has_one_pi.size(); k++) {
                const Element2d& el = mesh_.Element(has_one_pi[k]);
                const MeshPoint& p1 = mesh_.Point(el[0]);
                const MeshPoint& p2 = mesh_.Point(el[1]);
                const MeshPoint& p3 = mesh_.Point(el[2]);
                bad2 += CalcTriangleBadness_2(Point2d(p1), Point2d(p2), Point2d(p3), nv_z);

                // Vec3d hnv = Cross(p2 - p1, p3 - p1);
                double hnv_z = (p2.X() - p1.X()) * (p3.Y() - p1.Y()) - (p2.Y() - p1.Y()) * (p3.X() - p1.X());
                if (hnv_z * nv_z < 0) {
                    bad2 += 1e10;
                }
                if (nv_z < 0.5) {
                    bad2 += 1e10;
                }
            }
            bad2 /= has_one_pi.size();

            mesh_.Point(pi1) = p1_save;
            mesh_.Point(pi2) = p2_save;

            bool should = (bad2 < bad1 && bad2 < 1e4);

            if (should) {
                mesh_.Point(pi1) = pnew;
                for (size_t k = 0; k < elem_idx_2.size(); k++) {
                    Element2d& el = mesh_.Element(elem_idx_2[k]);
                    if (el.IsDeleted()) continue;
                    elements_on_node.Add(pi1, elem_idx_2[k]);

                    bool haspi1 = 0;
                    for (size_t l = 0; l < 3; l++) {
                        if (el[l] == pi1) haspi1 = 1;
                    }
                    if (haspi1) continue;

                    for (size_t l = 0; l < 3; l++) {
                        if (el[l] == pi2) {
                            el[l] = pi1;
                        }
                        fixed[el[l]] = true;
                    }
                }
                for (size_t k = 0; k < has_both_pi.size(); k++) {
                    mesh_.Element(has_both_pi[k]).Delete();
                }
            }
        }
    }
    mesh_.Compress();
}

}  // namespace meshit
