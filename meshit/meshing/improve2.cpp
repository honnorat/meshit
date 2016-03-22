#include "improve2.hpp"
#include "global.hpp"

namespace meshit {

    class Neighbour
    {
        int nr[3];
        int orient[3];

     public:
        Neighbour() { }

        void SetNr(int side, int anr)
        {
            nr[side] = anr;
        }

        int GetNr(int side)
        {
            return nr[side];
        }

        void SetOrientation(size_t side, int aorient)
        {
            orient[side] = aorient;
        }

        int GetOrientation(size_t side)
        {
            return orient[side];
        }
    };

    class trionedge
    {
     public:
        size_t tnr;
        size_t sidenr;

        trionedge() : tnr{0}, sidenr{0} { }

        trionedge(size_t atnr, size_t asidenr) : tnr{atnr}, sidenr{asidenr} { }
    };

    void MeshOptimize2d::EdgeSwapping(Mesh& mesh, int usemetric)
    {
        if (!faceindex) {
            if (usemetric)
                MESHIT_LOG_DEBUG("Edgeswapping, metric");
            else
                MESHIT_LOG_DEBUG("Edgeswapping, topological");

            for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++) {
                EdgeSwapping(mesh, usemetric);
            }

            faceindex = 0;
            mesh.CalcSurfacesOfNode();
            return;
        }

        std::vector<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace(faceindex, seia);

        std::vector<Neighbour> neighbors(mesh.GetNSE());
        INDEX_2_HASHTABLE<trionedge> other(seia.size() + 2);

        std::vector<char> swapped(mesh.GetNSE());
        std::vector<int> pdef(mesh.GetNP());
        std::vector<double> pangle(mesh.GetNP());

        static const double minangle[] = {0, 1.481, 2.565, 3.627, 4.683, 5.736, 7, 9};

        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++) {
                pangle[sel[j]] = 0.0;
            }
        }
        // pangle = 0;

        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                POINTTYPE typ = mesh.Point(sel[j]).Type();
                if (typ == FIXEDPOINT || typ == EDGEPOINT) {
                    pangle[sel[j]] += Angle(
                            mesh.Point(sel[(j + 1) % 3]) - mesh.Point(sel[j]),
                            mesh.Point(sel[(j + 2) % 3]) - mesh.Point(sel[j]));
                }
            }
        }

        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                PointIndex pi = sel[j];
                if (mesh.Point(pi).Type() == INNERPOINT || mesh.Point(pi).Type() == SURFACEPOINT) {
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
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++)
                pdef[sel[j]]++;
        }

        for (size_t i = 0; i < seia.size(); i++) {
            for (size_t j = 0; j < 3; j++) {
                neighbors[seia[i]].SetNr(j, -1);
                neighbors[seia[i]].SetOrientation(j, 0);
            }
        }

        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);

            for (size_t j = 0; j < 3; j++) {
                PointIndex pi1 = sel.PNumMod(j + 2);
                PointIndex pi2 = sel.PNumMod(j + 3);

                INDEX_2 edge(pi1, pi2);
                edge.Sort();

                if (mesh.IsSegment(pi1, pi2))
                    continue;

                INDEX_2 ii2(pi1, pi2);
                if (other.Used(ii2)) {
                    size_t i2 = other.Get(ii2).tnr;
                    size_t j2 = other.Get(ii2).sidenr;
                    neighbors[seia[i]].SetNr(j, i2);
                    neighbors[seia[i]].SetOrientation(j, j2);
                    neighbors[i2].SetNr(j2, seia[i]);
                    neighbors[i2].SetOrientation(j2, j);
                } else {
                    other.Set(INDEX_2(pi2, pi1), trionedge(seia[i], j));
                }
            }
        }
        for (size_t i = 0; i < seia.size(); i++) {
            swapped[seia[i]] = 0;
        }
        int t = 4;
        bool done = false;
        while (!done && t >= 2) {
            for (size_t i = 0; i < seia.size(); i++) {
                SurfaceElementIndex t1 = seia[i];

                if (mesh.SurfaceElement(t1).IsDeleted())
                    continue;

                if (mesh.SurfaceElement(t1).GetIndex() != faceindex)
                    continue;

                for (size_t o1 = 0; o1 < 3; o1++) {
                    bool should;

                    SurfaceElementIndex t2 = neighbors[t1].GetNr(o1);
                    size_t o2 = neighbors[t1].GetOrientation(o1);

                    if (t2 == -1) continue;
                    if (swapped[t1] || swapped[t2]) continue;

                    PointIndex pi1 = mesh.SurfaceElement(t1).PNumMod(o1 + 1 + 1);
                    PointIndex pi2 = mesh.SurfaceElement(t1).PNumMod(o1 + 1 + 2);
                    PointIndex pi3 = mesh.SurfaceElement(t1).PNumMod(o1 + 1);
                    PointIndex pi4 = mesh.SurfaceElement(t2).PNumMod(o2 + 1);

                    bool allowswap = true;

                    Vec3d auxvec1 = mesh.Point(pi3) - mesh.Point(pi4);
                    Vec3d auxvec2 = mesh.Point(pi1) - mesh.Point(pi4);

                    allowswap &= fabs(1. - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

                    if (!allowswap)
                        continue;

                    // normal of new
                    Vec3d nv1 = Cross(auxvec1, auxvec2);

                    auxvec1 = mesh.Point(pi4) - mesh.Point(pi3);
                    auxvec2 = mesh.Point(pi2) - mesh.Point(pi3);
                    allowswap &= fabs(1. - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

                    if (!allowswap)
                        continue;

                    Vec3d nv2 = Cross(auxvec1, auxvec2);

                    // normals of original
                    Vec3d nv3 = Cross(mesh.Point(pi1) - mesh.Point(pi4), mesh.Point(pi2) - mesh.Point(pi4));
                    Vec3d nv4 = Cross(mesh.Point(pi2) - mesh.Point(pi3), mesh.Point(pi1) - mesh.Point(pi3));

                    nv3 *= -1;
                    nv4 *= -1;
                    nv3.Normalize();
                    nv4.Normalize();

                    nv1.Normalize();
                    nv2.Normalize();

                    double critval = cos(M_PI / 6);  // 30 degree
                    allowswap &= (nv1.Z() > critval) &&
                                 (nv2.Z() > critval) &&
                                 (nv3.Z() > critval) &&
                                 (nv4.Z() > critval);

                    double horder = Dist(mesh.Point(pi1), mesh.Point(pi2));

                    if (allowswap &&
                        nv1.Length() > 1e-3 * horder * horder &&
                        nv2.Length() > 1e-3 * horder * horder) {
                        if (!usemetric) {
                            int e = pdef[pi1] + pdef[pi2] - pdef[pi3] - pdef[pi4];
                            double d = Dist2(mesh.Point(pi1), mesh.Point(pi2)) -
                                       Dist2(mesh.Point(pi3), mesh.Point(pi4));

                            should = e >= t && (e > 2 || d > 0);
                        } else {
                            double loch = mesh.GetH(mesh.Point(pi1));
                            double bad1 = CalcTriangleBadness(mesh.Point(pi4), mesh.Point(pi3), mesh.Point(pi1),
                                                              metricweight, loch);
                            double bad2 = CalcTriangleBadness(mesh.Point(pi3), mesh.Point(pi4), mesh.Point(pi2),
                                                              metricweight, loch);
                            double bad3 = CalcTriangleBadness(mesh.Point(pi1), mesh.Point(pi2), mesh.Point(pi3),
                                                              metricweight, loch);
                            double bad4 = CalcTriangleBadness(mesh.Point(pi2), mesh.Point(pi1), mesh.Point(pi4),
                                                              metricweight, loch);
                            should = (bad1 + bad2) < (bad3 + bad4);
                        }

                        if (should) {
                            // do swapping !
                            done = true;

                            mesh.SurfaceElement(t1).PNum(1) = pi1;
                            mesh.SurfaceElement(t1).PNum(2) = pi4;
                            mesh.SurfaceElement(t1).PNum(3) = pi3;

                            mesh.SurfaceElement(t2).PNum(1) = pi2;
                            mesh.SurfaceElement(t2).PNum(2) = pi3;
                            mesh.SurfaceElement(t2).PNum(3) = pi4;

                            pdef[pi1]--;
                            pdef[pi2]--;
                            pdef[pi3]++;
                            pdef[pi4]++;

                            swapped[t1] = 1;
                            swapped[t2] = 1;
                        }
                    }
                }
            }
            t--;
        }

        mesh.SetNextTimeStamp();
    }

    void MeshOptimize2d::CombineImprove(Mesh& mesh)
    {
        if (!faceindex) {
            MESHIT_LOG_DEBUG("Combine improve");

            for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++) {
                CombineImprove(mesh);
            }
            faceindex = 0;
            return;
        }

        std::vector<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace(faceindex, seia);

        size_t np = mesh.GetNP();

        TABLE<SurfaceElementIndex> elements_on_node(np);
        std::vector<SurfaceElementIndex> has_one_pi, has_both_pi;

        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                elements_on_node.Add(el[j], seia[i]);
            }
        }

        std::vector<bool> fixed(np, false);

        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < 3; j++) {
                PointIndex pi1 = sel.PNumMod(j + 2);
                PointIndex pi2 = sel.PNumMod(j + 3);
                if (mesh.IsSegment(pi1, pi2)) {
                    fixed[pi1] = true;
                    fixed[pi2] = true;
                }
            }
        }
        for (size_t i = 0; i < mesh.LockedPoints().size(); i++) {
            fixed[mesh.LockedPoints()[i]] = true;
        }

        for (size_t i = 0; i < seia.size(); i++) {
            SurfaceElementIndex sei = seia[i];
            Element2d& elem = mesh.SurfaceElement(sei);
            if (elem.IsDeleted()) continue;

            for (size_t j = 0; j < 3; j++) {
                PointIndex pi1 = elem[j];
                PointIndex pi2 = elem[(j + 1) % 3];

                if (pi1 < 0 || pi2 < 0)
                    continue;

                // more general
                if (fixed[pi2]) {
                    std::swap(pi1, pi2);
                }
                if (fixed[pi2]) {
                    continue;
                }

                INDEX_2 si2(pi1, pi2);
                si2.Sort();

                double nv_z = 0.0;
                has_one_pi.resize(0);
                has_both_pi.resize(0);

                std::vector<SurfaceElementIndex> elem_idx_1 = elements_on_node[pi1];
                std::vector<SurfaceElementIndex> elem_idx_2 = elements_on_node[pi2];

                for (size_t k = 0; k < elem_idx_1.size(); k++) {
                    const Element2d& el2 = mesh.SurfaceElement(elem_idx_1[k]);

                    if (el2.IsDeleted()) continue;

                    if (el2[0] == pi2 || el2[1] == pi2 || el2[2] == pi2) {
                        has_both_pi.push_back(elem_idx_1[k]);

                        const Point3d& p1 = mesh.Point(el2[0]);
                        const Point3d& p2 = mesh.Point(el2[1]);
                        const Point3d& p3 = mesh.Point(el2[2]);

                        // Vec3d nv = Cross(p2 -p1, p3 - p1);
                        nv_z = (p2.X() - p1.X()) * (p3.Y() - p1.Y()) - (p2.Y() - p1.Y()) * (p3.X() - p1.X());
                    } else {
                        has_one_pi.push_back(elem_idx_1[k]);
                    }
                }

                for (size_t k = 0; k < elem_idx_2.size(); k++) {
                    const Element2d& el2 = mesh.SurfaceElement(elem_idx_2[k]);
                    if (el2.IsDeleted()) continue;

                    if (el2[0] != pi1 && el2[1] != pi1 && el2[2] != pi1) {
                        has_one_pi.push_back(elem_idx_2[k]);
                    }
                }

                double bad1 = 0.0;
                for (size_t k = 0; k < has_one_pi.size(); k++) {
                    const Element2d& el = mesh.SurfaceElement(has_one_pi[k]);
                    bad1 += CalcTriangleBadness_2(mesh.Point(el[0]), mesh.Point(el[1]), mesh.Point(el[2]), nv_z);
                }

                for (size_t k = 0; k < has_both_pi.size(); k++) {
                    const Element2d& el = mesh.SurfaceElement(has_both_pi[k]);
                    bad1 += CalcTriangleBadness_2(mesh.Point(el[0]), mesh.Point(el[1]), mesh.Point(el[2]), nv_z);
                }
                bad1 /= (has_one_pi.size() + has_both_pi.size());

                MeshPoint p1 = mesh[pi1];
                MeshPoint p2 = mesh[pi2];

                MeshPoint pnew = p1;
                mesh[pi1] = pnew;
                mesh[pi2] = pnew;

                double bad2 = 0;
                for (size_t k = 0; k < has_one_pi.size(); k++) {
                    const Element2d& el = mesh.SurfaceElement(has_one_pi[k]);
                    const Point3d& p1 = mesh.Point(el[0]);
                    const Point3d& p2 = mesh.Point(el[1]);
                    const Point3d& p3 = mesh.Point(el[2]);
                    bad2 += CalcTriangleBadness_2(p1, p2, p3, nv_z);

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

                mesh[pi1] = p1;
                mesh[pi2] = p2;

                bool should = (bad2 < bad1 && bad2 < 1e4);

                if (should) {
                    mesh[pi1] = pnew;
                    for (size_t k = 0; k < elem_idx_2.size(); k++) {
                        Element2d& el = mesh.SurfaceElement(elem_idx_2[k]);
                        if (el.IsDeleted()) continue;
                        elements_on_node.Add(pi1, elem_idx_2[k]);

                        bool haspi1 = 0;
                        for (size_t l = 0; l < 3; l++) {
                            if (el[l] == pi1)
                                haspi1 = 1;
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
                        mesh.SurfaceElement(has_both_pi[k]).Delete();
                    }
                }
            }
        }

        mesh.Compress();
        mesh.SetNextTimeStamp();
    }
}  // namespace meshit
