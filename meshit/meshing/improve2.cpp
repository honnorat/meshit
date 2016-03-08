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

        void SetOrientation(int side, int aorient)
        {
            orient[side] = aorient;
        }

        int GetOrientation(int side)
        {
            return orient[side];
        }
    };

    class trionedge
    {
     public:
        int tnr;
        int sidenr;

        trionedge()
        {
            tnr = 0;
            sidenr = 0;
        }

        trionedge(int atnr, int asidenr)
        {
            tnr = atnr;
            sidenr = asidenr;
        }
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

        Array<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace(faceindex, seia);

        for (int i = 0; i < seia.size(); i++) {
            if (mesh.SurfaceElement(seia[i]).GetNP() != 3) {
                GenericImprove(mesh);
                return;
            }
        }

        Array<Neighbour> neighbors(mesh.GetNSE());
        INDEX_2_HASHTABLE<trionedge> other(seia.size() + 2);


        Array<char> swapped(mesh.GetNSE());
        Array<int, PointIndex::BASE> pdef(mesh.GetNP());
        Array<double, PointIndex::BASE> pangle(mesh.GetNP());

        static const double minangle[] = {0, 1.481, 2.565, 3.627, 4.683, 5.736, 7, 9};

        for (int i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++) {
                pangle[sel[j]] = 0.0;
            }
        }
        // pangle = 0;

        for (int i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++) {
                POINTTYPE typ = mesh.Point(sel[j]).Type();
                if (typ == FIXEDPOINT || typ == EDGEPOINT) {
                    pangle[sel[j]] += Angle(
                            mesh.Point(sel[(j + 1) % 3]) - mesh.Point(sel[j]),
                            mesh.Point(sel[(j + 2) % 3]) - mesh.Point(sel[j]));
                }
            }
        }

        for (int i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++) {
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

        for (int i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++)
                pdef[sel[j]]++;
        }

        for (int i = 0; i < seia.size(); i++) {
            for (int j = 0; j < 3; j++) {
                neighbors[seia[i]].SetNr(j, -1);
                neighbors[seia[i]].SetOrientation(j, 0);
            }
        }

        for (int i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);

            for (int j = 0; j < 3; j++) {
                PointIndex pi1 = sel.PNumMod(j + 2);
                PointIndex pi2 = sel.PNumMod(j + 3);

                INDEX_2 edge(pi1, pi2);
                edge.Sort();

                if (mesh.IsSegment(pi1, pi2))
                    continue;

                INDEX_2 ii2(pi1, pi2);
                if (other.Used(ii2)) {
                    int i2 = other.Get(ii2).tnr;
                    int j2 = other.Get(ii2).sidenr;
                    neighbors[seia[i]].SetNr(j, i2);
                    neighbors[seia[i]].SetOrientation(j, j2);
                    neighbors[i2].SetNr(j2, seia[i]);
                    neighbors[i2].SetOrientation(j2, j);
                } else {
                    other.Set(INDEX_2(pi2, pi1), trionedge(seia[i], j));
                }
            }
        }

        for (int i = 0; i < seia.size(); i++)
            swapped[seia[i]] = 0;

        int t = 4;
        int done = 0;
        while (!done && t >= 2) {
            for (int i = 0; i < seia.size(); i++) {
                SurfaceElementIndex t1 = seia[i];

                if (mesh.SurfaceElement(t1).IsDeleted())
                    continue;

                if (mesh.SurfaceElement(t1).GetIndex() != faceindex)
                    continue;

                for (int o1 = 0; o1 < 3; o1++) {
                    bool should;

                    SurfaceElementIndex t2 = neighbors[t1].GetNr(o1);
                    int o2 = neighbors[t1].GetOrientation(o1);

                    if (t2 == -1) continue;
                    if (swapped[t1] || swapped[t2]) continue;

                    PointIndex pi1 = mesh.SurfaceElement(t1).PNumMod(o1 + 1 + 1);
                    PointIndex pi2 = mesh.SurfaceElement(t1).PNumMod(o1 + 1 + 2);
                    PointIndex pi3 = mesh.SurfaceElement(t1).PNumMod(o1 + 1);
                    PointIndex pi4 = mesh.SurfaceElement(t2).PNumMod(o2 + 1);

                    PointGeomInfo gi1 = mesh.SurfaceElement(t1).GeomInfoPiMod(o1 + 1 + 1);
                    PointGeomInfo gi2 = mesh.SurfaceElement(t1).GeomInfoPiMod(o1 + 1 + 2);
                    PointGeomInfo gi3 = mesh.SurfaceElement(t1).GeomInfoPiMod(o1 + 1);
                    PointGeomInfo gi4 = mesh.SurfaceElement(t2).GeomInfoPiMod(o2 + 1);

                    bool allowswap = true;

                    Vec3d auxvec1 = mesh.Point(pi3) - mesh.Point(pi4);
                    Vec3d auxvec2 = mesh.Point(pi1) - mesh.Point(pi4);

                    allowswap =
                            allowswap && fabs(1. - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

                    if (!allowswap)
                        continue;

                    // normal of new
                    Vec3d nv1 = Cross(auxvec1, auxvec2);

                    auxvec1 = mesh.Point(pi4) - mesh.Point(pi3);
                    auxvec2 = mesh.Point(pi2) - mesh.Point(pi3);
                    allowswap =
                            allowswap && fabs(1. - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

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
                    allowswap = allowswap &&
                                (nv1.Z() > critval) &&
                                (nv2.Z() > critval) &&
                                (nv3.Z() > critval) &&
                                (nv4.Z() > critval);

                    double horder = Dist(mesh.Point(pi1), mesh.Point(pi2));

                    if (  // nv1 * nv2 >= 0 &&
                            nv1.Length() > 1e-3 * horder * horder &&
                            nv2.Length() > 1e-3 * horder * horder &&
                            allowswap) {
                        if (!usemetric) {
                            int e = pdef[pi1] + pdef[pi2] - pdef[pi3] - pdef[pi4];
                            double d =
                                    Dist2(mesh.Point(pi1), mesh.Point(pi2)) -
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

                            done = 1;

                            mesh.SurfaceElement(t1).PNum(1) = pi1;
                            mesh.SurfaceElement(t1).PNum(2) = pi4;
                            mesh.SurfaceElement(t1).PNum(3) = pi3;

                            mesh.SurfaceElement(t2).PNum(1) = pi2;
                            mesh.SurfaceElement(t2).PNum(2) = pi3;
                            mesh.SurfaceElement(t2).PNum(3) = pi4;

                            mesh.SurfaceElement(t1).GeomInfoPi(1) = gi1;
                            mesh.SurfaceElement(t1).GeomInfoPi(2) = gi4;
                            mesh.SurfaceElement(t1).GeomInfoPi(3) = gi3;

                            mesh.SurfaceElement(t2).GeomInfoPi(1) = gi2;
                            mesh.SurfaceElement(t2).GeomInfoPi(2) = gi3;
                            mesh.SurfaceElement(t2).GeomInfoPi(3) = gi4;

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

        Array<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace(faceindex, seia);

        for (int i = 0; i < seia.size(); i++) {
            if (mesh.SurfaceElement(seia[i]).GetNP() != 3) {
                std::cerr << "exit from CombineImprove " << i << std::endl;
                return;
            }
        }

        int surfnr = 0;
        if (faceindex)
            surfnr = mesh.GetFaceDescriptor(faceindex).SurfNr();

        double bad1, bad2;
        Vec3d nv;

        int np = mesh.GetNP();

        TABLE<SurfaceElementIndex, PointIndex::BASE> elementsonnode(np);
        Array<SurfaceElementIndex> hasonepi, hasbothpi;

        for (int i = 0; i < seia.size(); i++) {
            Element2d& el = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < el.GetNP(); j++)
                elementsonnode.Add(el[j], seia[i]);
        }

        Array<bool, PointIndex::BASE> fixed(np);
        fixed = false;

        for (int i = 0; i < seia.size(); i++) {
            Element2d& sel = mesh.SurfaceElement(seia[i]);
            for (size_t j = 0; j < sel.GetNP(); j++) {
                PointIndex pi1 = sel.PNumMod(j + 2);
                PointIndex pi2 = sel.PNumMod(j + 3);
                if (mesh.IsSegment(pi1, pi2)) {
                    fixed[pi1] = true;
                    fixed[pi2] = true;
                }
            }
        }

        for (int i = 0; i < mesh.LockedPoints().size(); i++)
            fixed[mesh.LockedPoints()[i]] = true;

        Array<Vec3d, PointIndex::BASE> normals(np);

        for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++) {
            if (elementsonnode[pi].size()) {
                Element2d& hel = mesh.SurfaceElement(elementsonnode[pi][0]);
                for (int k = 0; k < 3; k++) {
                    if (hel[k] == pi) {
                        GetNormalVector(surfnr, mesh.Point(pi), hel.GeomInfoPi(k + 1), normals[pi]);
                        break;
                    }
                }
            }
        }

        for (int i = 0; i < seia.size(); i++) {
            SurfaceElementIndex sei = seia[i];
            Element2d& elem = mesh.SurfaceElement(sei);
            if (elem.IsDeleted()) continue;

            for (int j = 0; j < 3; j++) {
                PointIndex pi1 = elem[j];
                PointIndex pi2 = elem[(j + 1) % 3];

                if (pi1 < PointIndex::BASE || pi2 < PointIndex::BASE)
                    continue;

                // more general
                if (fixed[pi2])
                    std::swap(pi1, pi2);

                if (fixed[pi2])
                    continue;

                double loch = mesh.GetH(mesh.Point(pi1));

                INDEX_2 si2(pi1, pi2);
                si2.Sort();

                hasonepi.resize(0);
                hasbothpi.resize(0);

                for (int k = 0; k < elementsonnode[pi1].size(); k++) {
                    const Element2d& el2 = mesh.SurfaceElement(elementsonnode[pi1][k]);

                    if (el2.IsDeleted()) continue;

                    if (el2[0] == pi2 || el2[1] == pi2 || el2[2] == pi2) {
                        hasbothpi.push_back(elementsonnode[pi1][k]);
                        nv = Cross(
                                Vec3d(mesh.Point(el2[0]), mesh.Point(el2[1])),
                                Vec3d(mesh.Point(el2[0]), mesh.Point(el2[2])));
                    } else {
                        hasonepi.push_back(elementsonnode[pi1][k]);
                    }
                }

                Element2d& hel = mesh.SurfaceElement(hasbothpi[0]);
                for (int k = 0; k < 3; k++)
                    if (hel[k] == pi1) {
                        GetNormalVector(surfnr, mesh.Point(pi1), hel.GeomInfoPi(k + 1), nv);
                        break;
                    }

                for (int k = 0; k < elementsonnode[pi2].size(); k++) {
                    const Element2d& el2 = mesh.SurfaceElement(elementsonnode[pi2][k]);
                    if (el2.IsDeleted()) continue;

                    if (el2[0] == pi1 || el2[1] == pi1 || el2[2] == pi1);
                    else
                        hasonepi.push_back(elementsonnode[pi2][k]);
                }

                bad1 = 0;
                for (int k = 0; k < hasonepi.size(); k++) {
                    const Element2d& el = mesh.SurfaceElement(hasonepi[k]);
                    bad1 += CalcTriangleBadness(
                            mesh.Point(el[0]), mesh.Point(el[1]), mesh.Point(el[2]), nv, -1, loch);
                }

                for (int k = 0; k < hasbothpi.size(); k++) {
                    const Element2d& el = mesh.SurfaceElement(hasbothpi[k]);
                    bad1 += CalcTriangleBadness(
                            mesh.Point(el[0]), mesh.Point(el[1]), mesh.Point(el[2]), nv, -1, loch);
                }
                bad1 /= (hasonepi.size() + hasbothpi.size());

                MeshPoint p1 = mesh[pi1];
                MeshPoint p2 = mesh[pi2];

                MeshPoint pnew = p1;
                mesh[pi1] = pnew;
                mesh[pi2] = pnew;

                bad2 = 0;
                for (int k = 0; k < hasonepi.size(); k++) {
                    Element2d& el = mesh.SurfaceElement(hasonepi[k]);
                    double err = CalcTriangleBadness(
                            mesh.Point(el[0]), mesh.Point(el[1]), mesh.Point(el[2]), nv, -1, loch);
                    bad2 += err;

                    Vec3d hnv = Cross(
                            mesh.Point(el[1]) - mesh.Point(el[0]),
                            mesh.Point(el[2]) - mesh.Point(el[0]));
                    if (hnv * nv < 0)
                        bad2 += 1e10;

                    for (int l = 0; l < 3; l++)
                        if ((normals[el[l]] * nv) < 0.5)
                            bad2 += 1e10;
                }
                bad2 /= hasonepi.size();

                mesh[pi1] = p1;
                mesh[pi2] = p2;

                bool should = (bad2 < bad1 && bad2 < 1e4);

                if (should) {
                    mesh[pi1] = pnew;
                    PointGeomInfo gi;

                    Element2d* el1p(NULL);
                    int l = 0;
                    while (mesh.SurfaceElement(elementsonnode[pi1][l]).IsDeleted() && l < elementsonnode.EntrySize(pi1))
                        l++;
                    if (l < elementsonnode.EntrySize(pi1))
                        el1p = &mesh.SurfaceElement(elementsonnode[pi1][l]);
                    else
                        std::cerr << "OOPS!" << std::endl;

                    for (size_t l = 0; l < el1p->GetNP(); l++) {
                        if ((*el1p)[l] == pi1) {
                            gi = el1p->GeomInfoPi(l + 1);
                        }
                    }
                    for (int k = 0; k < elementsonnode[pi2].size(); k++) {
                        Element2d& el = mesh.SurfaceElement(elementsonnode[pi2][k]);
                        if (el.IsDeleted()) continue;
                        elementsonnode.Add(pi1, elementsonnode[pi2][k]);

                        bool haspi1 = 0;
                        for (size_t l = 0; l < el.GetNP(); l++) {
                            if (el[l] == pi1)
                                haspi1 = 1;
                        }
                        if (haspi1) continue;

                        for (size_t l = 0; l < el.GetNP(); l++) {
                            if (el[l] == pi2) {
                                el[l] = pi1;
                                el.GeomInfoPi(l + 1) = gi;
                            }
                            fixed[el[l]] = true;
                        }
                    }
                    for (int k = 0; k < hasbothpi.size(); k++) {
                        mesh.SurfaceElement(hasbothpi[k]).Delete();
                    }
                }
            }
        }

        mesh.Compress();
        mesh.SetNextTimeStamp();
    }
}  // namespace meshit
