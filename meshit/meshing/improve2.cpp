#include <meshit.hpp>
#include "improve2.hpp"
#include "global.hpp"
#include "../general/ngexception.hpp"
#include "../gprim/geomfuncs.hpp"

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

    void MeshOptimize2d::EdgeSwapping(Mesh & mesh, int usemetric)
    {
        if (!faceindex) {
            if (usemetric)
                PrintMessage(3, "Edgeswapping, metric");
            else
                PrintMessage(3, "Edgeswapping, topological");

            for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++) {
                EdgeSwapping(mesh, usemetric);

                if (multithread.terminate)
                    throw NgException("Meshing stopped");
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
            const Element2d & sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++) {
                pangle[sel[j]] = 0.0;
            }
        }
        // pangle = 0;

        for (int i = 0; i < seia.size(); i++) {
            const Element2d & sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++) {
                POINTTYPE typ = mesh[sel[j]].Type();
                if (typ == FIXEDPOINT || typ == EDGEPOINT) {
                    pangle[sel[j]] +=
                            Angle(mesh[sel[(j + 1) % 3]] - mesh[sel[j]],
                            mesh[sel[(j + 2) % 3]] - mesh[sel[j]]);
                }
            }
        }

        for (int i = 0; i < seia.size(); i++) {
            const Element2d & sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < 3; j++) {
                PointIndex pi = sel[j];
                if (mesh[pi].Type() == INNERPOINT || mesh[pi].Type() == SURFACEPOINT)
                    pdef[pi] = -6;
                else {
                    for (int j = 0; j < 8; j++) {
                        if (pangle[pi] >= minangle[j]) {
                            pdef[pi] = -1 - j;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < seia.size(); i++) {
            const Element2d & sel = mesh.SurfaceElement(seia[i]);
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
            const Element2d & sel = mesh.SurfaceElement(seia[i]);

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
                }
                else {
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

                if (multithread.terminate)
                    throw NgException("Meshing stopped");

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

                    Vec<3> auxvec1 = mesh[pi3] - mesh[pi4];
                    Vec<3> auxvec2 = mesh[pi1] - mesh[pi4];

                    allowswap = allowswap && fabs(1. - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

                    if (!allowswap)
                        continue;

                    // normal of new
                    Vec<3> nv1 = Cross(auxvec1, auxvec2);

                    auxvec1 = mesh.Point(pi4) - mesh.Point(pi3);
                    auxvec2 = mesh.Point(pi2) - mesh.Point(pi3);
                    allowswap = allowswap && fabs(1. - (auxvec1 * auxvec2) / (auxvec1.Length() * auxvec2.Length())) > 1e-4;

                    if (!allowswap)
                        continue;

                    Vec<3> nv2 = Cross(auxvec1, auxvec2);

                    // normals of original
                    Vec<3> nv3 = Cross(mesh[pi1] - mesh[pi4], mesh[pi2] - mesh[pi4]);
                    Vec<3> nv4 = Cross(mesh[pi2] - mesh[pi3], mesh[pi1] - mesh[pi3]);

                    nv3 *= -1;
                    nv4 *= -1;
                    nv3.Normalize();
                    nv4.Normalize();

                    nv1.Normalize();
                    nv2.Normalize();

//                    Vec<3> nvp3 = Vec<3> (0, 0, 1);
//                    Vec<3> nvp4 = Vec<3> (0, 0, 1);

                    double critval = cos(M_PI / 6); // 30 degree
//                    allowswap = allowswap &&
//                            (nv1 * nvp3 > critval) &&
//                            (nv1 * nvp4 > critval) &&
//                            (nv2 * nvp3 > critval) &&
//                            (nv2 * nvp4 > critval) &&
//                            (nvp3 * nv3 > critval) &&
//                            (nvp4 * nv4 > critval);
                    allowswap = allowswap &&
                            (nv1(2) > critval) &&
                            (nv2(2) > critval) &&
                            (nv3(2) > critval) &&
                            (nv4(2) > critval);

                    double horder = Dist(mesh.Point(pi1), mesh.Point(pi2));

                    if (// nv1 * nv2 >= 0 &&
                            nv1.Length() > 1e-3 * horder * horder &&
                            nv2.Length() > 1e-3 * horder * horder &&
                            allowswap) {
                        if (!usemetric) {
                            int e = pdef[pi1] + pdef[pi2] - pdef[pi3] - pdef[pi4];
                            double d =
                                    Dist2(mesh.Point(pi1), mesh.Point(pi2)) -
                                    Dist2(mesh.Point(pi3), mesh.Point(pi4));

                            should = e >= t && (e > 2 || d > 0);
                        }
                        else {
                            double loch = mesh.GetH(mesh[pi1]);
                            should =
                                    CalcTriangleBadness(mesh.Point(pi4), mesh.Point(pi3), mesh.Point(pi1),
                                    metricweight, loch) +
                                    CalcTriangleBadness(mesh.Point(pi3), mesh.Point(pi4), mesh.Point(pi2),
                                    metricweight, loch) <
                                    CalcTriangleBadness(mesh.Point(pi1), mesh.Point(pi2), mesh.Point(pi3),
                                    metricweight, loch) +
                                    CalcTriangleBadness(mesh.Point(pi2), mesh.Point(pi1), mesh.Point(pi4),
                                    metricweight, loch);
                        }

                        if (allowswap) {
                            Element2d sw1(pi4, pi3, pi1);
                            Element2d sw2(pi3, pi4, pi2);

                            int legal1 =
                                    mesh.LegalTrig(mesh.SurfaceElement(t1)) +
                                    mesh.LegalTrig(mesh.SurfaceElement(t2));
                            int legal2 =
                                    mesh.LegalTrig(sw1) + mesh.LegalTrig(sw2);

                            if (legal1 < legal2) should = 1;
                            if (legal2 < legal1) should = 0;
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

    void MeshOptimize2d::CombineImprove(Mesh & mesh)
    {
        if (!faceindex) {
            PrintMessage(3, "Combine improve");

            for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++) {
                CombineImprove(mesh);

                if (multithread.terminate)
                    throw NgException("Meshing stopped");
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
        Vec<3> nv;

        int np = mesh.GetNP();
        //int nse = mesh.GetNSE();

        TABLE<SurfaceElementIndex, PointIndex::BASE> elementsonnode(np);
        Array<SurfaceElementIndex> hasonepi, hasbothpi;

        for (int i = 0; i < seia.size(); i++) {
            Element2d & el = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < el.GetNP(); j++)
                elementsonnode.Add(el[j], seia[i]);
        }

        Array<bool, PointIndex::BASE > fixed(np);
        fixed = false;

        for (int i = 0; i < seia.size(); i++) {
            Element2d & sel = mesh.SurfaceElement(seia[i]);
            for (int j = 0; j < sel.GetNP(); j++) {
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

        Array<Vec<3>, PointIndex::BASE> normals(np);

        for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++) {
            if (elementsonnode[pi].size()) {
                Element2d & hel = mesh.SurfaceElement(elementsonnode[pi][0]);
                for (int k = 0; k < 3; k++) {
                    if (hel[k] == pi) {
                        GetNormalVector(surfnr, mesh[pi], hel.GeomInfoPi(k + 1), normals[pi]);
                        break;
                    }
                }
            }
        }

        for (int i = 0; i < seia.size(); i++) {
            SurfaceElementIndex sei = seia[i];
            Element2d & elem = mesh.SurfaceElement(sei);
            if (elem.IsDeleted()) continue;

            for (int j = 0; j < 3; j++) {
                PointIndex pi1 = elem[j];
                PointIndex pi2 = elem[(j + 1) % 3];

                if (pi1 < PointIndex::BASE || pi2 < PointIndex::BASE)
                    continue;

                bool debugflag = 0;

                if (debugflag) {
                    std::cerr << "Combineimprove, face = " << faceindex
                            << "pi1 = " << pi1 << " pi2 = " << pi2 << std::endl;
                }

                // more general 
                if (fixed[pi2])
                    std::swap(pi1, pi2);

                if (fixed[pi2])
                    continue;

                double loch = mesh.GetH(mesh[pi1]);

                INDEX_2 si2(pi1, pi2);
                si2.Sort();

                hasonepi.resize(0);
                hasbothpi.resize(0);

                for (int k = 0; k < elementsonnode[pi1].size(); k++) {
                    const Element2d & el2 = mesh.SurfaceElement(elementsonnode[pi1][k]);

                    if (el2.IsDeleted()) continue;

                    if (el2[0] == pi2 || el2[1] == pi2 || el2[2] == pi2) {
                        hasbothpi.push_back(elementsonnode[pi1][k]);
                        nv = Cross(Vec3d(mesh[el2[0]], mesh[el2[1]]),
                                Vec3d(mesh[el2[0]], mesh[el2[2]]));
                    }
                    else {
                        hasonepi.push_back(elementsonnode[pi1][k]);
                    }
                }

                Element2d & hel = mesh.SurfaceElement(hasbothpi[0]);
                for (int k = 0; k < 3; k++)
                    if (hel[k] == pi1) {
                        GetNormalVector(surfnr, mesh[pi1], hel.GeomInfoPi(k + 1), nv);
                        break;
                    }

                for (int k = 0; k < elementsonnode[pi2].size(); k++) {
                    const Element2d & el2 = mesh.SurfaceElement(elementsonnode[pi2][k]);
                    if (el2.IsDeleted()) continue;

                    if (el2[0] == pi1 || el2[1] == pi1 || el2[2] == pi1)
                        ;
                    else
                        hasonepi.push_back(elementsonnode[pi2][k]);
                }

                bad1 = 0;
                int illegal1 = 0, illegal2 = 0;
                for (int k = 0; k < hasonepi.size(); k++) {
                    const Element2d & el = mesh.SurfaceElement(hasonepi[k]);
                    bad1 += CalcTriangleBadness(mesh[el[0]], mesh[el[1]], mesh[el[2]],
                            nv, -1, loch);
                    illegal1 += 1 - mesh.LegalTrig(el);
                }

                for (int k = 0; k < hasbothpi.size(); k++) {
                    const Element2d & el = mesh.SurfaceElement(hasbothpi[k]);
                    bad1 += CalcTriangleBadness(mesh[el[0]], mesh[el[1]], mesh[el[2]],
                            nv, -1, loch);
                    illegal1 += 1 - mesh.LegalTrig(el);
                }
                bad1 /= (hasonepi.size() + hasbothpi.size());

                MeshPoint p1 = mesh[pi1];
                MeshPoint p2 = mesh[pi2];

                MeshPoint pnew = p1;
                mesh[pi1] = pnew;
                mesh[pi2] = pnew;

                bad2 = 0;
                for (int k = 0; k < hasonepi.size(); k++) {
                    Element2d & el = mesh.SurfaceElement(hasonepi[k]);
                    double err =
                            CalcTriangleBadness(mesh[el[0]], mesh[el[1]], mesh[el[2]],
                            nv, -1, loch);
                    bad2 += err;

                    Vec<3> hnv = Cross(Vec3d(mesh[el[0]],
                            mesh[el[1]]),
                            Vec3d(mesh[el[0]],
                            mesh[el[2]]));
                    if (hnv * nv < 0)
                        bad2 += 1e10;

                    for (int l = 0; l < 3; l++)
                        if ((normals[el[l]] * nv) < 0.5)
                            bad2 += 1e10;

                    illegal2 += 1 - mesh.LegalTrig(el);
                }
                bad2 /= hasonepi.size();

                mesh[pi1] = p1;
                mesh[pi2] = p2;

                if (debugflag) {
                    std::cerr << "bad1 = " << bad1 << ", bad2 = " << bad2 << std::endl;
                }

                bool should = (bad2 < bad1 && bad2 < 1e4);
                if (bad2 < 1e4) {
                    if (illegal1 > illegal2) should = 1;
                    if (illegal2 > illegal1) should = 0;
                }

                if (should) {
                    
                    mesh[pi1] = pnew;
                    PointGeomInfo gi;

                    Element2d * el1p(NULL);
                    int l = 0;
                    while (mesh.SurfaceElement(elementsonnode[pi1][l]).IsDeleted() && l < elementsonnode.EntrySize(pi1))
                        l++;
                    if (l < elementsonnode.EntrySize(pi1))
                        el1p = &mesh.SurfaceElement(elementsonnode[pi1][l]);
                    else
                        std::cerr << "OOPS!" << std::endl;

                    for (l = 0; l < el1p->GetNP(); l++)
                        if ((*el1p)[l] == pi1) {
                            gi = el1p->GeomInfoPi(l + 1);
                        }

                    // std::cerr << "Connect point " << pi2 << " to " << pi1 << "\n";
                    for (int k = 0; k < elementsonnode[pi2].size(); k++) {
                        Element2d & el = mesh.SurfaceElement(elementsonnode[pi2][k]);
                        if (el.IsDeleted()) continue;
                        elementsonnode.Add(pi1, elementsonnode[pi2][k]);

                        bool haspi1 = 0;
                        for (l = 0; l < el.GetNP(); l++)
                            if (el[l] == pi1)
                                haspi1 = 1;
                        if (haspi1) continue;

                        for (int l = 0; l < el.GetNP(); l++) {
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
}