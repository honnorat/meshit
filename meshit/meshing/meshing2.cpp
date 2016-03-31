#include "meshing2.hpp"

#include "global.hpp"
#include "../gprim/geomtest3d.hpp"

namespace meshit
{
    Meshing2::Meshing2(const MeshingParameters& mp, const Box<3>& aboundingbox)
    {
        boundingbox = aboundingbox;

        LoadRules(NULL);
        // LoadRules ("rules/triangle.rls");

        adfront = new AdFront2(boundingbox);
        maxarea = -1;
    }

    Meshing2::~Meshing2()
    {
        delete adfront;
        for (size_t i = 0; i < rules.size(); i++) {
            delete rules[i];
        }
    }

    void Meshing2::AddPoint(const Point3d& p, PointIndex globind)
    {
        adfront->AddPoint(p, globind);
    }

    void Meshing2::AddBoundaryElement(int i1, int i2)
    {
        adfront->AddLine(i1 - 1, i2 - 1);
    }

    void Meshing2::StartMesh()
    {
        foundmap.resize(rules.size(), 0);
        canuse.resize(rules.size(), 0);
        ruleused.resize(rules.size(), 0);
    }

    void Meshing2::EndMesh()
    {
        for (size_t i = 0; i < ruleused.size(); i++) {
            MESHIT_LOG_DEBUG(std::setw(5) << ruleused[i] << " times used rule " << rules[i]->Name());
        }
    }

    void Meshing2::SetMaxArea(double amaxarea)
    {
        maxarea = amaxarea;
    }

    double Meshing2::CalcLocalH(const Point3d& /* p */, double gh) const
    {
        return gh;
    }

    void Meshing2::DefineTransformation(const Point3d& p1, const Point3d& p2)
    {
        globp1 = p1;
        ex = p2 - p1;
        ex /= ex.Length();
        ey.X() = -ex.Y();
        ey.Y() = ex.X();
        ey.Z() = 0.0;
    }

    void Meshing2::TransformToPlain(const Point3d& locpoint, Point2d& plainpoint, double h)
    {
        Vec3d p1p(globp1, locpoint);

        plainpoint.X() = (p1p * ex) / h;
        plainpoint.Y() = (p1p * ey) / h;
    }

    void Meshing2::TransformFromPlain(Point2d& plainpoint, Point3d& locpoint, double h)
    {
        locpoint.X() = globp1.X() + h * (plainpoint.X() * ex.X() + plainpoint.Y() * ey.X());
        locpoint.Y() = globp1.Y() + h * (plainpoint.X() * ex.Y() + plainpoint.Y() * ey.Y());
    }

    bool Meshing2::GenerateMesh(Mesh& mesh, const MeshingParameters& mp, double gh, int facenr)
    {
        std::vector<int> pindex, lindex;
        std::vector<int> delpoints;
        std::vector<uint32_t> dellines;
        std::vector<Element2d> locelements;

        bool found;
        int rulenr(-1);
        Point3d p1, p2;

        bool debugflag;

        double h, his, hshould;

        std::vector<Point3d> locpoints;
        std::vector<int> legalpoints;
        std::vector<Point2d> plainpoints;
        std::vector<INDEX_2> loclines;
        int cntelem = 0, trials = 0, nfaces = 0;
        int qualclass;

        // test for 3d overlaps
        Box3dTree surfeltree(static_cast<Point3d>(boundingbox.PMin()),
                             static_cast<Point3d>(boundingbox.PMax()));

        std::vector<size_t> intersecttrias;
        std::vector<Point3d> critpoints;

        // test for doubled edges

        StartMesh();

        std::vector<int> trigsonnode;
        std::vector<int> illegalpoint;

        trigsonnode.resize(mesh.GetNP(), 0);
        illegalpoint.resize(mesh.GetNP(), 0);

        double meshedarea = 0;

        std::vector<SurfaceElementIndex> seia;
        mesh.GetSurfaceElementsOfFace(facenr, seia);
        for (size_t i = 0; i < seia.size(); i++) {
            const Element2d& sel = mesh.SurfaceElement(seia[i]);

            if (sel.IsDeleted()) continue;

            Box<3> box;
            box.Set(mesh[sel[0]]);
            box.Add(mesh[sel[1]]);
            box.Add(mesh[sel[2]]);
            surfeltree.Insert(box, seia[i]);
        }

        if (maxarea > 0)
            meshedarea = mesh.SurfaceArea();

        adfront->SetStartFront();

        int plotnexttrial = 999;

        double meshedarea_before = meshedarea;

        while (!adfront->Empty()) {
            locpoints.clear();
            loclines.clear();
            pindex.clear();
            lindex.clear();
            delpoints.clear();
            dellines.clear();
            locelements.clear();

            // plot statistics
            if (trials > plotnexttrial) {
                MESHIT_LOG_DEBUG(nfaces << " faces, " << trials << " trials, " << mesh.GetNSE() << " elements.");
                plotnexttrial += 1000;
            }

            nfaces = adfront->GetNFL();
            trials++;

            if (trials % 1000 == 0) {
                for (size_t i = 0; i < canuse.size(); i++) {
                    MESHIT_LOG_DEBUG(" "
                                     << std::setw(4) << foundmap[i] << "/"
                                     << std::setw(4) << canuse[i] << "/"
                                     << std::setw(4) << ruleused[i] <<
                                     " map/can/use rule " << rules[i]->Name());
                }
            }

            int baselineindex = adfront->SelectBaseLine(p1, p2, qualclass);

            found = true;
            his = Dist(p1, p2);

            Point3d pmid = Center(p1, p2);
            hshould = CalcLocalH(pmid, mesh.GetH(pmid));
            if (gh < hshould) hshould = gh;

            mesh.RestrictLocalH(pmid, hshould);

            h = hshould;

            double hinner = (3 + qualclass) * std::max(his, hshould);

            adfront->GetLocals(baselineindex, locpoints, loclines, pindex, lindex, 2 * hinner);

            if (qualclass > mp.giveup_tol2d) {
                MESHIT_LOG_WARNING("give up with qualclass " << qualclass <<
                                   " : number of frontlines = " << adfront->GetNFL());
                break;
            }

            PointIndex gpi1 = adfront->GetGlobalIndex(pindex[loclines[0].I1() - 1]);
            PointIndex gpi2 = adfront->GetGlobalIndex(pindex[loclines[0].I2() - 1]);

            debugflag = false;

            if (debugparam.haltface && debugparam.haltfacenr == facenr) {
                debugflag = true;
                MESHIT_LOG_DEBUG("set debugflag");
            }

            if (debugparam.haltlargequalclass && qualclass > 50)
                debugflag = true;

            // problem recognition !
            if (found
                && gpi1 < static_cast<PointIndex>(illegalpoint.size())
                && gpi2 < static_cast<PointIndex>(illegalpoint.size())
                && (illegalpoint[gpi1] || illegalpoint[gpi2])) {
                found = false;
            }

            size_t oldnl = 0;
            size_t oldnp = 0;

            if (found) {
                DefineTransformation(p1, p2);

                plainpoints.resize(locpoints.size());

                if (debugflag) {
                    MESHIT_LOG_DEBUG("3d->2d transformation");
                }

                for (size_t i = 0; i < locpoints.size(); i++) {
                    TransformToPlain(locpoints[i], plainpoints[i], h);
                }

                legalpoints.assign(plainpoints.size(), 1);
                double avy = 0;
                for (size_t i = 0; i < plainpoints.size(); i++) {
                    avy += plainpoints[i].Y();
                }
                avy *= 1. / plainpoints.size();


                for (size_t i = 0; i < plainpoints.size(); i++) {
                    if (pindex[i] == -1) {
                        legalpoints[i] = 0;
                    }
                    if (plainpoints[i].Y() < -1e-10 * avy)  // changed
                    {
                        legalpoints[i] = 0;
                    }
                }

                size_t maxlegalpoint = locpoints.size();
                size_t maxlegalline = loclines.size();

                oldnl = loclines.size();
                oldnp = plainpoints.size();

                rulenr = ApplyRules(plainpoints, legalpoints, maxlegalpoint,
                                    loclines, maxlegalline, locelements,
                                    dellines, qualclass, mp);
                if (!rulenr) {
                    found = false;
                    if (debugflag || debugparam.haltnosuccess)
                        MESHIT_LOG_WARNING("no rule found");
                }
            }

            if (found) {
                for (size_t i = 0; i < locelements.size() && found; i++) {
                    const Element2d& el = locelements[i];
                    for (size_t j = 1; j <= 3; j++) {
                        if (el.PNum(j) <= static_cast<PointIndex>(oldnp) && pindex[el.PNum(j) - 1] == -1) {
                            found = false;
                            MESHIT_LOG_ERROR("meshing2, index missing");
                        }
                    }
                }
            }
            if (found) {
                locpoints.resize(plainpoints.size());

                for (size_t i = oldnp; i < plainpoints.size(); i++) {
                    TransformFromPlain(plainpoints[i], locpoints[i], h);
                }

                double violateminh = 3 + 0.1 * qualclass * qualclass;
                double minh = 1e8;
                double newedgemaxh = 0;
                for (size_t i = oldnl; i < loclines.size(); i++) {
                    double eh = Dist(locpoints[loclines[i].I1() - 1],
                                     locpoints[loclines[i].I2() - 1]);
                    if (eh > newedgemaxh)
                        newedgemaxh = eh;
                }

                for (size_t i = 0; i < locelements.size(); i++) {
                    Point3d pmin = locpoints[locelements[i].PNum(1) - 1];
                    Point3d pmax = pmin;
                    for (size_t j = 2; j <= 3; j++) {
                        const Point3d& hp =
                            locpoints[locelements[i].PNum(j) - 1];
                        pmin.SetToMin(hp);
                        pmax.SetToMax(hp);
                    }
                    double eh = mesh.GetMinH(pmin, pmax);
                    if (eh < minh)
                        minh = eh;
                }

                for (size_t i = 0; i < locelements.size(); i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        if (Dist2(locpoints[locelements[i].PNum(j) - 1], pmid) > hinner * hinner)
                            found = false;
                    }
                }

                static double maxviolate = 0;
                if (newedgemaxh / minh > maxviolate) {
                    maxviolate = newedgemaxh / minh;
                }

                if (newedgemaxh > violateminh * minh) {
                    found = 0;
                    loclines.resize(oldnl);
                    locpoints.resize(oldnp);

                    if (debugflag || debugparam.haltnosuccess)
                        MESHIT_LOG_ERROR("meshing2, maxh too large");
                }
            }

            if (found && mp.check_overlap) {

                Point3d hullmin(1e10, 1e10, 1e10);
                Point3d hullmax(-1e10, -1e10, -1e10);

                for (size_t i = 0; i < locelements.size(); i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        const Point3d& p = locpoints[locelements[i].PNum(j) - 1];
                        hullmin.SetToMin(p);
                        hullmax.SetToMax(p);
                    }
                }
                hullmin += Vec3d(-his, -his, -his);
                hullmax += Vec3d(his, his, his);

                surfeltree.GetIntersecting(hullmin, hullmax, intersecttrias);

                critpoints.resize(0);
                for (size_t i = oldnp; i < locpoints.size(); i++) {
                    critpoints.push_back(locpoints[i]);
                }
                for (size_t i = 0; i < locelements.size(); i++) {
                    const Element2d& tri = locelements[i];
                    const Point3d& tp1 = locpoints[tri.PNum(1) - 1];
                    const Point3d& tp2 = locpoints[tri.PNum(2) - 1];
                    const Point3d& tp3 = locpoints[tri.PNum(3) - 1];

                    Vec3d tv1(tp1, tp2);
                    Vec3d tv2(tp1, tp3);

                    double lam1, lam2;
                    for (lam1 = 0.2; lam1 <= 0.8; lam1 += 0.2) {
                        for (lam2 = 0.2; lam2 + lam1 <= 0.8; lam2 += 0.2) {
                            Point3d hp = tp1 + lam1 * tv1 + lam2 * tv2;
                            critpoints.push_back(hp);
                        }
                    }
                }

                for (size_t i = 0; i < critpoints.size(); i++) {
                    const Point3d& p = critpoints[i];

                    for (size_t jj = 0; jj < intersecttrias.size(); jj++) {
                        SurfaceElementIndex j = intersecttrias[jj];
                        const Element2d& el = mesh.SurfaceElement(j);

                        Point3d tp1, tp2, tp3;

                        tp1 = mesh.Point(el.PNum(1));
                        tp2 = mesh.Point(el.PNum(2));
                        tp3 = mesh.Point(el.PNum(3));

                        Vec3d e1(tp1, tp2);
                        Vec3d e2(tp1, tp3);
                        Vec3d n = Cross(e1, e2);
                        n /= n.Length();
                        double lam1, lam2, lam3;
                        lam3 = n * Vec3d(tp1, p);
                        LocalCoordinates(e1, e2, Vec3d(tp1, p), lam1, lam2);

                        if (fabs(lam3) < 0.1 * hshould && lam1 > 0 && lam2 > 0 && (lam1 + lam2) < 1.0) {
                            for (int k = 1; k <= 5; k++) {
                                adfront->IncrementClass(lindex[0]);
                            }
                            found = false;

                            if (debugflag || debugparam.haltnosuccess)
                                MESHIT_LOG_WARNING("overlapping");
                        }
                    }
                }
            }

            if (found) {
                // check, whether new front line already exists

                for (size_t i = oldnl; i < loclines.size(); i++) {
                    int nlgpi1 = loclines[i].I1();
                    int nlgpi2 = loclines[i].I2();
                    if (nlgpi1 <= static_cast<INDEX>(pindex.size()) &&
                        nlgpi2 <= static_cast<INDEX>(pindex.size())) {
                        nlgpi1 = adfront->GetGlobalIndex(pindex[nlgpi1 - 1]);
                        nlgpi2 = adfront->GetGlobalIndex(pindex[nlgpi2 - 1]);

                        int exval = adfront->ExistsLine(nlgpi1, nlgpi2);
                        if (exval) {
                            MESHIT_LOG_ERROR("ERROR: new line exits, val = " << exval);
                            found = false;
                        }
                    }
                }
            }

            if (found) {
                // everything is ok, perform mesh update

                ruleused[rulenr - 1]++;

                pindex.resize(locpoints.size());

                for (size_t i = oldnp; i < locpoints.size(); i++) {
                    PointIndex globind = mesh.AddPoint(locpoints[i]);
                    pindex[i] = adfront->AddPoint(locpoints[i], globind);
                }

                for (size_t i = oldnl; i < loclines.size(); i++) {
                    if (pindex[loclines[i].I1() - 1] == -1 ||
                        pindex[loclines[i].I2() - 1] == -1) {
                        std::cerr << "pindex is 0" << std::endl;
                    }
                    adfront->AddLine(pindex[loclines[i].I1() - 1],
                                     pindex[loclines[i].I2() - 1]);
                }
                for (size_t i = 0; i < locelements.size(); i++) {
                    Element2d mtri;
                    mtri = locelements[i];
                    mtri.SetIndex(facenr);

                    for (size_t j = 1; j <= 3; j++) {
                        mtri.PNum(j) = locelements[i].PNum(j) =
                            adfront->GetGlobalIndex(pindex[locelements[i].PNum(j) - 1]);
                    }

                    mesh.AddSurfaceElement(mtri);
                    cntelem++;

                    Box<3> box;
                    box.Set(mesh[mtri[0]]);
                    box.Add(mesh[mtri[1]]);
                    box.Add(mesh[mtri[2]]);
                    surfeltree.Insert(box, mesh.GetNSE() - 1);

                    const Point3d& sep1 = mesh.Point(mtri.PNum(1));
                    const Point3d& sep2 = mesh.Point(mtri.PNum(2));
                    const Point3d& sep3 = mesh.Point(mtri.PNum(3));

                    double trigarea = Cross(Vec3d(sep1, sep2),
                                            Vec3d(sep1, sep3)).Length() / 2;

                    meshedarea += trigarea;

                    if (maxarea > 0 && meshedarea - meshedarea_before > maxarea) {
                        std::cerr << "meshed area = " << meshedarea - meshedarea_before << std::endl
                        << "maximal area = " << maxarea << std::endl
                        << "GIVING UP" << std::endl;
                        return false;
                    }

                    for (size_t j = 1; j <= 3; j++) {

                        int gpi = locelements[i].PNum(j);
                        int oldts = trigsonnode.size();
                        if (gpi >= oldts) {
                            trigsonnode.resize(gpi + 1);
                            illegalpoint.resize(gpi + 1);
                            for (int k = oldts; k <= gpi; k++) {
                                trigsonnode[k] = 0;
                                illegalpoint[k] = 0;
                            }
                        }

                        trigsonnode[gpi]++;

                        if (trigsonnode[gpi] > 20) {
                            illegalpoint[gpi] = 1;
                            std::cerr << "illegal point: " << gpi << std::endl;
                        }

                        static int mtonnode = 0;
                        if (trigsonnode[gpi] > mtonnode)
                            mtonnode = trigsonnode[gpi];
                    }
                }

                for (size_t i = 0; i < dellines.size(); i++) {
                    adfront->DeleteLine(lindex[dellines[i] - 1]);
                }

                if (debugparam.haltsuccess || debugflag) {

                    std::cout << "success of rule" << rules[rulenr - 1]->Name() << std::endl;
                    std::cerr << "trials = " << trials << std::endl;
                    std::cerr << "locpoints " << std::endl;
                    for (size_t i = 0; i < pindex.size(); i++) {
                        std::cerr << adfront->GetGlobalIndex(pindex[i]) << std::endl;
                    }
                    std::cerr << "old number of lines = " << oldnl << std::endl;
                    for (size_t i = 0; i < loclines.size(); i++) {
                        std::cerr << "line ";
                        for (size_t j = 1; j <= 2; j++) {
                            int hi = 0;
                            if (loclines[i].I(j) >= 1 &&
                                loclines[i].I(j) <= static_cast<INDEX>(pindex.size())) {
                                hi = adfront->GetGlobalIndex(pindex[loclines[i].I(j) - 1]);
                            }
                            std::cerr << hi << " ";
                        }
                        std::cerr << " : "
                        << plainpoints[loclines[i].I1() - 1] << " - "
                        << plainpoints[loclines[i].I2() - 1] << " 3d: "
                        << locpoints[loclines[i].I1() - 1] << " - "
                        << locpoints[loclines[i].I2() - 1] << std::endl;
                    }
                }
            } else {
                adfront->IncrementClass(lindex[0]);

                if (debugparam.haltnosuccess || debugflag) {
                    std::cerr << "Problem with seg " << gpi1 << " - " << gpi2 << ", class = " << qualclass << std::endl;

                    for (size_t i = 0; i < loclines.size(); i++) {
                        std::cerr << "line ";
                        for (size_t j = 1; j <= 2; j++) {
                            int hi = 0;
                            if (loclines[i].I(j) >= 1 &&
                                loclines[i].I(j) <= static_cast<INDEX>(pindex.size()))
                                hi = adfront->GetGlobalIndex(pindex[loclines[i].I(j) - 1]);

                            std::cerr << hi << " ";
                        }
                        std::cerr << " : "
                        << plainpoints[loclines[i].I1() - 1] << " - "
                        << plainpoints[loclines[i].I2() - 1] << " 3d: "
                        << locpoints[loclines[i].I1() - 1] << " - "
                        << locpoints[loclines[i].I2() - 1] << std::endl;
                    }
                }
            }
        }

        MESHIT_LOG_DEBUG("Surface meshing done");

        adfront->PrintOpenSegments(std::cout);

        EndMesh();

        return (!adfront->Empty());
    }

}  // namespace meshit
