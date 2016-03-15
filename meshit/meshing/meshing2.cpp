#include "meshing2.hpp"

#include "global.hpp"
#include "../gprim/geomtest3d.hpp"

namespace meshit {

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
        foundmap.resize(rules.size());
        canuse.resize(rules.size());
        ruleused.resize(rules.size());

        foundmap = 0;
        canuse = 0;
        ruleused = 0;
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
        ey.Z() = 0;
    }

    void Meshing2::TransformToPlain(const Point3d& locpoint, Point2d& plainpoint, double h, int& zone)
    {
        Vec3d p1p(globp1, locpoint);

        p1p /= h;
        plainpoint.X() = p1p * ex;
        plainpoint.Y() = p1p * ey;
        zone = 0;
    }

    int Meshing2::TransformFromPlain(
            Point2d& plainpoint,
            Point3d& locpoint,
            double h)
    {
        Vec3d p1p;

        p1p = plainpoint.X() * ex + plainpoint.Y() * ey;
        p1p *= h;
        locpoint = globp1 + p1p;
        return 0;
    }

    int Meshing2::BelongsToActiveChart(const Point3d& p)
    {
        return 1;
    }

    int Meshing2::ComputePointGeomInfo(const Point3d& p)
    {
        return 0;
    }

    int Meshing2::IsLineVertexOnChart(const Point3d& p1, const Point3d& p2)
    {
        return 1;
    }

    void Meshing2::GetChartBoundary(Array<Point2d>& points, Array<Point3d>& points3d, Array<INDEX_2>& lines) const
    {
        points.resize(0);
        points3d.resize(0);
        lines.resize(0);
    }

    double Meshing2::Area() const
    {
        return -1;
    }

    MESHING2_RESULT Meshing2::GenerateMesh(Mesh& mesh, const MeshingParameters& mp, double gh, int facenr)
    {
        Array<int> pindex, lindex;
        Array<int> delpoints, dellines;

        Array<Element2d> locelements;

        int z1, z2, oldnp(-1);
        bool found;
        int rulenr(-1);
        Point3d p1, p2;

        bool debugflag;

        double h, his, hshould;

        Array<Point3d> locpoints;
        Array<int> legalpoints;
        Array<Point2d> plainpoints;
        Array<int> plainzones;
        Array<INDEX_2> loclines;
        int cntelem = 0, trials = 0, nfaces = 0;
        int oldnl = 0;
        int qualclass;

        // test for 3d overlaps
        Box3dTree surfeltree(boundingbox.PMin(), boundingbox.PMax());

        std::vector<size_t> intersecttrias;
        Array<Point3d> critpoints;

        // test for doubled edges

        StartMesh();

        Array<Point2d> chartboundpoints;
        Array<Point3d> chartboundpoints3d;
        Array<INDEX_2> chartboundlines;

        // illegal points: points with more then 50 elements per node
        int maxlegalpoint(-1), maxlegalline(-1);
        Array<int> trigsonnode;
        Array<int> illegalpoint;

        trigsonnode.resize(mesh.GetNP());
        illegalpoint.resize(mesh.GetNP());

        trigsonnode = 0;
        illegalpoint = 0;

        double totalarea = Area();
        double meshedarea = 0;

        Array<SurfaceElementIndex> seia;
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

        if (totalarea > 0 || maxarea > 0)
            meshedarea = mesh.SurfaceArea();

        adfront->SetStartFront();

        int plotnexttrial = 999;

        double meshedarea_before = meshedarea;

        while (!adfront->Empty()) {

            locpoints.resize(0);
            loclines.resize(0);
            pindex.resize(0);
            lindex.resize(0);
            delpoints.resize(0);
            dellines.resize(0);
            locelements.resize(0);

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

            if (qualclass > mp.giveuptol2d) {
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
            if (found && (gpi1 < illegalpoint.size()) && (gpi2 < illegalpoint.size())) {
                if (illegalpoint[gpi1] || illegalpoint[gpi2]) {
                    found = false;
                }
            }

            Point2d p12d, p22d;

            if (found) {
                oldnp = locpoints.size();
                oldnl = loclines.size();

                if (debugflag)
                    std::cerr << "define new transformation" << std::endl;

                DefineTransformation(p1, p2);

                plainpoints.resize(locpoints.size());
                plainzones.resize(locpoints.size());

                if (debugflag) {
                    MESHIT_LOG_DEBUG("3d->2d transformation");
                }

                for (int i = 0; i < locpoints.size(); i++) {
                    TransformToPlain(locpoints[i], plainpoints[i], h, plainzones[i]);
                }

                p12d = plainpoints[0];
                p22d = plainpoints[1];

                for (int i = 1; i < loclines.size(); i++) // don't remove first line
                {
                    z1 = plainzones[loclines[i].I1() - 1];
                    z2 = plainzones[loclines[i].I2() - 1];

                    // one inner point, one outer
                    if ((z1 >= 0) != (z2 >= 0)) {
                        int innerp = (z1 >= 0) ? 1 : 2;
                        if (IsLineVertexOnChart(locpoints[loclines[i].I1() - 1],
                                                locpoints[loclines[i].I2() - 1])) {

                            // use one end of line
                            int pini, pouti;
                            Vec2d v;

                            pini = loclines[i].I(innerp);
                            pouti = loclines[i].I(3 - innerp);

                            Point2d pin(plainpoints[pini - 1]);
                            Point2d pout(plainpoints[pouti - 1]);
                            v = pout - pin;
                            double len = v.Length();
                            if (len <= 1e-6)
                                std::cerr << "WARNING(js): inner-outer: short vector" << std::endl;
                            else
                                v /= len;

                            Point2d newpout = pin + 1000 * v;
                            newpout = pout;

                            plainpoints.push_back(newpout);
                            Point3d pout3d = locpoints[pouti - 1];
                            locpoints.push_back(pout3d);

                            plainzones.push_back(0);
                            pindex.push_back(-1);
                            oldnp++;
                            loclines[i - 1].I(3 - innerp) = oldnp;
                        } else {
                            // remove line
                            loclines.Delete(i - 1);
                            lindex.Delete(i - 1);
                            oldnl--;
                            i--;
                        }
                    } else if ((z1 > 0 && z2 > 0 && (z1 != z2)) || ((z1 < 0) && (z2 < 0))) {
                        loclines.Delete(i - 1);
                        lindex.Delete(i - 1);
                        oldnl--;
                        i--;
                    }
                }

                legalpoints.resize(plainpoints.size());
                for (int i = 0; i < legalpoints.size(); i++) {
                    legalpoints[i] = 1;
                }
                double avy = 0;
                for (int i = 0; i < plainpoints.size(); i++) {
                    avy += plainpoints[i].Y();
                }
                avy *= 1. / plainpoints.size();


                for (int i = 0; i < plainpoints.size(); i++) {
                    if (plainzones[i] < 0) {
                        plainpoints[i] = Point2d(1e4, 1e4);
                        legalpoints[i] = 0;
                    }
                    if (pindex[i] == -1) {
                        legalpoints[i] = 0;
                    }
                    if (plainpoints[i].Y() < -1e-10 * avy) // changed
                    {
                        legalpoints[i] = 0;
                    }
                }

                GetChartBoundary(chartboundpoints, chartboundpoints3d, chartboundlines);

                oldnp = plainpoints.size();
                maxlegalpoint = locpoints.size();
                maxlegalline = loclines.size();

                if (mp.checkchartboundary) {
                    for (int i = 0; i < chartboundpoints.size(); i++) {
                        plainpoints.push_back(chartboundpoints[i]);
                        locpoints.push_back(chartboundpoints3d[i]);
                        legalpoints.push_back(0);
                    }

                    for (int i = 0; i < chartboundlines.size(); i++) {
                        INDEX_2 line(chartboundlines[i].I1() + oldnp,
                                     chartboundlines[i].I2() + oldnp);
                        loclines.push_back(line);
                    }
                }

                oldnl = loclines.size();
                oldnp = plainpoints.size();
            }

            if (found) {
                rulenr = ApplyRules(plainpoints, legalpoints, maxlegalpoint,
                                    loclines, maxlegalline, locelements,
                                    dellines, qualclass, mp);
                if (!rulenr) {
                    found = false;
                    if (debugflag || debugparam.haltnosuccess)
                        MESHIT_LOG_WARNING("no rule found");
                }
            }

            for (size_t i = 0; i < locelements.size() && found; i++) {
                const Element2d& el = locelements[i];
                for (size_t j = 1; j <= 3; j++) {
                    if (el.PNum(j) <= oldnp && pindex[el.PNum(j) - 1] == -1) {
                        found = false;
                        MESHIT_LOG_ERROR("meshing2, index missing");
                    }
                }
            }

            if (found) {
                locpoints.resize(plainpoints.size());

                for (size_t i = oldnp; i < plainpoints.size(); i++) {
                    int err = TransformFromPlain(plainpoints[i], locpoints[i], h);

                    if (err) {
                        found = false;

                        if (debugflag || debugparam.haltnosuccess)
                            MESHIT_LOG_ERROR("meshing2, Backtransformation failed");

                        break;
                    }
                }
            }

            if (found) {
                double violateminh = 3 + 0.1 * qualclass * qualclass;
                double minh = 1e8;
                double newedgemaxh = 0;
                for (size_t i = oldnl; i < loclines.size(); i++) {
                    double eh = Dist(locpoints[loclines[i].I1() - 1],
                                     locpoints[loclines[i].I2() - 1]);

                    // Markus (brute force method to avoid bad elements on geometries like \_/ )
                    //if(eh > 4.*mesh.GetH(locpoints.Get(loclines[i].I1()))) found = 0;
                    //if(eh > 4.*mesh.GetH(locpoints.Get(loclines[i].I2()))) found = 0;
                    // Markus end

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

            if (found && mp.checkoverlap) {

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

                        if (fabs(lam3) < 0.1 * hshould &&
                            lam1 > 0 && lam2 > 0 && (lam1 + lam2) < 1) {

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

                for (int i = oldnl; i < loclines.size(); i++) {
                    int nlgpi1 = loclines[i].I1();
                    int nlgpi2 = loclines[i].I2();
                    if (nlgpi1 <= pindex.size() && nlgpi2 <= pindex.size()) {
                        nlgpi1 = adfront->GetGlobalIndex(pindex[nlgpi1 - 1]);
                        nlgpi2 = adfront->GetGlobalIndex(pindex[nlgpi2 - 1]);

                        int exval = adfront->ExistsLine(nlgpi1, nlgpi2);
                        if (exval) {
                            std::cout << "ERROR: new line exits, val = " << exval << std::endl;
                            std::cerr << "ERROR: new line exits, val = " << exval << std::endl;
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
                        return MESHING2_GIVEUP;
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
                    for (int i = 0; i < pindex.size(); i++) {
                        std::cerr << adfront->GetGlobalIndex(pindex[i]) << std::endl;
                    }
                    std::cerr << "old number of lines = " << oldnl << std::endl;
                    for (size_t i = 0; i < loclines.size(); i++) {
                        std::cerr << "line ";
                        for (size_t j = 1; j <= 2; j++) {
                            int hi = 0;
                            if (loclines[i].I(j) >= 1 &&
                                loclines[i].I(j) <= pindex.size())
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
            } else {
                adfront->IncrementClass(lindex[0]);

                if (debugparam.haltnosuccess || debugflag) {
                    std::cerr << "Problem with seg " << gpi1 << " - " << gpi2
                    << ", class = " << qualclass << std::endl;

                    for (size_t i = 0; i < loclines.size(); i++) {
                        std::cerr << "line ";
                        for (size_t j = 1; j <= 2; j++) {
                            int hi = 0;
                            if (loclines[i].I(j) >= 1 &&
                                loclines[i].I(j) <= pindex.size())
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

        if (!adfront->Empty())
            return MESHING2_GIVEUP;

        return MESHING2_OK;
    }
}
