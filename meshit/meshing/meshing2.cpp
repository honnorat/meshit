#include "meshing2.hpp"

#include "global.hpp"
#include "../gprim/geomtest3d.hpp"

namespace meshit
{
    Meshing2::Meshing2(const Box3d& boundingbox)
    {
        LoadRules(NULL);
        // LoadRules ("rules/triangle.rls");

        adfront = new AdFront2(boundingbox);
        max_area = -1;
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

    void Meshing2::EndMesh()
    {
        for (size_t i = 0; i < ruleused.size(); i++) {
            MESHIT_LOG_DEBUG(std::setw(5) << ruleused[i] << " times used rule " << rules[i]->Name());
        }
    }

    void Meshing2::SetMaxArea(double amaxarea)
    {
        max_area = amaxarea;
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
        int rulenr{-1};
        Point3d p1, p2;

        bool debugflag;

        double h, hshould;

        std::vector<Point3d> locpoints;
        std::vector<int> legalpoints;
        std::vector<Point2d> plainpoints;
        std::vector<INDEX_2> loclines;
        int trials = 0, nfaces = 0;
        int qualclass;

        std::vector<int> elements_on_node(mesh.GetNP(), 0);
        std::vector<bool> illegal_point(mesh.GetNP(), false);
        double meshed_area = 0.0;

        if (max_area > 0)
            meshed_area = mesh.SurfaceArea();

        foundmap.resize(rules.size(), 0);
        canuse.resize(rules.size(), 0);
        ruleused.resize(rules.size(), 0);

        adfront->SetStartFront();

        int plotnexttrial = 999;

        double meshedarea_before = meshed_area;

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

            Point3d pmid = Center(p1, p2);
            hshould = CalcLocalH(pmid, mesh.GetH(pmid));
            if (gh < hshould) hshould = gh;

            mesh.RestrictLocalH(pmid, hshould);

            h = hshould;

            double hinner = (3 + qualclass) * std::max(Dist(p1, p2), hshould);

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
            found = gpi1 >= static_cast<PointIndex>(illegal_point.size()) ||
                    gpi2 >= static_cast<PointIndex>(illegal_point.size()) ||
                    (!illegal_point[gpi1] && !illegal_point[gpi2]);

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
                    for (size_t j = 0; j < 3; j++) {
                        if (el.PointID(j) <= static_cast<PointIndex>(oldnp) && pindex[el.PointID(j) - 1] == -1) {
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
                    Point3d pmin = locpoints[locelements[i].PointID(0) - 1];
                    Point3d pmax = pmin;
                    for (size_t j = 1; j < 3; j++) {
                        const Point3d& hp = locpoints[locelements[i].PointID(j) - 1];
                        pmin.SetToMin(hp);
                        pmax.SetToMax(hp);
                    }
                    double eh = mesh.GetMinH(pmin, pmax);
                    if (eh < minh)
                        minh = eh;
                }

                for (size_t i = 0; i < locelements.size(); i++) {
                    for (size_t j = 0; j < 3; j++) {
                        if (Dist2(locpoints[locelements[i].PointID(j) - 1], pmid) > hinner * hinner)
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

                    for (size_t j = 0; j < 3; j++) {
                        mtri.PointID(j) = locelements[i].PointID(j) =
                            adfront->GetGlobalIndex(pindex[locelements[i].PointID(j) - 1]);
                    }

                    mesh.AddSurfaceElement(mtri);

                    const Point3d& sep1 = mesh.Point(mtri.PointID(0));
                    const Point3d& sep2 = mesh.Point(mtri.PointID(1));
                    const Point3d& sep3 = mesh.Point(mtri.PointID(2));

                    double trigarea = Cross(Vec3d(sep1, sep2),
                                            Vec3d(sep1, sep3)).Length() / 2;

                    meshed_area += trigarea;

                    if (max_area > 0 && meshed_area - meshedarea_before > max_area) {
                        std::cerr << "meshed area = " << meshed_area - meshedarea_before << std::endl
                        << "maximal area = " << max_area << std::endl
                        << "GIVING UP" << std::endl;
                        return false;
                    }

                    for (size_t j = 0; j < 3; j++) {

                        PointIndex gpi = locelements[i].PointID(j);
                        size_t oldts = elements_on_node.size();
                        if (gpi >= static_cast<PointIndex>(oldts)) {
                            elements_on_node.resize(gpi + 1, 0);
                            illegal_point.resize(gpi + 1, false);
                        }
                        elements_on_node[gpi]++;

                        if (elements_on_node[gpi] > 20) {
                            illegal_point[gpi] = true;
                            std::cerr << "illegal point: " << gpi << ": trigs_on_node = " << elements_on_node[gpi] <<
                            std::endl;
                        }

                        static int mtonnode = 0;
                        if (elements_on_node[gpi] > mtonnode)
                            mtonnode = elements_on_node[gpi];
                    }
                }
                for (size_t i = 0; i < dellines.size(); i++) {
                    adfront->DeleteLine(lindex[dellines[i] - 1]);
                }
            } else {
                adfront->IncrementClass(lindex[0]);
            }
        }

        MESHIT_LOG_DEBUG("Surface meshing done");

        adfront->PrintOpenSegments(std::cout);

        EndMesh();

        return (!adfront->Empty());
    }

}  // namespace meshit
