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

#include "mesh_generator.hpp"

#include "../gprim/geomtest3d.hpp"

namespace meshit {

MeshGenerator::MeshGenerator(Mesh& amesh, const Box2d& boundingbox)
    : mesh_{amesh}, boundingbox_{boundingbox}, adfront{new AdFront2(boundingbox)}, max_area{-1.0}
{
    LoadRules(NULL);
    // LoadRules ("rules/triangle.rls");
}

MeshGenerator::~MeshGenerator()
{
    delete adfront;
    for (size_t i = 0; i < rules.size(); i++) {
        delete rules[i];
    }
}

void MeshGenerator::Reset()
{
    if (adfront) delete adfront;
    adfront = new AdFront2(boundingbox_);
    max_area = -1.0;
}

void MeshGenerator::AddPoint(const Point2d& p, PointIndex globind)
{
    adfront->AddPoint(p, globind);
}

void MeshGenerator::AddBoundaryElement(PointIndex i1, PointIndex i2)
{
    adfront->AddLine(i1, i2);
}

void MeshGenerator::SetMaxArea(double amaxarea)
{
    max_area = amaxarea;
}

void MeshGenerator::DefineTransformation(const Point2d& p1, const Point2d& p2)
{
    glob_p1 = p1;
    ex.X() = p2.X() - p1.X();
    ex.Y() = p2.Y() - p1.Y();
    ex /= ex.Length();
    ey.X() = -ex.Y();
    ey.Y() = ex.X();
}

void MeshGenerator::TransformToPlain(const Point2d& locpoint, Point2d& plainpoint, double h)
{
    Vec2d p1p(glob_p1, locpoint);

    plainpoint.X() = (p1p * ex) / h;
    plainpoint.Y() = (p1p * ey) / h;
}

void MeshGenerator::TransformFromPlain(Point2d& plainpoint, Point2d& locpoint, double h)
{
    locpoint.X() = glob_p1.X() + h * (plainpoint.X() * ex.X() + plainpoint.Y() * ey.X());
    locpoint.Y() = glob_p1.Y() + h * (plainpoint.X() * ex.Y() + plainpoint.Y() * ey.Y());
}

bool MeshGenerator::GenerateMesh(const MeshingParameters& mp, double gh, DomainIndex facenr)
{
    std::vector<PointIndex> pindex;
    std::vector<FrontLineIndex> lindex;
    std::vector<int> delpoints;
    std::vector<uint32_t> dellines;
    std::vector<Element2d> locelements;
    std::vector<Point2d> locpoints;
    std::vector<bool> legalpoints;
    std::vector<Point2d> plainpoints;
    std::vector<IndexPair> loclines;
    uint32_t trials = 0;
    uint32_t plot_next_trial = 999;
    uint32_t qualclass;

    std::vector<int> elements_on_node(mesh_.GetNbPoints(), 0);
    std::vector<bool> illegal_point(mesh_.GetNbPoints(), false);
    double meshed_area = 0.0;

    if (max_area > 0) meshed_area = mesh_.SurfaceArea();

    foundmap.resize(rules.size(), 0);
    canuse.resize(rules.size(), 0);
    ruleused.resize(rules.size(), 0);

    adfront->SetStartFront();


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
        if (++trials > plot_next_trial) {
            MESHIT_LOG_DEBUG_CONT("Trial #" << trials << " : ");
            MESHIT_LOG_DEBUG(adfront->GetNFL() << " faces, " << mesh_.GetNbElements() << " elements.");
            for (size_t i = 0; i < canuse.size(); i++) {
                MESHIT_LOG_DEBUG("  map:" << std::setw(4) << foundmap[i] <<
                                 " / can:" << std::setw(4) << canuse[i] <<
                                 " / use:" << std::setw(4) << ruleused[i] <<
                                 " rule " << rules[i]->Name());
            }
            plot_next_trial += 1000;
        }

        Point2d p1, p2;
        int baselineindex = adfront->SelectBaseLine(p1, p2, qualclass);

        Point2d pmid = Center(p1, p2);
        double pdist = Dist(p1, p2);
        double hshould = std::min(gh, mesh_.GetH(pmid));
        double hinner = (3 + qualclass) * std::max(pdist, hshould);

        mesh_.RestrictLocalH(pmid, hshould);

        adfront->GetLocals(baselineindex, locpoints, loclines, pindex, lindex, 2 * hinner);

        if (qualclass > mp.giveup_tol2d) {
            MESHIT_LOG_WARNING("give up with qualclass "
                               << qualclass << " : number of frontlines = " << adfront->GetNFL());
            break;
        }

        PointIndex gpi1 = adfront->GetGlobalIndex(pindex[loclines[0].I1() - 1]);
        PointIndex gpi2 = adfront->GetGlobalIndex(pindex[loclines[0].I2() - 1]);

        // problem recognition !
        bool found = gpi1 >= static_cast<PointIndex>(illegal_point.size()) ||
                     gpi2 >= static_cast<PointIndex>(illegal_point.size()) ||
                     (!illegal_point[gpi1] && !illegal_point[gpi2]);

        size_t oldnl = 0;
        size_t oldnp = 0;
        int rulenr = -1;

        if (found) {
            DefineTransformation(p1, p2);

            plainpoints.resize(locpoints.size());

            for (size_t i = 0; i < locpoints.size(); i++) {
                TransformToPlain(locpoints[i], plainpoints[i], hshould);
            }

            legalpoints.assign(plainpoints.size(), true);
            double avy = 0;
            for (size_t i = 0; i < plainpoints.size(); i++) {
                avy += plainpoints[i].Y();
            }
            avy *= 1. / plainpoints.size();

            for (size_t i = 0; i < plainpoints.size(); i++) {
                if (pindex[i] == CONST<PointIndex>::undefined) {
                    legalpoints[i] = false;
                }
                if (plainpoints[i].Y() < -1e-10 * avy)  // changed
                {
                    legalpoints[i] = false;
                }
            }

            size_t maxlegalpoint = locpoints.size();
            size_t maxlegalline = loclines.size();

            oldnl = loclines.size();
            oldnp = plainpoints.size();

            rulenr = ApplyRules(plainpoints,
                                legalpoints, maxlegalpoint,
                                loclines, maxlegalline,
                                locelements, dellines, qualclass, mp);
            if (!rulenr) {
                found = false;
            }
        }

        if (found) {
            for (size_t i = 0; i < locelements.size() && found; i++) {
                const Element2d& el = locelements[i];
                for (size_t j = 0; j < 3; j++) {
                    if (el.PointID(j) <= static_cast<PointIndex>(oldnp) &&
                        pindex[el.PointID(j) - 1] == CONST<PointIndex>::undefined) {
                        found = false;
                        MESHIT_LOG_ERROR("meshing2, index missing");
                    }
                }
            }
        }
        if (found) {
            locpoints.resize(plainpoints.size());

            for (size_t i = oldnp; i < plainpoints.size(); i++) {
                TransformFromPlain(plainpoints[i], locpoints[i], hshould);
            }

            double violateminh = 3 + 0.1 * qualclass * qualclass;
            double minh = 1e8;
            double newedgemaxh = 0;
            for (size_t i = oldnl; i < loclines.size(); i++) {
                double eh = Dist(locpoints[loclines[i].I1() - 1],
                                 locpoints[loclines[i].I2() - 1]);
                if (eh > newedgemaxh) newedgemaxh = eh;
            }

            for (size_t i = 0; i < locelements.size(); i++) {
                Point2d pmin = locpoints[locelements[i].PointID(0) - 1];
                Point2d pmax = pmin;
                for (size_t j = 1; j < 3; j++) {
                    const Point2d& hp = locpoints[locelements[i].PointID(j) - 1];
                    pmin.SetToMin(hp);
                    pmax.SetToMax(hp);
                }
                minh = std::min(minh, mesh_.GetMinH(pmin, pmax));
            }

            for (size_t i = 0; i < locelements.size(); i++) {
                found &= Dist2(locpoints[locelements[i].PointID(0) - 1], pmid) <= hinner * hinner;
                found &= Dist2(locpoints[locelements[i].PointID(1) - 1], pmid) <= hinner * hinner;
                found &= Dist2(locpoints[locelements[i].PointID(2) - 1], pmid) <= hinner * hinner;
            }

            static double maxviolate = 0;
            if (newedgemaxh / minh > maxviolate) {
                maxviolate = newedgemaxh / minh;
            }

            if (newedgemaxh > violateminh * minh) {
                found = false;
                loclines.resize(oldnl);
                locpoints.resize(oldnp);
            }
        }

        if (found) {
            // check, whether new front line already exists
            for (size_t i = oldnl; i < loclines.size(); i++) {
                PointIndex nlgpi1 = loclines[i].I1();
                PointIndex nlgpi2 = loclines[i].I2();
                if (nlgpi1 <= static_cast<PointIndex>(pindex.size()) &&
                    nlgpi2 <= static_cast<PointIndex>(pindex.size())) {
                    nlgpi1 = adfront->GetGlobalIndex(pindex[nlgpi1 - 1]);
                    nlgpi2 = adfront->GetGlobalIndex(pindex[nlgpi2 - 1]);

                    bool exval = adfront->LineExists(nlgpi1, nlgpi2);
                    if (exval) {
                        MESHIT_LOG_ERROR("new line (" << nlgpi1 << "," << nlgpi2 << ") exits !");
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
                PointIndex globind = mesh_.AddPoint(locpoints[i]);
                pindex[i] = adfront->AddPoint(locpoints[i], globind);
            }

            for (size_t i = oldnl; i < loclines.size(); i++) {
                PointIndex pi1 = pindex[loclines[i].I1() - 1];
                PointIndex pi2 = pindex[loclines[i].I2() - 1];
                adfront->AddLine(pi1, pi2);
            }
            for (size_t i = 0; i < locelements.size(); i++) {
                Element2d mtri;
                mtri = locelements[i];
                mtri.SetFaceID(facenr);

                for (size_t j = 0; j < 3; j++) {
                    mtri.PointID(j) = locelements[i].PointID(j) =
                        adfront->GetGlobalIndex(pindex[locelements[i].PointID(j) - 1]);
                }

                mesh_.AddElement(mtri);

                const MeshPoint& sep1 = mesh_.Point(mtri.PointID(0));
                const MeshPoint& sep2 = mesh_.Point(mtri.PointID(1));
                const MeshPoint& sep3 = mesh_.Point(mtri.PointID(2));

                double trigarea =
                    Cross(Vec2d(Point2d(sep1), Point2d(sep2)), Vec2d(Point2d(sep1), Point2d(sep3))) / 2;

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
                        std::cerr << "illegal point: " << gpi << ": trigs_on_node = " << elements_on_node[gpi]
                        << std::endl;
                    }

                    static int mtonnode = 0;
                    if (elements_on_node[gpi] > mtonnode) mtonnode = elements_on_node[gpi];
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
    for (size_t i = 0; i < ruleused.size(); i++) {
        MESHIT_LOG_DEBUG(std::setw(5) << ruleused[i] << " times used rule " << rules[i]->Name());
    }

    adfront->PrintOpenSegments(std::cout);

    return (!adfront->Empty());
}

}  // namespace meshit
