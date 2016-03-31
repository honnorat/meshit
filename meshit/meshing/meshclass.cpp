#include <stdexcept>
#include <cassert>

#include "meshtool.hpp"
#include "../gprim/geomtest3d.hpp"
#include "../geom2d/geometry2d.hpp"
#include "meshing2.hpp"
#include "../interface/writeuser.hpp"

namespace meshit
{
    Mesh::Mesh()
        : surfarea(*this)
    {
        segment_ht = nullptr;
        lochfunc = nullptr;
        hglob_ = 1e10;
        hmin_ = 0;
        numvertices = 0;
    }

    Mesh::~Mesh()
    {
        delete lochfunc;
        delete segment_ht;

        for (size_t i = 0; i < materials.size(); i++) {
            delete[] materials[i];
        }
    }

    Mesh& Mesh::operator=(const Mesh& mesh2)
    {
        points = mesh2.points;
        segments = mesh2.segments;
        surf_elements = mesh2.surf_elements;
        lockedpoints = mesh2.lockedpoints;
        facedecoding = mesh2.facedecoding;

        return *this;
    }

    void Mesh::BuildFromSplineGeometry(SplineGeometry& geometry, MeshingParameters& mp)
    {
        MESHIT_LOG_DEBUG("Generate Mesh from spline geometry");

        // Take grading from MeshingParameters in priority, or from geometry if not user-defined.
        double grading = mp.grading;
        if (grading < 0) {
            grading = geometry.GetGrading();
        }
        if (grading < 0.01) {
            MESHIT_LOG_WARNING("grading is too small: " << grading << ". We reset it to 0.05");
            grading = 0.01;
        }
        geometry.SetGrading(grading);

        double h = mp.maxh;
        Box<2> bbox = geometry.GetBoundingBox();
        if (bbox.Diam() < h) {
            h = bbox.Diam();
            mp.maxh = h;
        }
        Point3d pmin(bbox.PMin()[0], bbox.PMin()[1], -bbox.Diam());
        Point3d pmax(bbox.PMax()[0], bbox.PMax()[1], bbox.Diam());

        SetLocalH(pmin, pmax, grading);
        hmin_ = mp.minh;
        hglob_ = h;

        geometry.PartitionBoundary(mp, h, *this);

        size_t maxdomnr = 0;
        for (size_t si = 0; si < GetNSeg(); si++) {
            if (segments[si].domin > maxdomnr) maxdomnr = segments[si].domin;
            if (segments[si].domout > maxdomnr) maxdomnr = segments[si].domout;
        }

        facedecoding.resize(0);
        for (size_t i = 1; i <= maxdomnr; i++) {
            facedecoding.push_back(FaceDescriptor(i));
        }

        CalcLocalH();

        size_t bnp = GetNP();  // boundary points

        for (size_t domnr = 1; domnr <= maxdomnr; domnr++) {
            if (geometry.GetDomainMaxh(domnr) > 0) {
                h = geometry.GetDomainMaxh(domnr);
            }
            MESHIT_LOG_DEBUG("Meshing domain " << domnr << " / " << maxdomnr);

            size_t oldnf = GetNSE();

            Meshing2 meshing(mp, Box<3>(pmin, pmax));

            std::vector<int> compress(bnp, -1);
            int cnt = 0;
            for (size_t pi = 0; pi < bnp; pi++) {
                meshing.AddPoint(points[pi], pi);
                cnt++;
                compress[pi] = cnt;
            }
            for (size_t si = 0; si < GetNSeg(); si++) {
                if (segments[si].domin == domnr) {
                    meshing.AddBoundaryElement(compress[segments[si][0]],
                                               compress[segments[si][1]]);
                }
                if (segments[si].domout == domnr) {
                    meshing.AddBoundaryElement(compress[segments[si][1]],
                                               compress[segments[si][0]]);
                }
            }

            mp.check_overlap = 0;
            meshing.GenerateMesh(*this, mp, h, domnr);

            for (size_t sei = oldnf; sei < GetNSE(); sei++) {
                surf_elements[sei].SetIndex(domnr);
            }

            // astrid
            char* material;
            geometry.GetMaterial(domnr, material);
            if (material) {
                SetMaterial(domnr, material);
            }
        }

        int hsteps = mp.optsteps2d;

        mp.optimize2d = "smm";
        mp.optsteps2d = hsteps / 2;
        Optimize2d(*this, mp);

        mp.optimize2d = "Smm";
        mp.optsteps2d = (hsteps + 1) / 2;
        Optimize2d(*this, mp);

        mp.optsteps2d = hsteps;

        Compress();
    }

    size_t Mesh::AddPoint(const Point3d& p, POINTTYPE type)
    {
        size_t pi = points.size();
        points.push_back(MeshPoint(p, type));

        return pi;
    }

    void Mesh::AddSegment(const Segment& s)
    {
        PointIndex maxn = std::max(s[0], s[1]);

        if (maxn < static_cast<PointIndex>(points.size())) {
            if (points[s[0]].Type() > EDGEPOINT) points[s[0]].SetType(EDGEPOINT);
            if (points[s[1]].Type() > EDGEPOINT) points[s[1]].SetType(EDGEPOINT);
        }
        segments.push_back(s);
    }

    void Mesh::AddSurfaceElement(const Element2d& el)
    {
        PointIndex maxn = el[0];
        for (size_t i = 1; i < 3; i++) {
            if (el[i] > maxn) maxn = el[i];
        }
        maxn += 1;

        if (static_cast<size_t>(maxn) <= points.size()) {
            for (size_t i = 0; i < 3; i++) {
                if (points[el[i]].Type() > SURFACEPOINT)
                    points[el[i]].SetType(SURFACEPOINT);
            }
        }

        size_t si = surf_elements.size();
        surf_elements.push_back(el);

        if (el.index > facedecoding.size()) {
            MESHIT_LOG_ERROR("has no facedecoding: fd.size = " << facedecoding.size() << ", ind = " << el.index);
        }

        Element2d& bref = surf_elements.back();
        bref.next = facedecoding[bref.index - 1].firstelement;
        facedecoding[bref.index - 1].firstelement = si;

        if (surfarea.Valid()) {
            surfarea.Add(el);
        }
    }

    void Mesh::Export(const std::string& filetype, const std::string& filename) const
    {
        WriteUserFormat(filetype, const_cast<const Mesh&> (*this), filename);
    }

    void Mesh::Export(const std::string& filetype, std::ostream& os) const
    {
        WriteUserFormat(filetype, const_cast<const Mesh&> (*this), os);
    }

    void Mesh::Export(const std::string& filename) const
    {
        WriteGmsh2Format(const_cast<const Mesh&> (*this), filename);
    }

    void Mesh::Export(std::ostream& os) const
    {
        WriteGmsh2Format(const_cast<const Mesh&> (*this), os);
    }

    void Mesh::Save(const std::string& filename) const
    {
        std::ofstream outfile(filename.c_str());
        Save(outfile);
    }

    void Mesh::Save(std::ostream& outfile) const
    {
        double scale = 1;    // globflags.GetNumFlag ("scale", 1);

        outfile.setf(std::ios::fixed, std::ios::floatfield);
        outfile.setf(std::ios::showpoint);
        outfile.precision(15);

        outfile << "mesh2d" << std::endl;
        outfile << "# surfnr    bcnr      np      p1      p2      p3" << std::endl;
        outfile << "surface_elements" << std::endl << GetNSE() << std::endl;

        for (size_t sei = 0; sei < GetNSE(); sei++) {
            if (surf_elements[sei].GetIndex()) {
                outfile << std::setw(8) << GetFaceDescriptor(surf_elements[sei].GetIndex()).SurfNr() + 1;
                outfile << std::setw(8) << GetFaceDescriptor(surf_elements[sei].GetIndex()).BCProperty();
            } else {
                outfile << "       0       0";
            }

            const Element2d& sel = surf_elements[sei];
            outfile << std::setw(8) << 3;
            outfile << std::setw(8) << sel[0];
            outfile << std::setw(8) << sel[1];
            outfile << std::setw(8) << sel[2] << std::endl;
        }

        outfile << "\n\n";
        outfile << "# surfid      p1      p2  dom_in dom_out";
        outfile << "   ednr1   ednr2                 dist1                 dist2\n";
        outfile << "edge_segments" << std::endl << GetNSeg() << std::endl;

        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            outfile << std::setw(8) << seg.si;
            outfile << std::setw(8) << seg[0];
            outfile << std::setw(8) << seg[1];
            outfile << std::setw(8) << seg.domin;
            outfile << std::setw(8) << seg.domout;
            outfile << std::setw(8) << seg.edgenr;
            outfile << std::setw(8) << seg.epgeominfo[1].edgenr;  // geometry dependent
            outfile << std::setw(22) << seg.epgeominfo[0].dist;  // splineparameter (2D)
            outfile << std::setw(22) << seg.epgeominfo[1].dist;
            outfile << std::endl;
        }

        outfile << "\n\n";
        outfile << "#                    X                     Y                     Z\n";
        outfile << "points" << std::endl << GetNP() << std::endl;
        for (size_t pi = 0; pi < GetNP(); pi++) {
            outfile << std::setw(22) << points[pi].X() / scale;
            outfile << std::setw(22) << points[pi].Y() / scale;
            outfile << std::setw(22) << points[pi].Z() / scale << std::endl;
        }

        int cntmat = 0;
        for (size_t i = 0; i < materials.size(); i++) {
            if (materials[i] && strlen(materials[i])) {
                cntmat++;
            }
        }
        if (cntmat) {
            outfile << "materials" << std::endl << cntmat << std::endl;
            for (size_t i = 0; i < materials.size(); i++) {
                if (materials[i] && strlen(materials[i]))
                    outfile << i + 1 << " " << materials[i] << std::endl;
            }
        }
    }

    void Mesh::Load(const std::string& filename)
    {
        std::ifstream infile(filename.c_str());
        if (!infile.good())
            throw std::runtime_error("mesh file not found");

        Load(infile);
    }

    void Mesh::Load(std::istream& infile)
    {
        char str[100];
        int n;

        double scale = 1;    // globflags.GetNumFlag ("scale", 1);

        facedecoding.resize(0);

        bool endmesh = false;

        while (infile.good() && !endmesh) {
            infile >> str;
            if (strcmp(str, "surface_elements") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " surface elements");

                for (int i = 1; i <= n; i++) {
                    int surfnr, bcp, nep, faceind = 0;

                    infile >> surfnr >> bcp;

                    if (surfnr < 1) {
                        MESHIT_LOG_ERROR("Invalid surface index number: " << surfnr);
                        throw std::runtime_error("Aborting.");
                    }
                    surfnr--;

                    for (size_t j = 0; j < facedecoding.size(); j++) {
                        if (GetFaceDescriptor(j + 1).SurfNr() == static_cast<size_t>(surfnr) &&
                            GetFaceDescriptor(j + 1).BCProperty() == static_cast<size_t>(bcp)) {
                            faceind = j + 1;
                        }
                    }
                    if (!faceind) {
                        facedecoding.push_back(FaceDescriptor(surfnr));
                        faceind = facedecoding.size();
                        GetFaceDescriptor(faceind).SetBCProperty(bcp);
                    }

                    infile >> nep;
                    if (!nep) nep = 3;

                    if (nep != 3) {
                        MESHIT_LOG_FATAL("Mesh::Load: undefined element type. nep = " << nep << ". Aborting");
                        exit(1);
                    }

                    Element2d tri;
                    tri.SetIndex(faceind);

                    for (int j = 1; j <= nep; j++) {
                        infile >> tri.PNum(j);
                    }
                    AddSurfaceElement(tri);
                }
            }
            if (strcmp(str, "edge_segments") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    Segment seg;
                    int hi;
                    infile >> seg.si >> hi >> seg[0] >> seg[1] >> hi >> hi
                    >> seg.surfnr1 >> seg.surfnr2
                    >> seg.edgenr
                    >> seg.epgeominfo[1].edgenr
                    >> seg.epgeominfo[0].dist
                    >> seg.epgeominfo[1].dist;

                    seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;

                    seg.domin = seg.surfnr1;
                    seg.domout = seg.surfnr2;

                    seg.surfnr1--;
                    seg.surfnr2--;

                    AddSegment(seg);
                }
            }
            if (strcmp(str, "points") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " points");
                for (int i = 1; i <= n; i++) {
                    Point3d p;
                    infile >> p.X() >> p.Y() >> p.Z();
                    p.X() *= scale;
                    p.Y() *= scale;
                    p.Z() *= scale;
                    AddPoint(p);
                }
            }
            if (strcmp(str, "materials") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " materials");
                for (int i = 1; i <= n; i++) {
                    int nr;
                    std::string mat;
                    infile >> nr >> mat;
                    SetMaterial(nr, mat.c_str());
                }
            }
            if (strcmp(str, "endmesh") == 0) {
                endmesh = true;
            }
//            strcpy(str, "");
        }

        CalcSurfacesOfNode();
    }

    void Mesh::CalcSurfacesOfNode()
    {
        surfaces_on_node.SetSize(GetNP());

        delete segment_ht;
        segment_ht = new INDEX_2_CLOSED_HASHTABLE<int>(3 * GetNSeg() + 1);

        for (size_t sei = 0; sei < GetNSE(); sei++) {
            const Element2d& sel = surf_elements[sei];
            if (sel.IsDeleted()) continue;

            int si = sel.GetIndex();

            for (size_t j = 0; j < 3; j++) {
                PointIndex pi = sel[j];
                std::vector<int> surf_idx = surfaces_on_node[pi];
                bool found = 0;
                for (size_t k = 0; k < surf_idx.size(); k++) {
                    if (surf_idx[k] == si) {
                        found = 1;
                        break;
                    }
                }
                if (!found)
                    surfaces_on_node.Add(pi, si);
            }
        }

        for (size_t i = 0; i < segments.size(); i++) {
            const Segment& seg = segments[i];
            MeshPoint& mp1 = points[seg[0]];
            MeshPoint& mp2 = points[seg[1]];
            if (mp1.Type() == INNERPOINT || mp1.Type() == SURFACEPOINT) mp1.SetType(EDGEPOINT);
            if (mp2.Type() == INNERPOINT || mp2.Type() == SURFACEPOINT) mp2.SetType(EDGEPOINT);
        }

        for (size_t i = 0; i < lockedpoints.size(); i++) {
            points[lockedpoints[i]].SetType(FIXEDPOINT);
        }

        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = segments[i];
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();

            segment_ht->Set(i2, i);
        }
    }

    void Mesh::SetLocalH(const Point3d& pmin, const Point3d& pmax, double grading)
    {
        Point3d c = Center(pmin, pmax);
        double d = 0.5 * std::max(pmax.X() - pmin.X(), std::max(pmax.Y() - pmin.Y(), pmax.Z() - pmin.Z()));
        Point3d pmin2 = c - Vec3d(d, d, d);
        Point3d pmax2 = c + Vec3d(d, d, d);

        delete lochfunc;
        lochfunc = new LocalH(pmin2, pmax2, grading);
    }

    void Mesh::RestrictLocalH(const Point3d& p, double hloc)
    {
        if (hloc < hmin_)
            hloc = hmin_;

        assert(lochfunc);
        lochfunc->SetH(p, hloc);
    }

    void Mesh::RestrictLocalHLine(const Point3d& p1, const Point3d& p2, double hloc)
    {
        if (hloc < hmin_) {
            hloc = hmin_;
        }
        int steps = static_cast<int>(Dist(p1, p2) / hloc) + 2;
        Vec3d v(p1, p2);

        for (int i = 0; i <= steps; i++) {
            Point3d p = p1 + (static_cast<double>(i) / static_cast<double>(steps) * v);
            RestrictLocalH(p, hloc);
        }
    }

    double Mesh::GetH(const Point3d& p) const
    {
        double hmin = hglob_;
        if (lochfunc) {
            double hl = lochfunc->GetH(p);
            if (hl < hglob_)
                hmin = hl;
        }
        return hmin;
    }

    double Mesh::GetMinH(const Point3d& pmin, const Point3d& pmax)
    {
        double hmin = hglob_;
        if (lochfunc) {
            double hl = lochfunc->GetMinH(pmin, pmax);
            if (hl < hmin)
                hmin = hl;
        }
        return hmin;
    }

    double Mesh::AverageH(size_t surfnr) const
    {
        double maxh = 0, minh = 1e10;
        double hsum = 0;
        size_t n = 0;
        for (size_t i = 0; i < GetNSE(); i++) {
            const Element2d& el = SurfaceElement(i);
            if (surfnr == 0 || el.GetIndex() == surfnr) {
                for (size_t j = 1; j <= 3; j++) {
                    double hi = Dist(Point(el.PNumMod(j)),
                                     Point(el.PNumMod(j + 1)));

                    hsum += hi;

                    if (hi > maxh) maxh = hi;
                    if (hi < minh) minh = hi;
                    n++;
                }
            }
        }

        MESHIT_LOG_DEBUG("minh = " << minh << " avh = " << (hsum / n) << " maxh = " << maxh);
        return (hsum / n);
    }

    void Mesh::CalcLocalH()
    {
        assert(lochfunc);

        MESHIT_LOG_DEBUG("CalcLocalH: " << GetNP() << " points, "
                         << GetNSE() << " surface elements.");

        for (size_t i = 0; i < GetNSE(); i++) {
            const Element2d& el = surf_elements[i];
            double hel = -1;
            for (size_t j = 1; j <= 3; j++) {
                const Point3d& p1 = points[el.PNumMod(j)];
                const Point3d& p2 = points[el.PNumMod(j + 1)];
                double hedge = Dist(p1, p2);
                if (hedge > hel) hel = hedge;
            }
            if (hel > 0) {
                const Point3d& p1 = points[el.PNum(1)];
                const Point3d& p2 = points[el.PNum(2)];
                const Point3d& p3 = points[el.PNum(3)];
                lochfunc->SetH(Center(p1, p2, p3), hel);
            }
        }

        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = segments[i];
            const Point3d& p1 = points[seg[0]];
            const Point3d& p2 = points[seg[1]];
            lochfunc->SetH(Center(p1, p2), Dist(p1, p2));
        }
    }

    void Mesh::RestrictLocalH(resthtype rht, size_t nr, double loc_h)
    {
        switch (rht) {
            case RESTRICTH_FACE: {
                for (size_t i = 0; i < GetNSE(); i++) {
                    const Element2d& sel = SurfaceElement(i);
                    if (sel.GetIndex() == nr)
                        RestrictLocalH(RESTRICTH_SURFACEELEMENT, i + 1, loc_h);
                }
                break;
            }
            case RESTRICTH_EDGE: {
                for (size_t i = 0; i < GetNSeg(); i++) {
                    const Segment& seg = LineSegment(i);
                    if (seg.edgenr == static_cast<int>(nr))
                        RestrictLocalH(RESTRICTH_SEGMENT, i + 1, loc_h);
                }
                break;
            }
            case RESTRICTH_POINT: {
                RestrictLocalH(points[nr - 1], loc_h);
                break;
            }

            case RESTRICTH_SURFACEELEMENT: {
                const Element2d& sel = SurfaceElement(nr - 1);
                Point3d p = Center(Point(sel.PNum(1)),
                                   Point(sel.PNum(2)),
                                   Point(sel.PNum(3)));
                RestrictLocalH(p, loc_h);
                break;
            }
            case RESTRICTH_SEGMENT: {
                const Segment& seg = LineSegment(nr + 1);
                RestrictLocalHLine(Point(seg[0]), Point(seg[1]), loc_h);
                break;
            }
        }
    }

    void Mesh::GetBox(Point3d& pmin, Point3d& pmax) const
    {
        if (points.size() == 0) {
            pmin = pmax = Point3d(0, 0, 0);
            return;
        }

        pmin = Point3d(1e10, 1e10, 1e10);
        pmax = Point3d(-1e10, -1e10, -1e10);
        for (size_t pi = 0; pi < points.size(); pi++) {
            pmin.SetToMin(points[pi]);
            pmax.SetToMax(points[pi]);
        }
        if (pmin.X() > 0.5e10) {
            pmin = pmax = Point3d(0, 0, 0);
        }
    }

    void Mesh::AddLockedPoint(PointIndex pi)
    {
        lockedpoints.push_back(pi);
    }

    void Mesh::ClearLockedPoints()
    {
        lockedpoints.resize(0);
    }

    void Mesh::Compress()
    {
        std::vector<PointIndex> op2np(GetNP());
        std::vector<MeshPoint> hpoints;
        BitArrayChar pused(GetNP());

        for (size_t i = 0; i < surf_elements.size(); i++) {
            if (surf_elements[i].IsDeleted()) {
                surf_elements.erase(surf_elements.begin() + i);
                i--;
            }
        }
        for (size_t i = 0; i < segments.size(); i++) {
            if (segments[i][0] <= -1) {
                segments.erase(segments.begin() + i);
                i--;
            }
        }
        pused.Clear();

        for (size_t i = 0; i < surf_elements.size(); i++) {
            const Element2d& el = surf_elements[i];
            for (size_t j = 0; j < 3; j++) {
                pused.Set(el[j]);
            }
        }

        for (size_t i = 0; i < segments.size(); i++) {
            const Segment& seg = segments[i];
            pused.Set(seg[0]);
            pused.Set(seg[1]);
        }

        for (size_t i = 0; i < lockedpoints.size(); i++) {
            pused.Set(lockedpoints[i]);
        }

        int npi = -1;

        for (size_t pi = 0; pi < points.size(); pi++) {
            if (pused.Test(pi)) {
                npi++;
                op2np[pi] = npi;
                hpoints.push_back(points[pi]);
            }
            else
                op2np[pi] = -1;
        }

        points.resize(0);
        for (size_t i = 0; i < hpoints.size(); i++) {
            points.push_back(hpoints[i]);
        }

        for (size_t i = 0; i < surf_elements.size(); i++) {
            Element2d& el = SurfaceElement(i);
            for (size_t j = 0; j < 3; j++) {
                el[j] = op2np[el[j]];
            }
        }

        for (size_t i = 0; i < segments.size(); i++) {
            Segment& seg = segments[i];
            seg[0] = op2np[seg[0]];
            seg[1] = op2np[seg[1]];
        }

        for (size_t i = 0; i < lockedpoints.size(); i++) {
            lockedpoints[i] = op2np[lockedpoints[i]];
        }

        for (size_t i = 0; i < facedecoding.size(); i++) {
            facedecoding[i].firstelement = -1;
        }
        for (int i = surf_elements.size() - 1; i >= 0; i--) {
            int ind = surf_elements[i].GetIndex();
            surf_elements[i].next = facedecoding[ind - 1].firstelement;
            facedecoding[ind - 1].firstelement = i;
        }

        CalcSurfacesOfNode();
    }

    int Mesh::CheckOverlappingBoundary()
    {
        Point3d pmin, pmax;
        GetBox(pmin, pmax);
        Box3dTree setree(pmin, pmax);
        std::vector<size_t> inters;

        bool overlap = false;

        for (size_t i = 0; i < GetNSE(); i++) {
            const Element2d& tri = SurfaceElement(i);

            Point3d tpmin(Point(tri[0]));
            Point3d tpmax(tpmin);

            for (size_t k = 1; k < 3; k++) {
                tpmin.SetToMin(Point(tri[k]));
                tpmax.SetToMax(Point(tri[k]));
            }
            Vec3d diag(tpmin, tpmax);

            tpmax = tpmax + 0.1 * diag;
            tpmin = tpmin - 0.1 * diag;

            setree.Insert(tpmin, tpmax, i + 1);
        }

        for (size_t i = 0; i < GetNSE(); i++) {
            const Element2d& tri1 = SurfaceElement(i);

            Point3d tpmin(Point(tri1[0]));
            Point3d tpmax(tpmin);

            for (size_t k = 1; k < 3; k++) {
                tpmin.SetToMin(Point(tri1[k]));
                tpmax.SetToMax(Point(tri1[k]));
            }

            setree.GetIntersecting(tpmin, tpmax, inters);

            for (size_t j = 0; j < inters.size(); j++) {
                const Element2d& tri2 = SurfaceElement(inters[j] - 1);

                const meshit::Point3d* trip1[3], * trip2[3];
                for (size_t k = 1; k <= 3; k++) {
                    trip1[k - 1] = &Point(tri1.PNum(k));
                    trip2[k - 1] = &Point(tri2.PNum(k));
                }

                if (IntersectTriangleTriangle(&trip1[0], &trip2[0])) {
                    overlap = 1;
                    MESHIT_LOG_WARNING("Intersecting elements " << i + 1 << " and " << inters[j]);
                    MESHIT_LOG_DEBUG(" el1 = " << tri1);
                    MESHIT_LOG_DEBUG(" el2 = " << tri2);

                    for (size_t k = 1; k <= 3; k++)
                        MESHIT_LOG_DEBUG_CONT(tri1.PNum(k) << "  ");
                    MESHIT_LOG_DEBUG("");
                    for (size_t k = 1; k <= 3; k++)
                        MESHIT_LOG_DEBUG_CONT(tri2.PNum(k) << "  ");
                    MESHIT_LOG_DEBUG("");

                    for (size_t k = 0; k <= 2; k++)
                        MESHIT_LOG_DEBUG_CONT(*trip1[k] << "   ");
                    MESHIT_LOG_DEBUG("");
                    for (size_t k = 0; k <= 2; k++)
                        MESHIT_LOG_DEBUG_CONT(*trip2[k] << "   ");
                    MESHIT_LOG_DEBUG("");

                    MESHIT_LOG_DEBUG("Face1 = " << GetFaceDescriptor(tri1.GetIndex()));
                    MESHIT_LOG_DEBUG("Face2 = " << GetFaceDescriptor(tri2.GetIndex()));
                }
            }
        }
        return overlap;
    }

    bool Mesh::PointContainedIn2DElement(const Point3d& p,
                                         double lami[3],
                                         const int element) const
    {
        const double eps = 1e-6;

        const Element2d& el = SurfaceElement(element);
        const Point3d& p1 = Point(el.PNum(1));
        const Point3d& p2 = Point(el.PNum(2));
        const Point3d& p3 = Point(el.PNum(3));

        Vec3d col1(p1, p2);
        Vec3d col2(p1, p3);
        Vec3d col3 = Cross(col1, col2);
        Vec3d rhs(p1, p);
        Vec3d sol;

        SolveLinearSystem(col1, col2, col3, rhs, sol);

        if (sol.X() >= -eps && sol.Y() >= -eps && sol.X() + sol.Y() <= 1 + eps) {
            lami[0] = sol.X();
            lami[1] = sol.Y();
            lami[2] = sol.Z();
            return true;
        }
        return false;
    }

    void Mesh::RebuildSurfaceElementLists()
    {
        for (size_t i = 0; i < facedecoding.size(); i++) {
            facedecoding[i].firstelement = -1;
        }
        for (int i = surf_elements.size() - 1; i >= 0; i--) {
            int ind = surf_elements[i].GetIndex();
            surf_elements[i].next = facedecoding[ind - 1].firstelement;
            facedecoding[ind - 1].firstelement = i;
        }
    }

    void Mesh::GetSurfaceElementsOfFace(size_t facenr, std::vector<SurfaceElementIndex>& sei) const
    {
        sei.clear();

        SurfaceElementIndex si = facedecoding[facenr - 1].firstelement;
        while (si != -1) {
            const Element2d& se = SurfaceElement(si);
            if (se.GetIndex() == facenr && se.PNum(1) >= 0 && !se.IsDeleted()) {
                sei.push_back(si);
            }
            si = se.next;
        }
    }

    void Mesh::ComputeNVertices()
    {
        numvertices = 0;
        for (size_t i = 0; i < surf_elements.size(); i++) {
            const Element2d& el = SurfaceElement(i);
            for (size_t j = 1; j <= 3; j++) {
                if (el.PNum(j) > static_cast<PointIndex>(numvertices)) {
                    numvertices = static_cast<size_t>(el.PNum(j));
                }
            }
        }
        numvertices += 1;
        MESHIT_LOG_INFO("Mesh::ComputeNVertices: numvertices = " << numvertices << " numpoints = " << GetNP());
    }

    void Mesh::SetNP(size_t np)
    {
        points.resize(np);
    }

    void Mesh::SetMaterial(size_t domnr, const char* mat)
    {
        if (domnr > materials.size()) {
            size_t olds = materials.size();
            materials.resize(domnr);
            for (size_t i = olds; i < domnr; i++) {
                materials[i] = 0;
            }
        }
        materials[domnr - 1] = new char[strlen(mat) + 1];
        strcpy(materials[domnr - 1], mat);
    }

    void Mesh::PrintMemInfo(std::ostream& ost) const
    {
        ost << "Mesh Mem:" << std::endl;

        ost << GetNP() << " Points, of size "
        << sizeof(Point3d) << " + " << sizeof(POINTTYPE) << " = "
        << GetNP() * (sizeof(Point3d) + sizeof(POINTTYPE)) << std::endl;

        ost << GetNSE() << " Surface elements, of size "
        << sizeof(Element2d) << " = "
        << GetNSE() * sizeof(Element2d) << std::endl;

        ost << "surfs on node:";
        surfaces_on_node.PrintMemInfo(ost);
    }

}  // namespace meshit
