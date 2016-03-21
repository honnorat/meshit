#include <stdexcept>
#include <cassert>

#include "meshtool.hpp"
#include "../gprim/geomtest3d.hpp"
#include "../geom2d/geometry2d.hpp"
#include "meshing2.hpp"
#include "../interface/writeuser.hpp"

namespace meshit {

    Mesh::Mesh()
            : surfarea(*this)
    {
        boundaryedges = NULL;
        surfelementht = NULL;
        segmentht = NULL;
        lochfunc = NULL;
        timestamp = NextTimeStamp();
        hglob = 1e10;
        hmin = 0;
        numvertices = -1;

        topology = new MeshTopology(*this);
        ident = new Identifications(*this);
    }

    Mesh::~Mesh()
    {
        delete lochfunc;
        delete boundaryedges;
        delete surfelementht;
        delete segmentht;
        delete topology;
        delete ident;

        for (size_t i = 0; i < materials.size(); i++) {
            delete[] materials[i];
        }
    }

    Mesh& Mesh::operator=(const Mesh& mesh2)
    {
        points = mesh2.points;
        segments = mesh2.segments;
        surfelements = mesh2.surfelements;
        lockedpoints = mesh2.lockedpoints;
        facedecoding = mesh2.facedecoding;

        return *this;
    }

    void Mesh::BuildFromSpline2D(SplineGeometry2d& geometry, MeshingParameters& mp)
    {
        MESHIT_LOG_DEBUG("Generate Mesh from spline geometry");

        // Take grading from MeshingParameters in priority, or from geometry if not user-defined.
        double grading = mp.grading;
        if (grading < 0) {
            grading = geometry.GetGrading();
        }

        double h = mp.maxh;
        Box<2> bbox = geometry.GetBoundingBox();
        if (bbox.Diam() < h) {
            h = bbox.Diam();
            mp.maxh = h;
        }
        Point3d pmin(bbox.PMin()[0], bbox.PMin()[1], -bbox.Diam());
        Point3d pmax(bbox.PMax()[0], bbox.PMax()[1], bbox.Diam());

        SetLocalH(pmin, pmax, grading);
        hglob = h;

        geometry.PartitionBoundary(mp, h, *this);

        int maxdomnr = 0;
        for (size_t si = 0; si < GetNSeg(); si++) {
            if (segments[si].domin > maxdomnr) maxdomnr = segments[si].domin;
            if (segments[si].domout > maxdomnr) maxdomnr = segments[si].domout;
        }

        facedecoding.resize(0);
        for (int i = 1; i <= maxdomnr; i++) {
            facedecoding.push_back(FaceDescriptor(i, 0, 0, i));
        }

        CalcLocalH();

        size_t bnp = GetNP();  // boundary points

        for (int domnr = 1; domnr <= maxdomnr; domnr++) {
            if (geometry.GetDomainMaxh(domnr) > 0) {
                h = geometry.GetDomainMaxh(domnr);
            }
            MESHIT_LOG_DEBUG("Meshing domain " << domnr << " / " << maxdomnr);

            int oldnf = GetNSE();

            Meshing2 meshing(mp, Box<3>(pmin, pmax));

            std::vector<int> compress(bnp, -1);
            int cnt = 0;
            for (size_t pi = 0; pi < bnp; pi++) {
                if (points[pi].GetLayer() == geometry.GetDomainLayer(domnr)) {
                    meshing.AddPoint(points[pi], pi);
                    cnt++;
                    compress[pi] = cnt;
                }
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
                surfelements[sei].SetIndex(domnr);
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

    PointIndex Mesh::AddPoint(const Point3d& p, int layer, POINTTYPE type)
    {
        timestamp = NextTimeStamp();

        PointIndex pi = points.size();
        points.push_back(MeshPoint(p, layer, type));

        return pi;
    }

    void Mesh::AddSegment(const Segment& s)
    {
        timestamp = NextTimeStamp();

        size_t maxn = std::max(s[0], s[1]);

        if (maxn < points.size()) {
            if (points[s[0]].Type() > EDGEPOINT) points[s[0]].SetType(EDGEPOINT);
            if (points[s[1]].Type() > EDGEPOINT) points[s[1]].SetType(EDGEPOINT);
        }
        segments.push_back(s);
    }

    void Mesh::AddSurfaceElement(const Element2d& el)
    {
        timestamp = NextTimeStamp();

        PointIndex maxn = el[0];
        for (size_t i = 1; i < 3; i++) {
            if (el[i] > maxn) maxn = el[i];
        }
        maxn += 1;

        if (maxn <= points.size()) {
            for (size_t i = 0; i < 3; i++) {
                if (points[el[i]].Type() > SURFACEPOINT)
                    points[el[i]].SetType(SURFACEPOINT);
            }
        }

        size_t si = surfelements.size();
        surfelements.push_back(el);

        if (el.index > facedecoding.size()) {
            MESHIT_LOG_ERROR("has no facedecoding: fd.size = " << facedecoding.size() << ", ind = " << el.index);
        }

        Element2d& bref = surfelements.back();
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
        int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");

        outfile << "mesh3d" << "\n";
        outfile << "dimension\n" << 2 << "\n";
        outfile << "geomtype\n" << 0 << "\n";
        outfile << "\n";
        outfile << "# surfnr    bcnr   domin  domout      np      p1      p2      p3" << "\n";
        outfile << "surfaceelements" << "\n";
        outfile << GetNSE() << "\n";

        for (size_t sei = 0; sei < GetNSE(); sei++) {
            if (surfelements[sei].GetIndex()) {
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).SurfNr() + 1;
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).BCProperty();
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).DomainIn();
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).DomainOut();
            } else {
                outfile << " 0 0 0";
            }

            Element2d sel = surfelements[sei];
            if (invertsurf)
                sel.Invert();

            outfile << " " << 3;
            for (size_t j = 0; j < 3; j++) {
                outfile << " " << sel[j];
            }
            outfile << "\n";
        }

        outfile << "\n" << "\n";
        outfile << "# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2   ";
        outfile << "ednr1   dist1   ednr2   dist2 \n";
        outfile << "edgesegmentsgi2" << "\n";
        outfile << GetNSeg() << "\n";

        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            outfile.width(8);
            outfile << seg.si;  // 2D: bc number, 3D: wievielte Kante
            outfile.width(8);
            outfile << 0;
            outfile.width(8);
            outfile << seg[0];
            outfile.width(8);
            outfile << seg[1];
            outfile << " ";
            outfile.width(8);
            outfile << 1;  // stl dreiecke
            outfile << " ";
            outfile.width(8);
            outfile << 1;  // <<std::endl;  // stl dreieck

            outfile << " ";
            outfile.width(8);
            outfile << seg.domin;
            outfile << " ";
            outfile.width(8);
            outfile << seg.domout;

            outfile << " ";
            outfile.width(8);
            outfile << seg.edgenr;
            outfile << " ";
            outfile.width(12);
            outfile.precision(16);
            outfile << seg.epgeominfo[0].dist;  // splineparameter (2D)
            outfile << " ";
            outfile.width(8);
            outfile.precision(16);
            outfile << seg.epgeominfo[1].edgenr;  // geometry dependent
            outfile << " ";
            outfile.width(12);
            outfile << seg.epgeominfo[1].dist;

            outfile << "\n";
        }

        outfile << "\n" << "\n";
        outfile << "#          X             Y             Z" << "\n";
        outfile << "points" << "\n";
        outfile << GetNP() << "\n";
        outfile.precision(16);
        outfile.setf(std::ios::fixed, std::ios::floatfield);
        outfile.setf(std::ios::showpoint);

        for (PointIndex pi = 0; pi < GetNP(); pi++) {
            outfile.width(22);
            outfile << points[pi].X() / scale << "  ";
            outfile.width(22);
            outfile << points[pi].Y() / scale << "  ";
            outfile.width(22);
            outfile << points[pi].Z() / scale << "\n";
        }

        if (ident->GetMaxNr() > 0) {
            outfile << "identifications\n";
            std::vector<INDEX_2> identpairs;
            int cnt = 0;
            for (int i = 1; i <= ident->GetMaxNr(); i++) {
                ident->GetPairs(i, identpairs);
                cnt += identpairs.size();
            }
            outfile << cnt << "\n";
            for (int i = 1; i <= ident->GetMaxNr(); i++) {
                ident->GetPairs(i, identpairs);
                for (int j = 0; j < identpairs.size(); j++) {
                    outfile.width(8);
                    outfile << identpairs[j].I1();
                    outfile.width(8);
                    outfile << identpairs[j].I2();
                    outfile.width(8);
                    outfile << i << "\n";
                }
            }

            outfile << "identificationtypes\n";
            outfile << ident->GetMaxNr() << "\n";
            for (size_t i = 0; i < ident->GetMaxNr(); i++) {
                int type = ident->GetType(i);
                outfile << " " << type;
            }
            outfile << "\n";
        }

        int cntmat = 0;
        for (int i = 0; i < materials.size(); i++) {
            if (materials[i] && strlen(materials[i]))
                cntmat++;
        }
        if (cntmat) {
            outfile << "materials" << std::endl;
            outfile << cntmat << std::endl;
            for (int i = 0; i < materials.size(); i++) {
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
        int invertsurf = 0;  // globflags.GetDefineFlag ("invertsurfacemesh");

        facedecoding.resize(0);

        bool endmesh = false;

        while (infile.good() && !endmesh) {
            infile >> str;

            if (strcmp(str, "geomtype") == 0) {
                int hi;
                infile >> hi;
            }
            if (strcmp(str, "surfaceelements") == 0 || strcmp(str, "surfaceelementsgi") == 0 ||
                strcmp(str, "surfaceelementsuv") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " surface elements");

                bool geominfo = strcmp(str, "surfaceelementsgi") == 0;
                bool uv = strcmp(str, "surfaceelementsuv") == 0;

                for (int i = 1; i <= n; i++) {
                    int surfnr, bcp, domin, domout, nep, faceind = 0;

                    infile >> surfnr >> bcp >> domin >> domout;
                    surfnr--;

                    bool invert_el = false;

                    for (int j = 1; j <= facedecoding.size(); j++) {
                        if (GetFaceDescriptor(j).SurfNr() == surfnr &&
                            GetFaceDescriptor(j).BCProperty() == bcp &&
                            GetFaceDescriptor(j).DomainIn() == domin &&
                            GetFaceDescriptor(j).DomainOut() == domout) {
                            faceind = j;
                        }
                    }
                    if (!faceind) {
                        facedecoding.push_back(FaceDescriptor(surfnr, domin, domout, 0));
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
                    if (geominfo) {
                        int hi;
                        for (int j = 1; j <= nep; j++) {
                            infile >> hi;
                        }
                    }
                    if (uv) {
                        int hi;
                        for (int j = 1; j <= nep; j++) {
                            infile >> hi >> hi;
                        }
                    }
                    if (invertsurf) tri.Invert();
                    if (invert_el) tri.Invert();

                    AddSurfaceElement(tri);
                }
            }
            if (strcmp(str, "edgesegments") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    Segment seg;
                    int hi;
                    infile >> seg.si >> hi >> seg[0] >> seg[1];
                    AddSegment(seg);
                }
            }
            if (strcmp(str, "edgesegmentsgi") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    Segment seg;
                    int hi;
                    infile >> seg.si >> hi >> seg[0] >> seg[1] >> hi >> hi;
                    AddSegment(seg);
                }
            }
            if (strcmp(str, "edgesegmentsgi2") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " curve elements");
                for (int i = 1; i <= n; i++) {
                    Segment seg;
                    int hi;
                    infile >> seg.si >> hi >> seg[0] >> seg[1] >> hi >> hi
                    >> seg.surfnr1 >> seg.surfnr2
                    >> seg.edgenr
                    >> seg.epgeominfo[0].dist
                    >> seg.epgeominfo[1].edgenr
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
            if (strcmp(str, "identifications") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " identifications");
                for (int i = 1; i <= n; i++) {
                    PointIndex pi1, pi2;
                    int ind;
                    infile >> pi1 >> pi2 >> ind;
                    ident->Add(pi1, pi2, ind);
                }
            }
            if (strcmp(str, "identificationtypes") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " identificationtypes");
                for (size_t i = 0; i < n; i++) {
                    int type;
                    infile >> type;
                    ident->SetType(i, Identifications::ID_TYPE(type));
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
            if (strcmp(str, "endmesh") == 0)
                endmesh = true;

            strcpy(str, "");
        }

        CalcSurfacesOfNode();

        topology->Update();
    }

    void Mesh::BuildBoundaryEdges(void)
    {
        delete boundaryedges;

        boundaryedges = new INDEX_2_CLOSED_HASHTABLE<int>
                (3 * (GetNSE() + GetNOpenElements()) + GetNSeg() + 1);

        for (size_t sei = 0; sei < GetNSE(); sei++) {
            const Element2d& sel = surfelements[sei];
            if (sel.IsDeleted()) continue;

            for (size_t j = 0; j < 3; j++) {
                INDEX_2 i2;
                i2.I1() = sel.PNumMod(j + 1);
                i2.I2() = sel.PNumMod(j + 2);
                i2.Sort();
                boundaryedges->Set(i2, 1);
            }
        }

        for (int i = 0; i < openelements.size(); i++) {
            const Element2d& sel = openelements[i];
            for (size_t j = 0; j < 3; j++) {
                INDEX_2 i2;
                i2.I1() = sel.PNumMod(j + 1);
                i2.I2() = sel.PNumMod(j + 2);
                i2.Sort();
                boundaryedges->Set(i2, 1);

                points[sel[j]].SetType(FIXEDPOINT);
            }
        }

        for (int i = 0; i < GetNSeg(); i++) {
            const Segment& seg = segments[i];
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();
            boundaryedges->Set(i2, 2);
        }
    }

    void Mesh::CalcSurfacesOfNode()
    {
        surfacesonnode.SetSize(GetNP());

        delete boundaryedges;
        boundaryedges = NULL;

        delete surfelementht;
        delete segmentht;

        surfelementht = new INDEX_3_CLOSED_HASHTABLE<int>(3 * GetNSE() + 1);
        segmentht = new INDEX_2_CLOSED_HASHTABLE<int>(3 * GetNSeg() + 1);

        for (size_t sei = 0; sei < GetNSE(); sei++) {
            const Element2d& sel = surfelements[sei];
            if (sel.IsDeleted()) continue;

            int si = sel.GetIndex();

            for (size_t j = 0; j < 3; j++) {
                PointIndex pi = sel[j];
                bool found = 0;
                for (int k = 0; k < surfacesonnode[pi].size(); k++) {
                    if (surfacesonnode[pi][k] == si) {
                        found = 1;
                        break;
                    }
                }
                if (!found)
                    surfacesonnode.Add(pi, si);
            }
        }

        for (size_t sei = 0; sei < GetNSE(); sei++) {
            const Element2d& sel = surfelements[sei];
            if (sel.IsDeleted()) continue;

            INDEX_3 i3;
            i3.I1() = sel.PNum(1);
            i3.I2() = sel.PNum(2);
            i3.I3() = sel.PNum(3);
            i3.Sort();
            surfelementht->Set(i3, sei); // war das wichtig ???    sel.GetIndex());
        }

        for (int i = 0; i < segments.size(); i++) {
            const Segment& seg = segments[i];
            for (int j = 1; j <= 2; j++) {
                PointIndex hi = (j == 1) ? seg[0] : seg[1];
                if (points[hi].Type() == INNERPOINT || points[hi].Type() == SURFACEPOINT) {
                    points[hi].SetType(EDGEPOINT);
                }
            }
        }

        for (int i = 0; i < lockedpoints.size(); i++) {
            points[lockedpoints[i]].SetType(FIXEDPOINT);
        }

        for (int i = 0; i < GetNSeg(); i++) {
            const Segment& seg = segments[i];
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();

            segmentht->Set(i2, i);
        }
    }

    void Mesh::FindOpenElements(int dom)
    {
        size_t np = GetNP();
        size_t nse = GetNSE();

        std::vector<bool> hasface(GetNFD());

        for (int i = 0; i < GetNFD(); i++) {
            int domin = GetFaceDescriptor(i + 1).DomainIn();
            int domout = GetFaceDescriptor(i + 1).DomainOut();
            hasface[i] = (dom == 0 && (domin != 0 || domout != 0)) ||
                         (dom != 0 && (domin == dom || domout == dom));
        }

        std::vector<PointIndex> numonpoint(np, 0);

        for (size_t sii = 0; sii < nse; sii++) {
            int ind = surfelements[sii].GetIndex();

            if (hasface[ind - 1]) {
                const Element2d& hel = surfelements[sii];
                int mini = 0;
                for (size_t j = 1; j < 3; j++) {
                    if (hel[j] < hel[mini]) {
                        mini = j;
                    }
                }
                numonpoint[hel[mini]]++;
            }
        }

        TABLE<SurfaceElementIndex> sels_on_point(numonpoint);
        for (size_t sii = 0; sii < nse; sii++) {
            int ind = surfelements[sii].GetIndex();

            if (hasface[ind - 1]) {
                const Element2d& hel = surfelements[sii];
                int mini = 0;
                for (size_t j = 1; j < 3; j++) {
                    if (hel[j] < hel[mini]) {
                        mini = j;
                    }
                }
                sels_on_point.Add(hel[mini], sii);
            }
        }

        Element2d hel;

        INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht(100);
        openelements.resize(0);

        for (PointIndex pi = 0; pi < points.size(); pi++) {
            if (sels_on_point[pi].size()) {

                faceht.SetSize(2 * sels_on_point[pi].size());

                FlatArray<SurfaceElementIndex> row = sels_on_point[pi];
                for (size_t ii = 0; ii < row.size(); ii++) {
                    hel = SurfaceElement(row[ii]);
                    int ind = hel.GetIndex();

                    if (GetFaceDescriptor(ind).DomainIn() &&
                        (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())) {
                        hel.NormalizeNumbering();
                        if (hel.PNum(1) == pi) {
                            INDEX_3 i3(hel[0], hel[1], hel[2]);
                            INDEX_2 i2(GetFaceDescriptor(ind).DomainIn(), PointIndex{-1});
                            faceht.Set(i3, i2);
                        }
                    }
                    if (GetFaceDescriptor(ind).DomainOut() &&
                        (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())) {
                        hel.Invert();
                        hel.NormalizeNumbering();
                        if (hel.PNum(1) == pi) {
                            INDEX_3 i3(hel[0], hel[1], hel[2]);
                            INDEX_2 i2(GetFaceDescriptor(ind).DomainOut(), PointIndex{-1});
                            faceht.Set(i3, i2);
                        }
                    }
                }
                for (size_t i = 0; i < faceht.Size(); i++) {
                    if (faceht.UsedPos(i)) {
                        INDEX_3 i3;
                        INDEX_2 i2;
                        faceht.GetData(i, i3, i2);
                        if (i2.I1() != -1) {
                            Element2d tri;
                            for (size_t l = 0; l < 3; l++) {
                                tri[l] = i3.I(l + 1);
                            }
                            tri.SetIndex(i2.I1());
                            openelements.push_back(tri);
                        }
                    }
                }
            }
        }

        int cnt3 = 0;
        for (size_t i = 0; i < openelements.size(); i++) {
            cnt3++;
        }
        int cnt4 = openelements.size() - cnt3;
        if (openelements.size() > 0) {
            MESHIT_LOG_WARNING(openelements.size() << " (" << cnt3 << " + " << cnt4 << ")" << " open elements");
        }

        BuildBoundaryEdges();

        for (size_t i = 0; i < openelements.size(); i++) {
            const Element2d& sel = openelements[i];
            if (boundaryedges) {
                for (size_t j = 1; j <= 3; j++) {
                    INDEX_2 i2;
                    i2.I1() = sel.PNumMod(j);
                    i2.I2() = sel.PNumMod(j + 1);
                    i2.Sort();
                    boundaryedges->Set(i2, 1);
                }
            }
            for (size_t j = 1; j <= 3; j++) {
                PointIndex pi = sel.PNum(j);
                if (pi < points.size()) {
                    points[pi].SetType(FIXEDPOINT);
                }
            }
        }
    }

    void Mesh::FindOpenSegments(int surfnr)
    {
        INDEX_2_HASHTABLE<INDEX_2> faceht(4 * GetNSE() + GetNSeg() + 1);

        MESHIT_LOG_DEBUG("Test Opensegments");
        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);

            if (surfnr == 0 || seg.si == surfnr) {
                INDEX_2 key(seg[0], seg[1]);
                INDEX_2 data(seg.si, -(i + 1));

                if (faceht.Used(key)) {
                    std::cerr << "ERROR: Segment " << seg << " already used" << std::endl;
                }

                faceht.Set(key, data);
            }
        }

        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);

            if (surfnr == 0 || seg.si == surfnr) {
                INDEX_2 key(seg[1], seg[0]);
                if (!faceht.Used(key)) {
                    std::cerr << "ERROR: Segment " << seg << " brother not used" << std::endl;
                }
            }
        }

        for (size_t i = 0; i < GetNSE(); i++) {
            const Element2d& el = SurfaceElement(i + 1);
            if (el.IsDeleted()) continue;

            if (surfnr == 0 || el.GetIndex() == surfnr) {
                for (size_t j = 1; j <= 3; j++) {
                    INDEX_2 seg(el.PNumMod(j), el.PNumMod(j + 1));
                    INDEX_2 data;

                    if (seg.I1() <= 0 || seg.I2() <= 0) {
                        std::cerr << "seg = " << seg << std::endl;
                    }
                    if (faceht.Used(seg)) {
                        data = faceht.Get(seg);
                        if (data.I1() == el.GetIndex()) {
                            data.I1() = 0;
                            faceht.Set(seg, data);
                        } else {
                            MESHIT_LOG_WARNING("hash table si not fitting for segment: " <<
                                               seg.I1() << "-" << seg.I2() << " other = " << data.I2());
                        }
                    } else {
                        std::swap(seg.I1(), seg.I2());
                        data.I1() = el.GetIndex();
                        data.I2() = i + 1;
                        faceht.Set(seg, data);
                    }
                }
            }
        }

        std::cerr << "open segments: " << std::endl;
        opensegments.resize(0);
        for (size_t i = 0; i < faceht.GetNBags(); i++) {
            for (size_t j = 0; j < faceht.GetBagSize(i); j++) {
                INDEX_2 i2;
                INDEX_2 data;
                faceht.GetData(i, j, i2, data);
                if (data.I1())  // surfnr
                {
                    Segment seg;
                    seg[0] = i2.I1();
                    seg[1] = i2.I2();
                    seg.si = data.I1();

                    std::cerr << seg[0] << " - " << seg[1]
                    << " len = " << Dist(Point(seg[0]), Point(seg[1]))
                    << std::endl;

                    opensegments.push_back(seg);
                }
            }
        }
        MESHIT_LOG_DEBUG(opensegments.size() << " open segments found");

        for (size_t i = 0; i < points.size(); i++) {
            points[i].SetType(SURFACEPOINT);
        }
        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            points[seg[0]].SetType(EDGEPOINT);
            points[seg[1]].SetType(EDGEPOINT);
        }
        for (size_t i = 0; i < GetNOpenSegments(); i++) {
            const Segment& seg = GetOpenSegment(i);
            points[seg[0]].SetType(EDGEPOINT);
            points[seg[1]].SetType(EDGEPOINT);
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
        if (hloc < hmin)
            hloc = hmin;

        assert(lochfunc);
        lochfunc->SetH(p, hloc);
    }

    void Mesh::RestrictLocalHLine(const Point3d& p1, const Point3d& p2, double hloc)
    {
        if (hloc < hmin) {
            hloc = hmin;
        }
        int steps = static_cast<int>(Dist(p1, p2) / hloc) + 2;
        Vec3d v(p1, p2);

        for (int i = 0; i <= steps; i++) {
            Point3d p = p1 + (static_cast<double>(i) / static_cast<double>(steps) * v);
            RestrictLocalH(p, hloc);
        }
    }

    void Mesh::SetMinimalH(double h)
    {
        hmin = h;
    }

    void Mesh::SetGlobalH(double h)
    {
        hglob = h;
    }

    double Mesh::MaxHDomain(int dom) const
    {
        if (maxhdomain.size())
            return maxhdomain[dom - 1];
        else
            return 1e10;
    }

    void Mesh::SetMaxHDomain(const std::vector<double>& mhd)
    {
        maxhdomain.resize(mhd.size());
        for (int i = 0; i < mhd.size(); i++) {
            maxhdomain[i] = mhd[i];
        }
    }

    double Mesh::GetH(const Point3d& p) const
    {
        double hmin = hglob;
        if (lochfunc) {
            double hl = lochfunc->GetH(p);
            if (hl < hglob)
                hmin = hl;
        }
        return hmin;
    }

    double Mesh::GetMinH(const Point3d& pmin, const Point3d& pmax)
    {
        double hmin = hglob;
        if (lochfunc) {
            double hl = lochfunc->GetMinH(pmin, pmax);
            if (hl < hmin)
                hmin = hl;
        }
        return hmin;
    }

    double Mesh::AverageH(int surfnr) const
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
            const Element2d& el = surfelements[i];
            double hel = -1;
            for (size_t j = 1; j <= 3; j++) {
                const Point3d& p1 = points[el.PNumMod(j)];
                const Point3d& p2 = points[el.PNumMod(j + 1)];

                if (!ident->UsedSymmetric(el.PNumMod(j), el.PNumMod(j + 1))) {
                    double hedge = Dist(p1, p2);
                    if (hedge > hel)
                        hel = hedge;
                }
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

            if (!ident->UsedSymmetric(seg[0], seg[1])) {
                lochfunc->SetH(Center(p1, p2), Dist(p1, p2));
            }
        }
    }

    void Mesh::CalcLocalHFromPointDistances()
    {
        MESHIT_LOG_DEBUG("Calculating local h from point distances");

        assert(lochfunc);

        for (size_t i = 0; i < GetNP(); i++) {
            for (size_t j = i + 1; j < GetNP(); j++) {
                const Point3d& p1 = points[i];
                const Point3d& p2 = points[j];
                double hl = Dist(p1, p2);
                RestrictLocalH(p1, hl);
                RestrictLocalH(p2, hl);
            }
        }
    }

    void Mesh::CalcLocalHFromSurfaceCurvature(double elperr)
    {
        MESHIT_LOG_DEBUG("Calculating local h from surface curvature");

        assert(lochfunc);

        INDEX_2_HASHTABLE<int> edges(3 * GetNP() + 2);
        INDEX_2_HASHTABLE<int> bedges(GetNSeg() + 2);
        int j;

        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();
            bedges.Set(i2, 1);
        }
        for (size_t i = 0; i < GetNSE(); i++) {
            const Element2d& sel = SurfaceElement(i);
            if (!sel.PNum(1))
                continue;
            for (j = 1; j <= 3; j++) {
                INDEX_2 i2(sel.PNumMod(j), sel.PNumMod(j + 1));
                i2.Sort();
                if (bedges.Used(i2)) continue;

                if (edges.Used(i2)) {
                    int other = edges.Get(i2);

                    const Element2d& elother = SurfaceElement(other - 1);

                    int pi3 = 1;
                    while ((sel.PNum(pi3) == i2.I1()) ||
                           (sel.PNum(pi3) == i2.I2()))
                        pi3++;
                    pi3 = sel.PNum(pi3);

                    int pi4 = 1;
                    while ((elother.PNum(pi4) == i2.I1()) ||
                           (elother.PNum(pi4) == i2.I2()))
                        pi4++;
                    pi4 = elother.PNum(pi4);

                    double rad = ComputeCylinderRadius(points[i2.I1() - 1],
                                                       points[i2.I2() - 1],
                                                       points[pi3 - 1],
                                                       points[pi4 - 1]);

                    RestrictLocalHLine(points[i2.I1() - 1], points[i2.I2() - 1], rad / elperr);
                }
                else
                    edges.Set(i2, i + 1);
            }
        }


        // Restrict h due to line segments

        for (size_t i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            const Point3d& p1 = Point(seg[0]);
            const Point3d& p2 = Point(seg[1]);
            RestrictLocalH(Center(p1, p2), Dist(p1, p2));
        }
    }

    void Mesh::RestrictLocalH(resthtype rht, int nr, double loch)
    {
        switch (rht) {
            case RESTRICTH_FACE: {
                for (size_t i = 0; i < GetNSE(); i++) {
                    const Element2d& sel = SurfaceElement(i);
                    if (sel.GetIndex() == nr)
                        RestrictLocalH(RESTRICTH_SURFACEELEMENT, i + 1, loch);
                }
                break;
            }
            case RESTRICTH_EDGE: {
                for (size_t i = 0; i < GetNSeg(); i++) {
                    const Segment& seg = LineSegment(i);
                    if (seg.edgenr == nr)
                        RestrictLocalH(RESTRICTH_SEGMENT, i + 1, loch);
                }
                break;
            }
            case RESTRICTH_POINT: {
                RestrictLocalH(points[nr - 1], loch);
                break;
            }

            case RESTRICTH_SURFACEELEMENT: {
                const Element2d& sel = SurfaceElement(nr - 1);
                Point3d p = Center(Point(sel.PNum(1)),
                                   Point(sel.PNum(2)),
                                   Point(sel.PNum(3)));
                RestrictLocalH(p, loch);
                break;
            }
            case RESTRICTH_SEGMENT: {
                const Segment& seg = LineSegment(nr + 1);
                RestrictLocalHLine(Point(seg[0]), Point(seg[1]), loch);
                break;
            }
        }
    }

    void Mesh::LoadLocalMeshSize(const char* meshsizefilename)
    {
        // Philippose - 10/03/2009
        // Improve error checking when loading and reading
        // the local mesh size file

        if (!meshsizefilename) return;

        std::ifstream msf(meshsizefilename);

        // Philippose - 09/03/2009
        // Adding print message information in case the specified 
        // does not exist, or does not load successfully due to 
        // other reasons such as access rights, etc...
        if (!msf) {
            MESHIT_LOG_ERROR("Error loading mesh size file: " << meshsizefilename << "....  Skipping!");
            return;
        }

        MESHIT_LOG_DEBUG("Load local mesh-size file: " << meshsizefilename);

        int nmsp = 0;
        int nmsl = 0;

        msf >> nmsp;
        if (!msf.good())
            throw std::runtime_error("Mesh-size file error: No points found\n");

        if (nmsp > 0)
            MESHIT_LOG_DEBUG("Number of mesh-size restriction points: " << nmsp);

        for (int i = 0; i < nmsp; i++) {
            Point3d pi;
            double hi;
            msf >> pi.X() >> pi.Y() >> pi.Z();
            msf >> hi;
            if (!msf.good())
                throw std::runtime_error("Mesh-size file error: Number of points don't match specified list size\n");
            RestrictLocalH(pi, hi);
        }

        msf >> nmsl;
        if (!msf.good())
            throw std::runtime_error("Mesh-size file error: No line definitions found\n");

        if (nmsl > 0)
            MESHIT_LOG_DEBUG("Number of mesh-size restriction lines: " << nmsl);

        for (int i = 0; i < nmsl; i++) {
            Point3d p1, p2;
            double hi;
            msf >> p1.X() >> p1.Y() >> p1.Z();
            msf >> p2.X() >> p2.Y() >> p2.Z();
            msf >> hi;
            if (!msf.good())
                throw std::runtime_error(
                        "Mesh-size file error: Number of line definitions don't match specified list size\n");
            RestrictLocalHLine(p1, p2, hi);
        }

        msf.close();
    }

    void Mesh::GetBox(Point3d& pmin, Point3d& pmax, int dom) const
    {
        if (points.size() == 0) {
            pmin = pmax = Point3d(0, 0, 0);
            return;
        }

        pmin = Point3d(1e10, 1e10, 1e10);
        pmax = Point3d(-1e10, -1e10, -1e10);
        if (dom <= 0) {
            for (size_t pi = 0; pi < points.size(); pi++) {
                pmin.SetToMin(points[pi]);
                pmax.SetToMax(points[pi]);
            }
        } else {
            for (size_t sei = 0; sei < GetNSE(); sei++) {
                const Element2d& el = surfelements[sei];
                if (el.IsDeleted()) continue;

                if (dom == -1 || el.GetIndex() == dom) {
                    for (size_t j = 0; j < 3; j++) {
                        pmin.SetToMin(points[el[j]]);
                        pmax.SetToMax(points[el[j]]);
                    }
                }
            }
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

        for (size_t i = 0; i < surfelements.size(); i++) {
            if (surfelements[i].IsDeleted()) {
                surfelements.erase(surfelements.begin() + i);
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

        for (size_t i = 0; i < surfelements.size(); i++) {
            const Element2d& el = surfelements[i];
            for (size_t j = 0; j < 3; j++) {
                pused.Set(el[j]);
            }
        }

        for (size_t i = 0; i < segments.size(); i++) {
            const Segment& seg = segments[i];
            pused.Set(seg[0]);
            pused.Set(seg[1]);
        }

        for (size_t i = 0; i < openelements.size(); i++) {
            const Element2d& el = openelements[i];
            for (size_t j = 0; j < 3; j++) {
                pused.Set(el[j]);
            }
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

        for (size_t i = 0; i < surfelements.size(); i++) {
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

        for (size_t i = 0; i < openelements.size(); i++) {
            Element2d& el = openelements[i];
            for (size_t j = 0; j < 3; j++) {
                el[j] = op2np[el[j]];
            }
        }

        for (size_t i = 0; i < lockedpoints.size(); i++) {
            lockedpoints[i] = op2np[lockedpoints[i]];
        }

        for (size_t i = 0; i < facedecoding.size(); i++) {
            facedecoding[i].firstelement = -1;
        }
        for (int i = surfelements.size() - 1; i >= 0; i--) {
            int ind = surfelements[i].GetIndex();
            surfelements[i].next = facedecoding[ind - 1].firstelement;
            facedecoding[ind - 1].firstelement = i;
        }

        CalcSurfacesOfNode();

        //  FindOpenElements();
        timestamp = NextTimeStamp();
    }

    int Mesh::CheckConsistentBoundary() const
    {
        size_t nf = GetNOpenElements();
        INDEX_2_HASHTABLE<int> edges(nf + 2);
        INDEX_2 i2, i2s, edge;
        int err = 0;

        for (size_t i = 0; i < nf; i++) {
            const Element2d& sel = OpenElement(i);

            for (size_t j = 1; j <= 3; j++) {
                i2.I1() = sel.PNumMod(j);
                i2.I2() = sel.PNumMod(j + 1);

                int sign = (i2.I2() > i2.I1()) ? 1 : -1;
                i2.Sort();
                if (!edges.Used(i2)) {
                    edges.Set(i2, 0);
                }
                edges.Set(i2, edges.Get(i2) + sign);
            }
        }

        for (size_t i = 0; i < edges.GetNBags(); i++) {
            for (size_t j = 0; j < edges.GetBagSize(i); j++) {
                int cnt = 0;
                edges.GetData(i, j, i2, cnt);
                if (cnt) {
                    MESHIT_LOG_ERROR("Edge " << i2.I1() << " - " << i2.I2() << " multiple times in surface mesh");
                    i2s = i2;
                    i2s.Sort();
                    for (size_t k = 0; k < nf; k++) {
                        const Element2d& sel = OpenElement(k);
                        for (size_t l = 1; l <= 3; l++) {
                            edge.I1() = sel.PNumMod(l);
                            edge.I2() = sel.PNumMod(l + 1);
                            edge.Sort();

                            if (edge == i2s) {
                                MESHIT_LOG_ERROR("  edge of element " << sel);
                            }
                        }
                    }
                    err = 2;
                }
            }
        }
        return err;
    }

    int Mesh::CheckOverlappingBoundary()
    {
        Point3d pmin, pmax;
        GetBox(pmin, pmax);
        Box3dTree setree(pmin, pmax);
        std::vector<size_t> inters;

        bool overlap = 0;
        bool incons_layers = 0;

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
            const Element2d& tri = SurfaceElement(i);

            Point3d tpmin(Point(tri[0]));
            Point3d tpmax(tpmin);

            for (size_t k = 1; k < 3; k++) {
                tpmin.SetToMin(Point(tri[k]));
                tpmax.SetToMax(Point(tri[k]));
            }

            setree.GetIntersecting(tpmin, tpmax, inters);

            for (size_t j = 0; j < inters.size(); j++) {
                const Element2d& tri2 = SurfaceElement(inters[j] - 1);

                if (points[tri[0]].GetLayer() != points[tri2[0]].GetLayer())
                    continue;

                if (points[tri[0]].GetLayer() != points[tri[1]].GetLayer() ||
                    points[tri[0]].GetLayer() != points[tri[2]].GetLayer()) {
                    incons_layers = 1;
                    MESHIT_LOG_WARNING("inconsistent layers in triangle");
                }

                const meshit::Point3d* trip1[3], * trip2[3];
                for (size_t k = 1; k <= 3; k++) {
                    trip1[k - 1] = &Point(tri.PNum(k));
                    trip2[k - 1] = &Point(tri2.PNum(k));
                }

                if (IntersectTriangleTriangle(&trip1[0], &trip2[0])) {
                    overlap = 1;
                    MESHIT_LOG_WARNING("Intersecting elements " << i + 1 << " and " << inters[j]);
                    MESHIT_LOG_DEBUG(" el1 = " << tri);
                    MESHIT_LOG_DEBUG(" el2 = " << tri2);

                    for (size_t k = 1; k <= 3; k++)
                        MESHIT_LOG_DEBUG_CONT(tri.PNum(k) << "  ");
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

                    MESHIT_LOG_DEBUG("Face1 = " << GetFaceDescriptor(tri.GetIndex()));
                    MESHIT_LOG_DEBUG("Face1 = " << GetFaceDescriptor(tri2.GetIndex()));
                }
            }
        }

        // bug 'fix'
        if (incons_layers) overlap = 0;

        return overlap;
    }

    int Mesh::GetNDomains() const
    {
        int ndom = 0;

        for (size_t k = 0; k < facedecoding.size(); k++) {
            if (facedecoding[k].DomainIn() > ndom)
                ndom = facedecoding[k].DomainIn();
            if (facedecoding[k].DomainOut() > ndom)
                ndom = facedecoding[k].DomainOut();
        }

        return ndom;
    }

    bool Mesh::PointContainedIn2DElement(const Point3d& p,
                                         double lami[3],
                                         const int element,  // 0-based
                                         bool consider3D) const
    {
        const double eps = 1e-6;

        const Element2d& el = SurfaceElement(element);
        const Point3d& p1 = Point(el.PNum(1));
        const Point3d& p2 = Point(el.PNum(2));
        const Point3d& p3 = Point(el.PNum(3));

        Vec3d col1 = p2 - p1;
        Vec3d col2 = p3 - p1;
        Vec3d col3 = Cross(col1, col2);
        Vec3d rhs = p - p1;
        Vec3d sol;

        SolveLinearSystem(col1, col2, col3, rhs, sol);

        if (sol.X() >= -eps && sol.Y() >= -eps &&
            sol.X() + sol.Y() <= 1 + eps) {
            if (!consider3D || (sol.Z() >= -eps && sol.Z() <= eps)) {
                lami[0] = sol.X();
                lami[1] = sol.Y();
                lami[2] = sol.Z();

                return true;
            }
        }
        return false;
    }

    void Mesh::RebuildSurfaceElementLists()
    {
        for (size_t i = 0; i < facedecoding.size(); i++) {
            facedecoding[i].firstelement = -1;
        }
        for (int i = surfelements.size() - 1; i >= 0; i--) {
            int ind = surfelements[i].GetIndex();
            surfelements[i].next = facedecoding[ind - 1].firstelement;
            facedecoding[ind - 1].firstelement = i;
        }
    }

    void Mesh::GetSurfaceElementsOfFace(int facenr, std::vector<SurfaceElementIndex>& sei) const
    {
        sei.resize(0);

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
        for (size_t i = 0; i < surfelements.size(); i++) {
            const Element2d& el = SurfaceElement(i);
            for (size_t j = 1; j <= 3; j++) {
                if (el.PNum(j) > numvertices)
                    numvertices = el.PNum(j);
            }
        }
        numvertices += 1;
    }

    size_t Mesh::GetNV() const
    {
        if (numvertices < 0)
            return GetNP();
        else
            return static_cast<size_t>(numvertices);
    }

    void Mesh::SetNP(int np)
    {
        points.resize(np);
        //  ptyps.SetSize(np);

        int mlold = mlbetweennodes.size();
        mlbetweennodes.resize(np);
        if (np > mlold) {
            for (int i = mlold; i < np; i++) {
                mlbetweennodes[i].I1() = -1;
                mlbetweennodes[i].I2() = -1;
            }
        }
        GetIdentifications().SetMaxPointNr(np - 1);
    }

    void Mesh::UpdateTopology()
    {
        topology->Update();
    }

    void Mesh::SetMaterial(int domnr, const char* mat)
    {
        if (domnr > materials.size()) {
            int olds = materials.size();
            materials.resize(domnr);
            for (int i = olds; i < domnr; i++) {
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
        surfacesonnode.PrintMemInfo(ost);

        ost << "boundaryedges: ";
        if (boundaryedges)
            boundaryedges->PrintMemInfo(ost);

        ost << "surfelementht: ";
        if (surfelementht)
            surfelementht->PrintMemInfo(ost);
    }
}  // namespace meshit
