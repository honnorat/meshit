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
        elementsearchtree = NULL;
        elementsearchtreets = NextTimeStamp();
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
        delete elementsearchtree;

        for (size_t i = 0; i < materials.size(); i++) {
            delete[] materials[i];
        }
        for (size_t i = 0; i < userdata_int.Size(); i++) {
            delete userdata_int[i];
        }
        for (size_t i = 0; i < userdata_double.Size(); i++) {
            delete userdata_double[i];
        }
        for (size_t i = 0; i < bcnames.size(); i++) {
            if (bcnames[i]) delete bcnames[i];
        }
    }

    Mesh& Mesh::operator=(const Mesh& mesh2)
    {
        points = mesh2.points;
        segments = mesh2.segments;
        surfelements = mesh2.surfelements;
        lockedpoints = mesh2.lockedpoints;
        facedecoding = mesh2.facedecoding;

        bcnames.resize(mesh2.bcnames.size());
        for (size_t i = 0; i < mesh2.bcnames.size(); i++) {
            if (mesh2.bcnames[i]) {
                bcnames[i] = new std::string(*mesh2.bcnames[i]);
            } else {
                bcnames[i] = nullptr;
            }
        }

        return *this;
    }

    void Mesh::BuildFromSpline2D(SplineGeometry2d& geometry, MeshingParameters& mp)
    {
        MESHIT_LOG_DEBUG("Generate Mesh from spline geometry");

        double h = mp.maxh;

        Box<2> bbox = geometry.GetBoundingBox();

        if (bbox.Diam() < h) {
            h = bbox.Diam();
            mp.maxh = h;
        }

        Point3d pmin(bbox.PMin()(0), bbox.PMin()(1), -bbox.Diam());
        Point3d pmax(bbox.PMax()(0), bbox.PMax()(1), bbox.Diam());

        SetLocalH(pmin, pmax, mp.grading);
        SetGlobalH(h);

        geometry.PartitionBoundary(mp, h, *this);

        // marks mesh points for hp-refinement
        for (size_t i = 0; i < geometry.GetNP(); i++) {
            if (geometry.GetPoint(i).hpref) {
                double mindist = 1e99;
                size_t mpi = 0;
                ::meshit::Point<2> gp = geometry.GetPoint(i);
                ::meshit::Point<3> gp3(gp(0), gp(1), 0);
                for (size_t pi = 0; pi < GetNP(); pi++) {
                    if (Dist2(gp3, points[pi]) < mindist) {
                        mpi = pi;
                        mindist = Dist2(gp3, points[pi]);
                    }
                }
                points[mpi].Singularity(1.);
            }
        }

        int maxdomnr = 0;
        for (size_t si = 0; si < GetNSeg(); si++) {
            if (segments[si].domin > maxdomnr) maxdomnr = segments[si].domin;
            if (segments[si].domout > maxdomnr) maxdomnr = segments[si].domout;
        }

        ClearFaceDescriptors();
        for (int i = 1; i <= maxdomnr; i++) {
            AddFaceDescriptor(FaceDescriptor(i, 0, 0, i));
        }

        int maxsegmentindex = 0;
        for (size_t si = 0; si < GetNSeg(); si++) {
            if (segments[si].si > maxsegmentindex) maxsegmentindex = segments[si].si;
        }

        SetNBCNames(maxsegmentindex);

        for (int sindex = 0; sindex < maxsegmentindex; sindex++) {
            SetBCName(sindex, geometry.GetBCName(sindex + 1));
        }

        for (size_t si = 0; si < GetNSeg(); si++) {
            segments[si].SetBCName(bcnames[segments[si].si - 1]);
        }

        CalcLocalH();

        int bnp = GetNP();  // boundary points
        int hquad = mp.quad;

        for (int domnr = 1; domnr <= maxdomnr; domnr++) {
            if (!geometry.GetDomainTensorMeshing(domnr))
                continue;

            // tensor product mesh
            Array<PointIndex> nextpi(bnp);
            Array<int> si1(bnp), si2(bnp);

            nextpi = -1;
            si1 = -1;
            si2 = -1;
            for (size_t si = 0; si < GetNSeg(); si++) {
                int p1 = -1, p2 = -2;

                if (segments[si].domin == domnr) {
                    p1 = segments[si][0];
                    p2 = segments[si][1];
                }
                if (segments[si].domout == domnr) {
                    p1 = segments[si][1];
                    p2 = segments[si][0];
                }

                if (p1 == -1) continue;

                nextpi[p1] = p2;  // counter-clockwise

                int index = segments[si].si;
                if (si1[p1] != index && si2[p1] != index) {
                    si2[p1] = si1[p1];
                    si1[p1] = index;
                }
                if (si1[p2] != index && si2[p2] != index) {
                    si2[p2] = si1[p2];
                    si1[p2] = index;
                }
            }

            PointIndex c1(0), c2, c3, c4;  // 4 corner points
            size_t nex = 1, ney = 1;

            for (PointIndex pi = 1; pi <= si2.size(); pi++) {
                if (si2[pi] != -1) {
                    c1 = pi;
                    break;
                }
            }
            for (c2 = nextpi[c1]; si2[c2] == -1; c2 = nextpi[c2], nex++) { }
            for (c3 = nextpi[c2]; si2[c3] == -1; c3 = nextpi[c3], ney++) { }
            for (c4 = nextpi[c3]; si2[c4] == -1; c4 = nextpi[c4]) { }

            Array<PointIndex> pts((nex + 1) * (ney + 1));  // x ... inner loop
            pts = -1;

            for (PointIndex pi = c1, i = 0; pi != c2; pi = nextpi[pi], i++) {
                pts[i] = pi;
            }
            for (PointIndex pi = c2, i = 0; pi != c3; pi = nextpi[pi], i++) {
                pts[(nex + 1) * i + nex] = pi;
            }
            for (PointIndex pi = c3, i = 0; pi != c4; pi = nextpi[pi], i++) {
                pts[(nex + 1) * (ney + 1) - i - 1] = pi;
            }
            for (PointIndex pi = c4, i = 0; pi != c1; pi = nextpi[pi], i++) {
                pts[(nex + 1) * (ney - i)] = pi;
            }

            for (PointIndex pix = nextpi[c1], ix = 0; pix != c2; pix = nextpi[pix], ix++) {
                for (PointIndex piy = nextpi[c2], iy = 0; piy != c3; piy = nextpi[piy], iy++) {
                    ::meshit::Point<3> p = points[pix] + (points[piy] - points[c2]);
                    pts[(nex + 1) * (iy + 1) + ix + 1] = AddPoint(p, 1, FIXEDPOINT);
                }
            }

            for (size_t i = 0; i < ney; i++) {
                for (size_t j = 0; j < nex; j++) {
                    Element2d el(QUAD);
                    el[0] = pts[i * (nex + 1) + j];
                    el[1] = pts[i * (nex + 1) + j + 1];
                    el[2] = pts[(i + 1) * (nex + 1) + j + 1];
                    el[3] = pts[(i + 1) * (nex + 1) + j];
                    el.SetIndex(domnr);
                    AddSurfaceElement(el);
                }
            }
        }

        for (int domnr = 1; domnr <= maxdomnr; domnr++) {
            if (geometry.GetDomainTensorMeshing(domnr))
                continue;

            if (geometry.GetDomainMaxh(domnr) > 0)
                h = geometry.GetDomainMaxh(domnr);

            MESHIT_LOG_DEBUG("Meshing domain " << domnr << " / " << maxdomnr);

            int oldnf = GetNSE();

            mp.quad = hquad || geometry.GetDomainQuadMeshing(domnr);

            Meshing2 meshing(mp, Box<3>(pmin, pmax));

            Array<int> compress(bnp);
            compress = -1;
            int cnt = 0;
            for (PointIndex pi = 0; pi < bnp; pi++) {
                if (points[pi].GetLayer() == geometry.GetDomainLayer(domnr)) {
                    meshing.AddPoint(points[pi], pi);
                    cnt++;
                    compress[pi] = cnt;
                }
            }
            PointGeomInfo gi;
            gi.trignum = 1;
            for (size_t si = 0; si < GetNSeg(); si++) {
                if (segments[si].domin == domnr) {
                    meshing.AddBoundaryElement(
                            compress[segments[si][0]],
                            compress[segments[si][1]], gi, gi);
                }
                if (segments[si].domout == domnr) {
                    meshing.AddBoundaryElement(
                            compress[segments[si][1]],
                            compress[segments[si][0]], gi, gi);
                }
            }

            mp.checkoverlap = 0;
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

        mp.quad = hquad;

        int hsteps = mp.optsteps2d;

        mp.optimize2d = "smcm";
        mp.optsteps2d = hsteps / 2;
        Optimize2d(*this, mp);

        mp.optimize2d = "Smcm";
        mp.optsteps2d = (hsteps + 1) / 2;
        Optimize2d(*this, mp);

        mp.optsteps2d = hsteps;

        Compress();
    }

    PointIndex Mesh::AddPoint(const Point3d& p, int layer)
    {
        return AddPoint(p, layer, INNERPOINT);
    }

    PointIndex Mesh::AddPoint(const Point3d& p, int layer, POINTTYPE type)
    {
        timestamp = NextTimeStamp();

        PointIndex pi = points.End();
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

        int maxn = el[0];
        for (size_t i = 1; i < el.GetNP(); i++) {
            if (el[i] > maxn) maxn = el[i];
        }
        maxn += 1;

        if (maxn <= points.size()) {
            for (size_t i = 0; i < el.GetNP(); i++) {
                if (points[el[i]].Type() > SURFACEPOINT)
                    points[el[i]].SetType(SURFACEPOINT);
            }
        }

        size_t si = surfelements.size();
        surfelements.push_back(el);

        if (el.index > facedecoding.size()) {
            MESHIT_LOG_ERROR("has no facedecoding: fd.size = " << facedecoding.size() << ", ind = " << el.index);
        }

        surfelements.Last().next = facedecoding[el.index - 1].firstelement;
        facedecoding[el.index - 1].firstelement = si;

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
        outfile << "dimension\n" << GetDimension() << "\n";
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

            outfile << " " << sel.GetNP();
            for (size_t j = 0; j < sel.GetNP(); j++) {
                outfile << " " << sel[j];
            }
            outfile << "\n";
        }

        outfile << "\n" << "\n";
        outfile << "# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2   ";
        outfile << "ednr1   dist1   ednr2   dist2 \n";
        outfile << "edgesegmentsgi2" << "\n";
        outfile << GetNSeg() << "\n";

        for (int i = 0; i < GetNSeg(); i++) {
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
            outfile << seg.geominfo[0].trignum;  // stl dreiecke
            outfile << " ";
            outfile.width(8);
            outfile << seg.geominfo[1].trignum;  // <<std::endl;  // stl dreieck

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
            Array<INDEX_2> identpairs;
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

        int cntbcnames = 0;
        for (int ii = 0; ii < bcnames.size(); ii++) {
            if (bcnames[ii]) cntbcnames++;
        }

        if (cntbcnames) {
            outfile << "\n\nbcnames" << std::endl << bcnames.size() << std::endl;
            for (int i = 0; i < bcnames.size(); i++) {
                outfile << i + 1 << "\t" << GetBCName(i) << std::endl;
            }
            outfile << std::endl << std::endl;
        }

        int cnt_sing = 0;
        for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
            if (points[pi].Singularity() >= 1.) cnt_sing++;
        }

        if (cnt_sing) {
            outfile << "singular_points" << std::endl << cnt_sing << std::endl;
            for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
                if (points[pi].Singularity() >= 1.)
                    outfile << int(pi) << "\t" << points[pi].Singularity() << std::endl;
            }
        }

        cnt_sing = 0;
        for (size_t si = 0; si < GetNSeg(); si++) {
            if (segments[si].singedge_left) cnt_sing++;
        }
        if (cnt_sing) {
            outfile << "singular_edge_left" << std::endl << cnt_sing << std::endl;
            for (size_t si = 0; si < GetNSeg(); si++) {
                if (segments[si].singedge_left) {
                    outfile << si << "\t" << segments[si].singedge_left << std::endl;
                }
            }
        }

        cnt_sing = 0;
        for (size_t si = 0; si < GetNSeg(); si++) {
            if (segments[si].singedge_right) cnt_sing++;
        }
        if (cnt_sing) {
            outfile << "singular_edge_right" << std::endl << cnt_sing << std::endl;
            for (size_t si = 0; si < GetNSeg(); si++) {
                if (segments[si].singedge_right) {
                    outfile << si << "\t" << segments[si].singedge_right << std::endl;
                }
            }
        }


        cnt_sing = 0;
        for (size_t sei = 0; sei < GetNSE(); sei++) {
            if (GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular) {
                cnt_sing++;
            }
        }

        if (cnt_sing) {
            outfile << "singular_face_inside" << std::endl << cnt_sing << std::endl;
            for (size_t sei = 0; sei < GetNSE(); sei++) {
                if (GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular) {
                    outfile << int(sei) << "\t"
                    << GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular << std::endl;
                }
            }
        }

        cnt_sing = 0;
        for (size_t sei = 0; sei < GetNSE(); sei++) {
            if (GetFaceDescriptor(surfelements[sei].GetIndex()).domout_singular) {
                cnt_sing++;
            }
        }
        if (cnt_sing) {
            outfile << "singular_face_outside" << std::endl << cnt_sing << std::endl;
            for (size_t sei = 0; sei < GetNSE(); sei++) {
                if (GetFaceDescriptor(surfelements[sei].GetIndex()).domout_singular) {
                    outfile << sei << "\t";
                    outfile << GetFaceDescriptor(surfelements[sei].GetIndex()).domout_singular << std::endl;
                }
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
                        faceind = AddFaceDescriptor(FaceDescriptor(surfnr, domin, domout, 0));
                        GetFaceDescriptor(faceind).SetBCProperty(bcp);
                    }

                    infile >> nep;
                    if (!nep) nep = 3;

                    Element2d tri(nep);
                    tri.SetIndex(faceind);

                    for (int j = 1; j <= nep; j++) {
                        infile >> tri.PNum(j);
                    }
                    if (geominfo) {
                        for (int j = 1; j <= nep; j++) {
                            infile >> tri.GeomInfoPi(j).trignum;
                        }
                    }
                    if (uv) {
                        for (int j = 1; j <= nep; j++) {
                            infile >> tri.GeomInfoPi(j).u >> tri.GeomInfoPi(j).v;
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
                    infile >> seg.si >> hi >> seg[0] >> seg[1]
                    >> seg.geominfo[0].trignum
                    >> seg.geominfo[1].trignum;
                    AddSegment(seg);
                }
            }
            if (strcmp(str, "edgesegmentsgi2") == 0) {
                int a;
                infile >> a;
                n = a;

                MESHIT_LOG_DEBUG(n << " curve elements");

                for (int i = 1; i <= n; i++) {
                    Segment seg;
                    int hi;
                    infile >> seg.si >> hi >> seg[0] >> seg[1]
                    >> seg.geominfo[0].trignum
                    >> seg.geominfo[1].trignum
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
            if (strcmp(str, "bcnames") == 0) {
                infile >> n;
                MESHIT_LOG_DEBUG(n << " bcnames");
                Array<int> bcnrs(n);
                SetNBCNames(n);
                for (int i = 1; i <= n; i++) {
                    std::string nextbcname;
                    infile >> bcnrs[i - 1] >> nextbcname;
                    bcnames[bcnrs[i - 1] - 1] = new std::string(nextbcname);
                }
                if (GetDimension() == 2) {
                    for (int i = 0; i < GetNSeg(); i++) {
                        Segment& seg = LineSegment(i);
                        if (seg.si <= n)
                            seg.SetBCName(bcnames[seg.si - 1]);
                        else
                            seg.SetBCName(0);
                    }
                } else {
                    for (size_t sei = 0; sei < GetNSE(); sei++) {
                        if (surfelements[sei].GetIndex()) {
                            int bcp = GetFaceDescriptor(surfelements[sei].GetIndex()).BCProperty();
                            if (bcp <= n)
                                GetFaceDescriptor(surfelements[sei].GetIndex()).SetBCName(bcnames[bcp - 1]);
                            else
                                GetFaceDescriptor(surfelements[sei].GetIndex()).SetBCName(0);
                        }
                    }

                }
            }
            if (strcmp(str, "singular_points") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    PointIndex pi;
                    double s;
                    infile >> pi;
                    infile >> s;
                    points[pi].Singularity(s);
                }
            }
            if (strcmp(str, "singular_edge_left") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    size_t si;
                    double s;
                    infile >> si;
                    infile >> s;
                    segments[si].singedge_left = s;
                }
            }
            if (strcmp(str, "singular_edge_right") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    size_t si;
                    double s;
                    infile >> si;
                    infile >> s;
                    segments[si].singedge_right = s;
                }
            }
            if (strcmp(str, "singular_face_inside") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    size_t sei;
                    double s;
                    infile >> sei;
                    infile >> s;
                    GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular = s;
                }
            }
            if (strcmp(str, "singular_face_outside") == 0) {
                infile >> n;
                for (int i = 1; i <= n; i++) {
                    size_t sei;
                    double s;
                    infile >> sei;
                    infile >> s;
                    GetFaceDescriptor(surfelements[sei].GetIndex()).domout_singular = s;
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

            // int si = sel.GetIndex();

            if (sel.GetNP() <= 4)
                for (size_t j = 0; j < sel.GetNP(); j++) {
                    INDEX_2 i2;
                    i2.I1() = sel.PNumMod(j + 1);
                    i2.I2() = sel.PNumMod(j + 2);
                    i2.Sort();
                    boundaryedges->Set(i2, 1);
                }
            else if (sel.GetType() == TRIG6) {
                for (int j = 0; j < 3; j++) {
                    INDEX_2 i2;
                    i2.I1() = sel[j];
                    i2.I2() = sel[(j + 1) % 3];
                    i2.Sort();
                    boundaryedges->Set(i2, 1);
                }
            }
            else
                std::cerr << "illegal elemenet for buildboundaryedges" << std::endl;
        }

        for (int i = 0; i < openelements.size(); i++) {
            const Element2d& sel = openelements[i];
            for (size_t j = 0; j < sel.GetNP(); j++) {
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

            for (size_t j = 0; j < sel.GetNP(); j++) {
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
                if (points[hi].Type() == INNERPOINT ||
                    points[hi].Type() == SURFACEPOINT)
                    points[hi].SetType(EDGEPOINT);
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

        Array<int> numonpoint(np);

        numonpoint = 0;

        Array<char> hasface(GetNFD());

        for (int i = 0; i < GetNFD(); i++) {
            int domin = GetFaceDescriptor(i + 1).DomainIn();
            int domout = GetFaceDescriptor(i + 1).DomainOut();
            hasface[i] =
                    (dom == 0 && (domin != 0 || domout != 0)) ||
                    (dom != 0 && (domin == dom || domout == dom));
        }

        numonpoint = 0;
        for (size_t sii = 0; sii < nse; sii++) {
            int ind = surfelements[sii].GetIndex();

            if (hasface[ind - 1]) {
                const Element2d& hel = surfelements[sii];
                int mini = 0;
                for (size_t j = 1; j < hel.GetNP(); j++) {
                    if (hel[j] < hel[mini]) {
                        mini = j;
                    }
                }
                numonpoint[hel[mini]]++;
            }
        }

        TABLE<SurfaceElementIndex> selsonpoint(numonpoint);
        for (size_t sii = 0; sii < nse; sii++) {
            int ind = surfelements[sii].GetIndex();

            if (hasface[ind - 1]) {
                const Element2d& hel = surfelements[sii];
                int mini = 0;
                for (size_t j = 1; j < hel.GetNP(); j++) {
                    if (hel[j] < hel[mini]) {
                        mini = j;
                    }
                }
                selsonpoint.Add(hel[mini], sii);
            }
        }

        Element2d hel;

        INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht(100);
        openelements.resize(0);

        for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
            if (selsonpoint[pi].size()) {

                faceht.SetSize(2 * selsonpoint[pi].size());

                FlatArray<SurfaceElementIndex> row = selsonpoint[pi];
                for (size_t ii = 0; ii < row.size(); ii++) {
                    hel = SurfaceElement(row[ii]);
                    if (hel.GetType() == TRIG6) hel.SetType(TRIG);
                    int ind = hel.GetIndex();

                    if (GetFaceDescriptor(ind).DomainIn() &&
                        (dom == 0 || dom == GetFaceDescriptor(ind).DomainIn())) {
                        hel.NormalizeNumbering();
                        if (hel.PNum(1) == pi) {
                            INDEX_3 i3(hel[0], hel[1], hel[2]);
                            INDEX_2 i2(GetFaceDescriptor(ind).DomainIn(),
                                       (hel.GetNP() == 3)
                                       ? PointIndex{-1}
                                       : hel.PNum(4));
                            faceht.Set(i3, i2);
                        }
                    }
                    if (GetFaceDescriptor(ind).DomainOut() &&
                        (dom == 0 || dom == GetFaceDescriptor(ind).DomainOut())) {
                        hel.Invert();
                        hel.NormalizeNumbering();
                        if (hel.PNum(1) == pi) {
                            INDEX_3 i3(hel[0], hel[1], hel[2]);
                            INDEX_2 i2(GetFaceDescriptor(ind).DomainOut(),
                                       (hel.GetNP() == 3)
                                       ? PointIndex{-1}
                                       : hel.PNum(4));
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
                            tri.SetType((i2.I2() == -1) ? TRIG : QUAD);
                            for (int l = 0; l < 3; l++) {
                                tri[l] = i3.I(l + 1);
                            }
                            tri.PNum(4) = i2.I2();
                            tri.SetIndex(i2.I1());
                            openelements.push_back(tri);
                        }
                    }
                }
            }
        }

        int cnt3 = 0;
        for (size_t i = 0; i < openelements.size(); i++) {
            if (openelements[i].GetNP() == 3) {
                cnt3++;
            }
        }
        int cnt4 = openelements.size() - cnt3;
        if (openelements.size() > 0) {
            MESHIT_LOG_WARNING(openelements.size() << " (" << cnt3 << " + " << cnt4 << ")" << " open elements");
        }

        BuildBoundaryEdges();

        for (size_t i = 0; i < openelements.size(); i++) {
            const Element2d& sel = openelements[i];

            if (boundaryedges)
                for (size_t j = 1; j <= sel.GetNP(); j++) {
                    INDEX_2 i2;
                    i2.I1() = sel.PNumMod(j);
                    i2.I2() = sel.PNumMod(j + 1);
                    i2.Sort();
                    boundaryedges->Set(i2, 1);
                }

            for (size_t j = 1; j <= 3; j++) {
                PointIndex pi = sel.PNum(j);
                if (pi < points.End()) {
                    points[pi].SetType(FIXEDPOINT);
                }
            }
        }
    }

    bool Mesh::HasOpenQuads() const
    {
        int no = GetNOpenElements();
        for (int i = 0; i < no; i++) {
            if (openelements[i].GetNP() == 4)
                return true;
        }
        return false;
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
            const Element2d& el = SurfaceElement(i+1);
            if (el.IsDeleted()) continue;

            if (surfnr == 0 || el.GetIndex() == surfnr) {
                for (size_t j = 1; j <= el.GetNP(); j++) {
                    INDEX_2 seg(el.PNumMod(j), el.PNumMod(j + 1));
                    INDEX_2 data;

                    if (seg.I1() <= 0 || seg.I2() <= 0)
                        std::cerr << "seg = " << seg << std::endl;

                    if (faceht.Used(seg)) {
                        data = faceht.Get(seg);
                        if (data.I1() == el.GetIndex()) {
                            data.I1() = 0;
                            faceht.Set(seg, data);
                        }
                        else {
                            MESHIT_LOG_WARNING("hash table si not fitting for segment: " <<
                                               seg.I1() << "-" << seg.I2() << " other = " << data.I2());
                        }
                    }
                    else {
                        std::swap(seg.I1(), seg.I2());
                        data.I1() = el.GetIndex();
                        data.I2() = i+1;

                        faceht.Set(seg, data);
                    }
                }
            }
        }

        std::cerr << "open segments: " << std::endl;
        opensegments.resize(0);
        for (size_t i = 0; i < faceht.GetNBags(); i++) {
            for (int j = 0; j < faceht.GetBagSize(i); j++) {
                INDEX_2 i2;
                INDEX_2 data;
                faceht.GetData(i, j, i2, data);
                if (data.I1()) // surfnr
                {
                    Segment seg;
                    seg[0] = i2.I1();
                    seg[1] = i2.I2();
                    seg.si = data.I1();

                    // find geomdata:
                    if (data.I2() > 0) {
                        // segment due to triangle
                        const Element2d& el = SurfaceElement(data.I2());
                        for (size_t k = 1; k <= el.GetNP(); k++) {
                            if (seg[0] == el.PNum(k))
                                seg.geominfo[0] = el.GeomInfoPi(k);
                            if (seg[1] == el.PNum(k))
                                seg.geominfo[1] = el.GeomInfoPi(k);
                        }
                        std::cerr << "trig seg: ";
                    }
                    else {
                        // segment due to line
                        const Segment& lseg = LineSegment(-(data.I2() - 1));
                        seg.geominfo[0] = lseg.geominfo[0];
                        seg.geominfo[1] = lseg.geominfo[1];
                        std::cerr << "line seg: ";
                    }

                    std::cerr << seg[0] << " - " << seg[1]
                    << " len = " << Dist(Point(seg[0]), Point(seg[1]))
                    << std::endl;

                    opensegments.push_back(seg);
                    if (seg.geominfo[0].trignum <= 0 || seg.geominfo[1].trignum <= 0) {
                        std::cerr << "Problem with open segment: " << seg << std::endl;
                    }

                }
            }
        }
        MESHIT_LOG_DEBUG(opensegments.size() << " open segments found");

        for (int i = 0; i < points.size(); i++) {
            points[i].SetType(SURFACEPOINT);
        }

        for (int i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            points[seg[0]].SetType(EDGEPOINT);
            points[seg[1]].SetType(EDGEPOINT);
        }
        for (int i = 1; i <= GetNOpenSegments(); i++) {
            const Segment& seg = GetOpenSegment(i);
            points[seg[0]].SetType(EDGEPOINT);
            points[seg[1]].SetType(EDGEPOINT);
        }
    }

    void Mesh::SetLocalH(const Point3d& pmin, const Point3d& pmax, double grading)
    {
        Point3d c = Center(pmin, pmax);
        double d = 0.5 * std::max(
                pmax.X() - pmin.X(), std::max(
                pmax.Y() - pmin.Y(),
                pmax.Z() - pmin.Z()));
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
        if (hloc < hmin)
            hloc = hmin;

        int steps = int(Dist(p1, p2) / hloc) + 2;
        Vec3d v(p1, p2);

        for (int i = 0; i <= steps; i++) {
            Point3d p = p1 + (double(i) / double(steps) * v);
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

    void Mesh::SetMaxHDomain(const Array<double>& mhd)
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
        int i, j, n;
        double hi, hsum;
        double maxh = 0, minh = 1e10;

        hsum = 0;
        n = 0;
        for (i = 1; i <= GetNSE(); i++) {
            const Element2d& el = SurfaceElement(i);
            if (surfnr == 0 || el.GetIndex() == surfnr) {
                for (j = 1; j <= 3; j++) {
                    hi = Dist(Point(el.PNumMod(j)),
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

        for (int i = 0; i < GetNSE(); i++) {
            const Element2d& el = surfelements[i];
            int j;

            if (el.GetNP() == 3) {
                double hel = -1;
                for (j = 1; j <= 3; j++) {
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
            else {
                {
                    const Point3d& p1 = points[el.PNum(1)];
                    const Point3d& p2 = points[el.PNum(2)];
                    lochfunc->SetH(Center(p1, p2), 2 * Dist(p1, p2));
                }
                {
                    const Point3d& p1 = points[el.PNum(3)];
                    const Point3d& p2 = points[el.PNum(4)];
                    lochfunc->SetH(Center(p1, p2), 2 * Dist(p1, p2));
                }
            }
        }

        for (int i = 0; i < GetNSeg(); i++) {
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

        for (PointIndex i = 0; i < GetNP(); i++) {
            for (PointIndex j = i + 1; j < GetNP(); j++) {
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
        int i, j;

        for (i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();
            bedges.Set(i2, 1);
        }
        for (i = 1; i <= GetNSE(); i++) {
            const Element2d& sel = SurfaceElement(i);
            if (!sel.PNum(1))
                continue;
            for (j = 1; j <= 3; j++) {
                INDEX_2 i2(sel.PNumMod(j), sel.PNumMod(j + 1));
                i2.Sort();
                if (bedges.Used(i2)) continue;

                if (edges.Used(i2)) {
                    int other = edges.Get(i2);

                    const Element2d& elother = SurfaceElement(other);

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
                    edges.Set(i2, i);
            }
        }


        // Restrict h due to line segments

        for (i = 0; i < GetNSeg(); i++) {
            const Segment& seg = LineSegment(i);
            const Point3d& p1 = Point(seg[0]);
            const Point3d& p2 = Point(seg[1]);
            RestrictLocalH(Center(p1, p2), Dist(p1, p2));
        }
    }

    void Mesh::RestrictLocalH(resthtype rht, int nr, double loch)
    {
        int i;
        switch (rht) {
            case RESTRICTH_FACE: {
                for (i = 1; i <= GetNSE(); i++) {
                    const Element2d& sel = SurfaceElement(i);
                    if (sel.GetIndex() == nr)
                        RestrictLocalH(RESTRICTH_SURFACEELEMENT, i, loch);
                }
                break;
            }
            case RESTRICTH_EDGE: {
                for (i = 0; i < GetNSeg(); i++) {
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
                const Element2d& sel = SurfaceElement(nr);
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
            for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
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

    void Mesh::GetBox(Point3d& pmin, Point3d& pmax, POINTTYPE ptyp) const
    {
        if (points.size() == 0) {
            pmin = pmax = Point3d(0, 0, 0);
            return;
        }

        pmin = Point3d(1e10, 1e10, 1e10);
        pmax = Point3d(-1e10, -1e10, -1e10);

        for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
            if (points[pi].Type() <= ptyp) {
                pmin.SetToMin(points[pi]);
                pmax.SetToMax(points[pi]);
            }
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
        Array<PointIndex> op2np(GetNP());
        Array<MeshPoint> hpoints;
        BitArrayChar pused(GetNP());

        for (int i = 0; i < surfelements.size(); i++) {
            if (surfelements[i].IsDeleted()) {
                surfelements.Delete(i);
                i--;
            }
        }
        for (int i = 0; i < segments.size(); i++) {
            if (segments[i][0] <= -1) {
                segments.Delete(i);
                i--;
            }
        }
        pused.Clear();

        for (int i = 0; i < surfelements.size(); i++) {
            const Element2d& el = surfelements[i];
            for (size_t j = 0; j < el.GetNP(); j++) {
                pused.Set(el[j]);
            }
        }

        for (int i = 0; i < segments.size(); i++) {
            const Segment& seg = segments[i];
            pused.Set(seg[0]);
            pused.Set(seg[1]);
        }

        for (int i = 0; i < openelements.size(); i++) {
            const Element2d& el = openelements[i];
            for (size_t j = 0; j < el.GetNP(); j++) {
                pused.Set(el[j]);
            }
        }

        for (int i = 0; i < lockedpoints.size(); i++) {
            pused.Set(lockedpoints[i]);
        }

        int npi = -1;

        for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
            if (pused.Test(pi)) {
                npi++;
                op2np[pi] = npi;
                hpoints.push_back(points[pi]);
            }
            else
                op2np[pi] = -1;
        }

        points.resize(0);
        for (int i = 0; i < hpoints.size(); i++) {
            points.push_back(hpoints[i]);
        }

        for (int i = 1; i <= surfelements.size(); i++) {
            Element2d& el = SurfaceElement(i);
            for (size_t j = 0; j < el.GetNP(); j++) {
                el[j] = op2np[el[j]];
            }
        }

        for (int i = 0; i < segments.size(); i++) {
            Segment& seg = segments[i];
            seg[0] = op2np[seg[0]];
            seg[1] = op2np[seg[1]];
        }

        for (int i = 0; i < openelements.size(); i++) {
            Element2d& el = openelements[i];
            for (size_t j = 0; j < el.GetNP(); j++) {
                el[j] = op2np[el[j]];
            }
        }

        for (int i = 0; i < lockedpoints.size(); i++) {
            lockedpoints[i] = op2np[lockedpoints[i]];
        }

        for (int i = 0; i < facedecoding.size(); i++) {
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
        int nf = GetNOpenElements();
        INDEX_2_HASHTABLE<int> edges(nf + 2);
        INDEX_2 i2, i2s, edge;
        int err = 0;

        for (int i = 1; i <= nf; i++) {
            const Element2d& sel = OpenElement(i);

            for (size_t j = 1; j <= sel.GetNP(); j++) {
                i2.I1() = sel.PNumMod(j);
                i2.I2() = sel.PNumMod(j + 1);

                int sign = (i2.I2() > i2.I1()) ? 1 : -1;
                i2.Sort();
                if (!edges.Used(i2))
                    edges.Set(i2, 0);
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
                    for (int k = 1; k <= nf; k++) {
                        const Element2d& sel = OpenElement(k);
                        for (size_t l = 1; l <= sel.GetNP(); l++) {
                            edge.I1() = sel.PNumMod(l);
                            edge.I2() = sel.PNumMod(l + 1);
                            edge.Sort();

                            if (edge == i2s)
                                MESHIT_LOG_ERROR("  edge of element " << sel);
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
        Array<size_t> inters;

        bool overlap = 0;
        bool incons_layers = 0;

        for (size_t i = 1; i <= GetNSE(); i++) {
            const Element2d& tri = SurfaceElement(i);

            Point3d tpmin(Point(tri[0]));
            Point3d tpmax(tpmin);

            for (size_t k = 1; k < tri.GetNP(); k++) {
                tpmin.SetToMin(Point(tri[k]));
                tpmax.SetToMax(Point(tri[k]));
            }
            Vec3d diag(tpmin, tpmax);

            tpmax = tpmax + 0.1 * diag;
            tpmin = tpmin - 0.1 * diag;

            setree.Insert(tpmin, tpmax, i);
        }

        for (size_t i = 1; i <= GetNSE(); i++) {
            const Element2d& tri = SurfaceElement(i);

            Point3d tpmin(Point(tri[0]));
            Point3d tpmax(tpmin);

            for (size_t k = 1; k < tri.GetNP(); k++) {
                tpmin.SetToMin(Point(tri[k]));
                tpmax.SetToMax(Point(tri[k]));
            }

            setree.GetIntersecting(tpmin, tpmax, inters);

            for (int j = 0; j < inters.size(); j++) {
                const Element2d& tri2 = SurfaceElement(inters[j]);

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
                    MESHIT_LOG_WARNING("Intersecting elements " << i << " and " << inters[j]);
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

        for (int k = 0; k < facedecoding.size(); k++) {
            if (facedecoding[k].DomainIn() > ndom)
                ndom = facedecoding[k].DomainIn();
            if (facedecoding[k].DomainOut() > ndom)
                ndom = facedecoding[k].DomainOut();
        }

        return ndom;
    }

    void Mesh::BuildElementSearchTree()
    {
        if (elementsearchtreets == GetTimeStamp()) return;

        //#pragma omp critical (buildsearchtree)
        {
            if (elementsearchtreets != GetTimeStamp()) {

                MESHIT_LOG_DEBUG("Rebuild element searchtree");

                delete elementsearchtree;
                elementsearchtree = NULL;

                int ne = GetNSE();

                if (ne) {
                    if (dimension == 2) {
                        Box<3> box(Box<3>::EMPTY_BOX);
                        for (size_t sei = 0; sei < ne; sei++) {
                            box.Add(points[surfelements[sei].PNums()]);
                        }
                        box.Increase(1.01 * box.Diam());
                        elementsearchtree = new Box3dTree(box);

                        for (size_t sei = 0; sei < ne; sei++) {
                            box.Set(points[surfelements[sei].PNums()]);
                            elementsearchtree->Insert(box, sei + 1);
                        }
                    }

                    elementsearchtreets = GetTimeStamp();
                }
            }
        }
    }

    bool Mesh::PointContainedIn2DElement(const Point3d& p,
                                         double lami[3],
                                         const int element,
                                         bool consider3D) const
    {
        Vec3d col1, col2, col3;
        Vec3d rhs, sol;
        const double eps = 1e-6;

        Array<Element2d> loctrigs;

        //SZ 
        if (SurfaceElement(element).GetType() == QUAD) {
            const Element2d& el = SurfaceElement(element);

            const Point3d& p1 = Point(el.PNum(1));
            const Point3d& p2 = Point(el.PNum(2));
            const Point3d& p3 = Point(el.PNum(3));
            const Point3d& p4 = Point(el.PNum(4));

            // Coefficients of Bilinear Mapping from Ref-Elem to global Elem
            // X = a + b x + c y + d x y 
            Vec3d a = p1;
            Vec3d b = p2 - a;
            Vec3d c = p4 - a;
            Vec3d d = p3 - a - b - c;

            double dxb = d.X() * b.Y() - d.Y() * b.X();
            double dxc = d.X() * c.Y() - d.Y() * c.X();
            double dxa = d.X() * a.Y() - d.Y() * a.X();
            double dxp = d.X() * p.Y() - d.Y() * p.X();

            const double eps = 1.E-12;
            double c0, c1, c2;

            lami[2] = 0.;

            if (fabs(d.X()) <= eps && fabs(d.Y()) <= eps) {
                // Solve Linear System
                lami[0] = (c.Y() * (p.X() - a.X()) - c.X() * (p.Y() - a.Y())) / (b.X() * c.Y() - b.Y() * c.X());
                lami[1] = (-b.Y() * (p.X() - a.X()) + b.X() * (p.Y() - a.Y())) / (b.X() * c.Y() - b.Y() * c.X());
            } else if (fabs(dxb) <= eps) {
                lami[1] = (dxp - dxa) / dxc;
                if (fabs(b.X() - d.X() * lami[1]) >= eps) {
                    lami[0] = (p.X() - a.X() - c.X() * lami[1]) / (b.X() + d.X() * lami[1]);
                } else {
                    lami[0] = (p.Y() - a.Y() - c.Y() * lami[1]) / (b.Y() + d.Y() * lami[1]);
                }
            } else if (fabs(dxc) <= eps) {
                lami[0] = (dxp - dxa) / dxb;
                if (fabs(c.X() - d.X() * lami[0]) >= eps) {
                    lami[1] = (p.X() - a.X() - b.X() * lami[0]) / (c.X() + d.X() * lami[0]);
                } else {
                    lami[1] = (p.Y() - a.Y() - b.Y() * lami[0]) / (c.Y() + d.Y() * lami[0]);
                }
            } else {
                // Solve quadratic equation
                if (fabs(d.X()) >= eps) {
                    c2 = d.X() * dxc;
                    c1 = d.X() * dxc - c.X() * dxb - d.X() * (dxp - dxa);
                    c0 = -b.X() * (dxp - dxa) - (a.X() - p.X()) * dxb;
                }
                else {
                    c2 = d.Y() * dxc;
                    c1 = d.Y() * dxc - c.Y() * dxb - d.Y() * (dxp - dxa);
                    c0 = -b.Y() * (dxp - dxa) - (a.Y() - p.Y()) * dxb;
                }

                double rt = c1 * c1 - 4 * c2 * c0;
                if (rt < 0.) return false;
                lami[1] = (-c1 + sqrt(rt)) / 2 / c2;
                if (lami[1] <= 1. && lami[1] >= 0.) {
                    lami[0] = (dxp - dxa - dxc * lami[1]) / dxb;
                    if (lami[0] <= 1. && lami[0] >= 0.)
                        return true;
                }

                lami[1] = (-c1 - sqrt(rt)) / 2 / c2;
                lami[0] = (dxp - dxa - dxc * lami[1]) / dxb;
            }

            if (lami[0] <= 1. + eps && lami[0] >= -eps && lami[1] <= 1. + eps && lami[1] >= -eps) {
                if (consider3D) {
                    Vec3d n = Cross(b, c);
                    lami[2] = 0;
                    for (int i = 1; i <= 3; i++) {
                        lami[2] += (p.X(i) - a.X(i) - lami[0] * b.X(i) - lami[1] * c.X(i)) * n.X(i);
                    }
                    if (lami[2] >= -eps && lami[2] <= eps)
                        return true;
                }
                else
                    return true;
            }

            return false;

        } else {
            loctrigs.resize(1);
            loctrigs[0] = SurfaceElement(element);

            for (int j = 0; j < loctrigs.size(); j++) {
                const Element2d& el = loctrigs[j];


                const Point3d& p1 = Point(el.PNum(1));
                const Point3d& p2 = Point(el.PNum(2));
                const Point3d& p3 = Point(el.PNum(3));

                col1 = p2 - p1;
                col2 = p3 - p1;
                col3 = Cross(col1, col2);
                rhs = p - p1;

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
            }
        }
        return false;
    }

    void Mesh::RebuildSurfaceElementLists()
    {
        for (int i = 0; i < facedecoding.size(); i++) {
            facedecoding[i].firstelement = -1;
        }
        for (int i = surfelements.size() - 1; i >= 0; i--) {
            int ind = surfelements[i].GetIndex();
            surfelements[i].next = facedecoding[ind - 1].firstelement;
            facedecoding[ind - 1].firstelement = i;
        }
    }

    void Mesh::GetSurfaceElementsOfFace(int facenr, Array<SurfaceElementIndex>& sei) const
    {
        /* Philippose - 01/10/2009
        Commented out the following lines, and activated the originally 
        commented out lines above because of a bug which causes corruption 
        of the variable "facedecoding" when a mesh is converted to second order
         */
        sei.resize(0);

        SurfaceElementIndex si = facedecoding[facenr - 1].firstelement;
        while (si != -1) {
            const Element2d& se = SurfaceElement(si);
            if (se.GetIndex() == facenr
                && se[0] >= 0
                && !se.IsDeleted()) {
                sei.push_back(si);
            }
            si = se.next;
        }
    }

    void Mesh::InitPointCurve(double red, double green, double blue) const
    {
        pointcurves_startpoint.push_back(pointcurves.size());
        pointcurves_red.push_back(red);
        pointcurves_green.push_back(green);
        pointcurves_blue.push_back(blue);
    }

    void Mesh::AddPointCurvePoint(const Point3d& pt) const
    {
        pointcurves.push_back(pt);
    }

    int Mesh::GetNumPointCurves(void) const
    {
        return pointcurves_startpoint.size();
    }

    int Mesh::GetNumPointsOfPointCurve(int curve) const
    {
        if (curve == pointcurves_startpoint.size() - 1)
            return (pointcurves.size() - pointcurves_startpoint.Last());
        else
            return (pointcurves_startpoint[curve + 1] - pointcurves_startpoint[curve]);
    }

    Point3d& Mesh::GetPointCurvePoint(int curve, int n) const
    {
        return pointcurves[pointcurves_startpoint[curve] + n];
    }

    void Mesh::GetPointCurveColor(int curve, double& red, double& green, double& blue) const
    {
        red = pointcurves_red[curve];
        green = pointcurves_green[curve];
        blue = pointcurves_blue[curve];
    }

    void Mesh::ComputeNVertices()
    {
        int i, j, nv;
        int nse = GetNSE();

        numvertices = 0;
        for (i = 1; i <= nse; i++) {
            const Element2d& el = SurfaceElement(i);
            nv = el.GetNV();
            for (j = 1; j <= nv; j++) {
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

    bool Mesh::PureTrigMesh(int faceindex) const
    {
        if (!faceindex) {
            for (int i = 1; i <= GetNSE(); i++) {
                if (SurfaceElement(i).GetNP() != 3)
                    return false;
            }
            return true;
        }

        for (int i = 1; i <= GetNSE(); i++) {
            if (SurfaceElement(i).GetIndex() == faceindex && SurfaceElement(i).GetNP() != 3)
                return false;
        }
        return true;
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

    void Mesh::SetNBCNames(int nbcn)
    {
        for (size_t i = 0; i < bcnames.size(); i++) {
            if (bcnames[i]) {
                delete bcnames[i];
            }
        }
        bcnames.resize(nbcn);
        bcnames = 0;
    }

    void Mesh::SetBCName(int bcnr, const std::string& abcname)
    {
        if (bcnames[bcnr]) delete bcnames[bcnr];
        if (abcname != "default")
            bcnames[bcnr] = new std::string(abcname);
        else
            bcnames[bcnr] = 0;
    }

    const std::string& Mesh::GetBCName(int bcnr) const
    {
        static std::string defaultstring = "default";

        if (!bcnames.size())
            return defaultstring;
        if (bcnames[bcnr])
            return *bcnames[bcnr];
        else
            return defaultstring;
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
