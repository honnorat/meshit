#include <meshit.hpp>

#include <stdexcept>

#include "curvedelems.hpp"
#include "clusters.hpp"
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
        mglevels = 1;
        elementsearchtree = NULL;
        elementsearchtreets = NextTimeStamp();
        majortimestamp = timestamp = NextTimeStamp();
        hglob = 1e10;
        hmin = 0;
        numvertices = -1;
        dimension = 3;

        topology = new MeshTopology(*this);
        curvedelems = new CurvedElements(*this);
        clusters = new AnisotropicClusters(*this);
        ident = new Identifications(*this);

        hpelements = NULL;
        coarsemesh = NULL;

        ps_startelement = 0;

        geomtype = NO_GEOM;

        bcnames.resize(0);
    }

    Mesh::~Mesh()
    {
        delete lochfunc;
        delete boundaryedges;
        delete surfelementht;
        delete segmentht;
        delete curvedelems;
        delete clusters;
        delete topology;
        delete ident;
        delete elementsearchtree;
        delete coarsemesh;
        delete hpelements;

        for (int i = 0; i < materials.size(); i++) {
            delete [] materials[i];
        }
        for (int i = 0; i < userdata_int.Size(); i++) {
            delete userdata_int[i];
        }
        for (int i = 0; i < userdata_double.Size(); i++) {
            delete userdata_double[i];
        }
        for (int i = 0; i < bcnames.size(); i++) {
            if (bcnames[i]) delete bcnames[i];
        }
    }

    Mesh & Mesh::operator=(const Mesh & mesh2)
    {
        points = mesh2.points;
        // eltyps = mesh2.eltyps;
        segments = mesh2.segments;
        surfelements = mesh2.surfelements;
        volelements = mesh2.volelements;
        lockedpoints = mesh2.lockedpoints;
        facedecoding = mesh2.facedecoding;
        dimension = mesh2.dimension;

        bcnames.resize(mesh2.bcnames.size());
        for (int i = 0; i < mesh2.bcnames.size(); i++) {
            if (mesh2.bcnames[i])
                bcnames[i] = new std::string(*mesh2.bcnames[i]);
            else
                bcnames[i] = 0;
        }

        return *this;
    }

    void Mesh::DeleteMesh()
    {
        points.resize(0);
        segments.resize(0);
        surfelements.resize(0);
        volelements.resize(0);
        lockedpoints.resize(0);
        surfacesonnode.SetSize(0);

        delete boundaryedges;
        boundaryedges = NULL;

        openelements.resize(0);
        facedecoding.resize(0);

        delete ident;
        ident = new Identifications(*this);
        delete topology;
        topology = new MeshTopology(*this);
        delete curvedelems;
        curvedelems = new CurvedElements(*this);
        delete clusters;
        clusters = new AnisotropicClusters(*this);

        for (int i = 0; i < bcnames.size(); i++) {
            if (bcnames[i]) delete bcnames[i];
        }

        timestamp = NextTimeStamp();
    }

    void Mesh::ClearSurfaceElements()
    {
        surfelements.resize(0);
        for (int i = 0; i < facedecoding.size(); i++) {
            facedecoding[i].firstelement = -1;
        }

        timestamp = NextTimeStamp();
    }

    void Mesh::BuildFromSpline2D(SplineGeometry2d & geometry, MeshingParameters & mp)
    {
        LOG_DEBUG("Generate Mesh from spline geometry");

        double h = mp.maxh;

        Box<2> bbox = geometry.GetBoundingBox();

        if (bbox.Diam() < h) {
            h = bbox.Diam();
            mp.maxh = h;
        }

        dimension = 2;

        Point3d pmin(bbox.PMin()(0), bbox.PMin()(1), -bbox.Diam());
        Point3d pmax(bbox.PMax()(0), bbox.PMax()(1), bbox.Diam());

        SetLocalH(pmin, pmax, mp.grading);
        SetGlobalH(h);

        geometry.PartitionBoundary(mp, h, *this);

        // marks mesh points for hp-refinement
        for (int i = 0; i < geometry.GetNP(); i++) {
            if (geometry.GetPoint(i).hpref) {
                double mindist = 1e99;
                PointIndex mpi(0);
                ::meshit::Point<2> gp = geometry.GetPoint(i);
                ::meshit::Point<3> gp3(gp(0), gp(1), 0);
                for (PointIndex pi = PointIndex::BASE; pi < GetNP() + PointIndex::BASE; pi++)
                    if (Dist2(gp3, points[pi]) < mindist) {
                        mpi = pi;
                        mindist = Dist2(gp3, points[pi]);
                    }
                points[mpi].Singularity(1.);
            }
        }

        int maxdomnr = 0;
        for (SegmentIndex si = 0; si < GetNSeg(); si++) {
            if (segments[si].domin > maxdomnr) maxdomnr = segments[si].domin;
            if (segments[si].domout > maxdomnr) maxdomnr = segments[si].domout;
        }

        ClearFaceDescriptors();
        for (int i = 1; i <= maxdomnr; i++) {
            AddFaceDescriptor(FaceDescriptor(i, 0, 0, i));
        }

        int maxsegmentindex = 0;
        for (SegmentIndex si = 0; si < GetNSeg(); si++) {
            if (segments[si].si > maxsegmentindex) maxsegmentindex = segments[si].si;
        }

        SetNBCNames(maxsegmentindex);

        for (int sindex = 0; sindex < maxsegmentindex; sindex++) {
            SetBCName(sindex, geometry.GetBCName(sindex + 1));
        }

        for (SegmentIndex si = 0; si < GetNSeg(); si++) {
            segments[si].SetBCName(bcnames[segments[si].si - 1]);
        }

        CalcLocalH(mp.grading);

        int bnp = GetNP(); // boundary points
        int hquad = mp.quad;

        for (int domnr = 1; domnr <= maxdomnr; domnr++) {

            if (!geometry.GetDomainTensorMeshing(domnr))
                continue;

            // tensor product mesh

            Array<PointIndex, PointIndex::BASE> nextpi(bnp);
            Array<int, PointIndex::BASE> si1(bnp), si2(bnp);
            PointIndex firstpi;

            nextpi = -1;
            si1 = -1;
            si2 = -1;
            for (SegmentIndex si = 0; si < GetNSeg(); si++) {
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

                nextpi[p1] = p2; // counter-clockwise

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

            PointIndex c1(0), c2, c3, c4; // 4 corner points
            int nex = 1, ney = 1;

            for (PointIndex pi = 1; pi <= si2.size(); pi++) {
                if (si2[pi] != -1) {
                    c1 = pi;
                    break;
                }
            }
            for (c2 = nextpi[c1]; si2[c2] == -1; c2 = nextpi[c2], nex++);
            for (c3 = nextpi[c2]; si2[c3] == -1; c3 = nextpi[c3], ney++);
            for (c4 = nextpi[c3]; si2[c4] == -1; c4 = nextpi[c4]);

            Array<PointIndex> pts((nex + 1) * (ney + 1)); // x ... inner loop
            pts = -1;

            for (PointIndex pi = c1, i = 0; pi != c2; pi = nextpi[pi], i++) {
                pts[i] = pi;
            }
            for (PointIndex pi = c2, i = 0; pi != c3; pi = nextpi[pi], i++) {
                pts[(nex + 1) * i + nex] = pi;
            }
            for (PointIndex pi = c3, i = 0; pi != c4; pi = nextpi[pi], i++) {
                pts[(nex + 1)*(ney + 1) - i - 1] = pi;
            }
            for (PointIndex pi = c4, i = 0; pi != c1; pi = nextpi[pi], i++) {
                pts[(nex + 1)*(ney - i)] = pi;
            }

            for (PointIndex pix = nextpi[c1], ix = 0; pix != c2; pix = nextpi[pix], ix++) {
                for (PointIndex piy = nextpi[c2], iy = 0; piy != c3; piy = nextpi[piy], iy++) {
                    ::meshit::Point<3> p = points[pix] + (points[piy] - points[c2]);
                    pts[(nex + 1)*(iy + 1) + ix + 1] = AddPoint(p, 1, FIXEDPOINT);
                }
            }

            for (int i = 0; i < ney; i++) {
                for (int j = 0; j < nex; j++) {
                    Element2d el(QUAD);
                    el[0] = pts[i * (nex + 1) + j];
                    el[1] = pts[i * (nex + 1) + j + 1];
                    el[2] = pts[(i + 1)*(nex + 1) + j + 1];
                    el[3] = pts[(i + 1)*(nex + 1) + j];
                    el.SetIndex(domnr);

                    AddSurfaceElement(el);
                }
            }
        }

        for (int domnr = 1; domnr <= maxdomnr; domnr++) {

            if (geometry.GetDomainTensorMeshing(domnr)) continue;

            if (geometry.GetDomainMaxh(domnr) > 0)
                h = geometry.GetDomainMaxh(domnr);

            LOG_DEBUG("Meshing domain " << domnr << " / " << maxdomnr);

            int oldnf = GetNSE();

            mp.quad = hquad || geometry.GetDomainQuadMeshing(domnr);

            Meshing2 meshing(mp, Box<3> (pmin, pmax));

            Array<int, PointIndex::BASE> compress(bnp);
            compress = -1;
            int cnt = 0;
            for (PointIndex pi = PointIndex::BASE; pi < bnp + PointIndex::BASE; pi++) {
                if (points[pi].GetLayer() == geometry.GetDomainLayer(domnr)) {
                    meshing.AddPoint(points[pi], pi);
                    cnt++;
                    compress[pi] = cnt;
                }
            }
            PointGeomInfo gi;
            gi.trignum = 1;
            for (SegmentIndex si = 0; si < GetNSeg(); si++) {
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

            for (SurfaceElementIndex sei = oldnf; sei < GetNSE(); sei++) {
                surfelements[sei].SetIndex(domnr);
            }

            // astrid
            char * material;
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
        SetNextMajorTimeStamp();
    }

    PointIndex Mesh::AddPoint(const Point3d & p, int layer)
    {
        return AddPoint(p, layer, INNERPOINT);
    }

    PointIndex Mesh::AddPoint(const Point3d & p, int layer, POINTTYPE type)
    {
        timestamp = NextTimeStamp();

        PointIndex pi = points.End();
        points.push_back(MeshPoint(p, layer, type));

        return pi;
    }

    SegmentIndex Mesh::AddSegment(const Segment & s)
    {
        timestamp = NextTimeStamp();

        int maxn = max2(s[0], s[1]);
        maxn += 1 - PointIndex::BASE;

        if (maxn <= points.size()) {
            if (points[s[0]].Type() > EDGEPOINT)
                points[s[0]].SetType(EDGEPOINT);
            if (points[s[1]].Type() > EDGEPOINT)
                points[s[1]].SetType(EDGEPOINT);
        }

        SegmentIndex si = segments.size();
        segments.push_back(s);

        return si;
    }

    SurfaceElementIndex Mesh::AddSurfaceElement(const Element2d & el)
    {
        timestamp = NextTimeStamp();

        int maxn = el[0];
        for (int i = 1; i < el.GetNP(); i++) {
            if (el[i] > maxn) maxn = el[i];
        }
        maxn += 1 - PointIndex::BASE;

        if (maxn <= points.size()) {
            for (int i = 0; i < el.GetNP(); i++) {
                if (points[el[i]].Type() > SURFACEPOINT)
                    points[el[i]].SetType(SURFACEPOINT);
            }
        }

        SurfaceElementIndex si = surfelements.size();
        surfelements.push_back(el);

        if (el.index > facedecoding.size())
            std::cerr << "has no facedecoding: fd.size = " << facedecoding.size() << ", ind = " << el.index << std::endl;

        surfelements.Last().next = facedecoding[el.index - 1].firstelement;
        facedecoding[el.index - 1].firstelement = si;

        if (SurfaceArea().Valid())
            SurfaceArea().Add(el);

        return si;
    }

    ElementIndex Mesh::AddVolumeElement(const Element & el)
    {
        int maxn = el[0];
        for (int i = 1; i < el.GetNP(); i++) {
            if (el[i] > maxn) maxn = el[i];
        }

        maxn += 1 - PointIndex::BASE;

        int ve = volelements.size();

        volelements.push_back(el);
        volelements.Last().flags.illegal_valid = 0;

        timestamp = NextTimeStamp();

        return ve;
    }

    void Mesh::Export(
            const std::string & filename,
            const std::string & filetype) const
    {
        WriteUserFormat(filetype, const_cast<const Mesh &> (*this), filename);
    }

    void Mesh::Save(const std::string & filename) const
    {

        std::ofstream outfile(filename.c_str());
        // ogzstream outfile( (filename+".gz") .c_str());

        Save(outfile);
    }

    void Mesh::Save(std::ostream & outfile) const
    {
        int i, j;

        double scale = 1; // globflags.GetNumFlag ("scale", 1);
        int inverttets = 0; // globflags.GetDefineFlag ("inverttets");
        int invertsurf = 0; // globflags.GetDefineFlag ("invertsurfacemesh");



        outfile << "mesh3d" << "\n";

        outfile << "dimension\n" << GetDimension() << "\n";

        outfile << "geomtype\n" << int(geomtype) << "\n";


        outfile << "\n";
        outfile << "# surfnr    bcnr   domin  domout      np      p1      p2      p3"
                << "\n";


        switch (geomtype) {
            case GEOM_STL:
                outfile << "surfaceelementsgi" << "\n";
                break;
            case GEOM_OCC: case GEOM_ACIS:
                outfile << "surfaceelementsuv" << "\n";
                break;
            default:
                outfile << "surfaceelements" << "\n";
        }

        outfile << GetNSE() << "\n";

        SurfaceElementIndex sei;
        for (sei = 0; sei < GetNSE(); sei++) {
            if (surfelements[sei].GetIndex()) {
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).SurfNr() + 1;
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).BCProperty();
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).DomainIn();
                outfile << " " << GetFaceDescriptor(surfelements[sei].GetIndex()).DomainOut();
            }
            else
                outfile << " 0 0 0";

            Element2d sel = surfelements[sei];
            if (invertsurf)
                sel.Invert();

            outfile << " " << sel.GetNP();
            for (j = 0; j < sel.GetNP(); j++) {
                outfile << " " << sel[j];
            }

            switch (geomtype) {
                case GEOM_STL:
                    for (j = 1; j <= sel.GetNP(); j++) {
                        outfile << " " << sel.GeomInfoPi(j).trignum;
                    }
                    break;
                case GEOM_OCC: case GEOM_ACIS:
                    for (j = 1; j <= sel.GetNP(); j++) {
                        outfile << " " << sel.GeomInfoPi(j).u;
                        outfile << " " << sel.GeomInfoPi(j).v;
                    }
                    break;
                default:
                    ;
            }
            outfile << "\n";
        }

        outfile << "\n" << "\n";
        outfile << "#  matnr      np      p1      p2      p3      p4" << "\n";
        outfile << "volumeelements" << "\n";
        outfile << GetNE() << "\n";

        for (ElementIndex ei = 0; ei < GetNE(); ei++) {
            outfile << volelements[ei].GetIndex();
            outfile << " " << volelements[ei].GetNP();

            Element el = volelements[ei];
            if (inverttets) el.Invert();

            for (j = 0; j < el.GetNP(); j++) {
                outfile << " " << el[j];
            }
            outfile << "\n";
        }

        outfile << "\n" << "\n";
        outfile << "# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2   ednr1   dist1   ednr2   dist2 \n";
        outfile << "edgesegmentsgi2" << "\n";
        outfile << GetNSeg() << "\n";

        for (i = 1; i <= GetNSeg(); i++) {
            const Segment & seg = LineSegment(i);
            outfile.width(8);
            outfile << seg.si; // 2D: bc number, 3D: wievielte Kante
            outfile.width(8);
            outfile << 0;
            outfile.width(8);
            outfile << seg[0];
            outfile.width(8);
            outfile << seg[1];
            outfile << " ";
            outfile.width(8);
            outfile << seg.geominfo[0].trignum; // stl dreiecke
            outfile << " ";
            outfile.width(8);
            outfile << seg.geominfo[1].trignum; // <<std::endl;  // stl dreieck

            if (dimension == 3) {
                outfile << " ";
                outfile.width(8);
                outfile << seg.surfnr1 + 1;
                outfile << " ";
                outfile.width(8);
                outfile << seg.surfnr2 + 1;
            }
            else {
                outfile << " ";
                outfile.width(8);
                outfile << seg.domin;
                outfile << " ";
                outfile.width(8);
                outfile << seg.domout;
            }

            outfile << " ";
            outfile.width(8);
            outfile << seg.edgenr;
            outfile << " ";
            outfile.width(12);
            outfile.precision(16);
            outfile << seg.epgeominfo[0].dist; // splineparameter (2D)
            outfile << " ";
            outfile.width(8);
            outfile.precision(16);
            outfile << seg.epgeominfo[1].edgenr; // geometry dependent
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

        for (PointIndex pi = PointIndex::BASE; pi < GetNP() + PointIndex::BASE; pi++) {
            outfile.width(22);
            outfile << points[pi](0) / scale << "  ";
            outfile.width(22);
            outfile << points[pi](1) / scale << "  ";
            outfile.width(22);
            outfile << points[pi](2) / scale << "\n";
        }

        if (ident -> GetMaxNr() > 0) {
            outfile << "identifications\n";
            Array<INDEX_2> identpairs;
            int cnt = 0;
            for (i = 1; i <= ident -> GetMaxNr(); i++) {
                ident -> GetPairs(i, identpairs);
                cnt += identpairs.size();
            }
            outfile << cnt << "\n";
            for (i = 1; i <= ident -> GetMaxNr(); i++) {
                ident -> GetPairs(i, identpairs);
                for (j = 1; j <= identpairs.size(); j++) {
                    outfile.width(8);
                    outfile << identpairs.Get(j).I1();
                    outfile.width(8);
                    outfile << identpairs.Get(j).I2();
                    outfile.width(8);
                    outfile << i << "\n";
                }
            }

            outfile << "identificationtypes\n";
            outfile << ident -> GetMaxNr() << "\n";
            for (i = 1; i <= ident -> GetMaxNr(); i++) {
                int type = ident -> GetType(i);
                outfile << " " << type;
            }
            outfile << "\n";
        }

        int cntmat = 0;
        for (i = 1; i <= materials.size(); i++) {
            if (materials.Get(i) && strlen(materials.Get(i)))
                cntmat++;
        }

        if (cntmat) {
            outfile << "materials" << std::endl;
            outfile << cntmat << std::endl;
            for (i = 1; i <= materials.size(); i++) {
                if (materials.Get(i) && strlen(materials.Get(i)))
                    outfile << i << " " << materials.Get(i) << std::endl;
            }
        }

        int cntbcnames = 0;
        for (int ii = 0; ii < bcnames.size(); ii++) {
            if (bcnames[ii]) cntbcnames++;
        }

        if (cntbcnames) {
            outfile << "\n\nbcnames" << std::endl << bcnames.size() << std::endl;
            for (i = 0; i < bcnames.size(); i++) {
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
        for (SegmentIndex si = 0; si < GetNSeg(); si++) {
            if (segments[si].singedge_left) cnt_sing++;
        }
        if (cnt_sing) {
            outfile << "singular_edge_left" << std::endl << cnt_sing << std::endl;
            for (SegmentIndex si = 0; si < GetNSeg(); si++) {
                if (segments[si].singedge_left) {
                    outfile << int(si) << "\t" << segments[si].singedge_left << std::endl;
                }
            }
        }

        cnt_sing = 0;
        for (SegmentIndex si = 0; si < GetNSeg(); si++) {
            if (segments[si].singedge_right) cnt_sing++;
        }
        if (cnt_sing) {
            outfile << "singular_edge_right" << std::endl << cnt_sing << std::endl;
            for (SegmentIndex si = 0; si < GetNSeg(); si++) {
                if (segments[si].singedge_right)
                    outfile << int(si) << "\t" << segments[si].singedge_right << std::endl;
            }
        }


        cnt_sing = 0;
        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
            if (GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular)
                cnt_sing++;
        }

        if (cnt_sing) {
            outfile << "singular_face_inside" << std::endl << cnt_sing << std::endl;
            for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
                if (GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular) {
                    outfile << int(sei) << "\t"
                            << GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular << std::endl;
                }
            }
        }

        cnt_sing = 0;
        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
            if (GetFaceDescriptor(surfelements[sei].GetIndex()).domout_singular) cnt_sing++;
        }
        if (cnt_sing) {
            outfile << "singular_face_outside" << std::endl << cnt_sing << std::endl;
            for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
                if (GetFaceDescriptor(surfelements[sei].GetIndex()).domout_singular) {
                    outfile << int(sei) << "\t"
                            << GetFaceDescriptor(surfelements[sei].GetIndex()).domout_singular << std::endl;
                }
            }
        }


        // Philippose - 09/07/2009
        // Add mesh face colours to Netgen Vol file format
        // The colours are saved in RGB triplets
        int cnt_facedesc = GetNFD();
        if (cnt_facedesc) {
            outfile << std::endl << std::endl << "#   Surfnr     Red     Green     Blue" << std::endl;
            outfile << "face_colours" << std::endl << cnt_facedesc << std::endl;

            outfile.precision(8);
            outfile.setf(std::ios::fixed, std::ios::floatfield);
            outfile.setf(std::ios::showpoint);

            for (i = 1; i <= cnt_facedesc; i++) {
                outfile.width(8);
                outfile << GetFaceDescriptor(i).SurfNr() + 1 << " ";
                outfile.width(12);
                outfile << GetFaceDescriptor(i).SurfColour().X() << " ";
                outfile.width(12);
                outfile << GetFaceDescriptor(i).SurfColour().Y() << " ";
                outfile.width(12);
                outfile << GetFaceDescriptor(i).SurfColour().Z();
                outfile << std::endl;
            }
        }

    }

    void Mesh::Load(const std::string & filename)
    {

        std::ifstream infile(filename.c_str());
        if (!infile.good())
            throw std::runtime_error("mesh file not found");

        Load(infile);
    }

    void Mesh::Load(std::istream & infile)
    {

        char str[100];
        int i, n;

        double scale = 1; // globflags.GetNumFlag ("scale", 1);
        int inverttets = 0; // globflags.GetDefineFlag ("inverttets");
        int invertsurf = 0; // globflags.GetDefineFlag ("invertsurfacemesh");


        facedecoding.resize(0);

        bool endmesh = false;

        while (infile.good() && !endmesh) {
            infile >> str;

            if (strcmp(str, "dimension") == 0) {
                infile >> dimension;
            }

            if (strcmp(str, "geomtype") == 0) {
                int hi;
                infile >> hi;
                geomtype = GEOM_TYPE(hi);
            }

            if (strcmp(str, "surfaceelements") == 0 || strcmp(str, "surfaceelementsgi") == 0 || strcmp(str, "surfaceelementsuv") == 0) {
                infile >> n;
                LOG_DEBUG(n << " surface elements");

                bool geominfo = strcmp(str, "surfaceelementsgi") == 0;
                bool uv = strcmp(str, "surfaceelementsuv") == 0;


                for (i = 1; i <= n; i++) {
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

            if (strcmp(str, "volumeelements") == 0) {
                infile >> n;
                LOG_DEBUG(n << " volume elements");
                for (i = 1; i <= n; i++) {
                    Element el;
                    int hi, nep;
                    infile >> hi;
                    if (hi == 0) hi = 1;
                    el.SetIndex(hi);
                    infile >> nep;
                    el.SetNP(nep);

                    for (int j = 0; j < nep; j++) {
                        infile >> (int&) (el[j]);
                    }
                    if (inverttets) {
                        el.Invert();
                    }
                    AddVolumeElement(el);
                }
            }


            if (strcmp(str, "edgesegments") == 0) {
                infile >> n;
                for (i = 1; i <= n; i++) {
                    Segment seg;
                    int hi;
                    infile >> seg.si >> hi >> seg[0] >> seg[1];
                    AddSegment(seg);
                }
            }



            if (strcmp(str, "edgesegmentsgi") == 0) {
                infile >> n;
                for (i = 1; i <= n; i++) {
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

                LOG_DEBUG(n << " curve elements");

                for (i = 1; i <= n; i++) {
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
                LOG_DEBUG(n << " points");
                for (i = 1; i <= n; i++) {
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
                LOG_DEBUG(n << " identifications");
                for (i = 1; i <= n; i++) {
                    PointIndex pi1, pi2;
                    int ind;
                    infile >> pi1 >> pi2 >> ind;
                    ident -> Add(pi1, pi2, ind);
                }
            }

            if (strcmp(str, "identificationtypes") == 0) {
                infile >> n;
                LOG_DEBUG(n << " identificationtypes");
                for (i = 1; i <= n; i++) {
                    int type;
                    infile >> type;
                    ident -> SetType(i, Identifications::ID_TYPE(type));
                }
            }

            if (strcmp(str, "materials") == 0) {
                infile >> n;
                LOG_DEBUG(n << " materials");
                for (i = 1; i <= n; i++) {
                    int nr;
                    std::string mat;
                    infile >> nr >> mat;
                    SetMaterial(nr, mat.c_str());
                }
            }

            if (strcmp(str, "bcnames") == 0) {
                infile >> n;
                LOG_DEBUG(n << " bcnames");
                Array<int, 0> bcnrs(n);
                SetNBCNames(n);
                for (i = 1; i <= n; i++) {
                    std::string nextbcname;
                    infile >> bcnrs[i - 1] >> nextbcname;
                    bcnames[bcnrs[i - 1] - 1] = new std::string(nextbcname);
                }

                if (GetDimension() == 2) {
                    for (i = 1; i <= GetNSeg(); i++) {
                        Segment & seg = LineSegment(i);
                        if (seg.si <= n)
                            seg.SetBCName(bcnames[seg.si - 1]);
                        else
                            seg.SetBCName(0);
                    }
                }
                else {
                    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
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
                for (i = 1; i <= n; i++) {
                    PointIndex pi;
                    double s;
                    infile >> pi;
                    infile >> s;
                    points[pi].Singularity(s);
                }
            }

            if (strcmp(str, "singular_edge_left") == 0) {
                infile >> n;
                for (i = 1; i <= n; i++) {
                    SegmentIndex si;
                    double s;
                    infile >> si;
                    infile >> s;
                    segments[si].singedge_left = s;
                }
            }
            if (strcmp(str, "singular_edge_right") == 0) {
                infile >> n;
                for (i = 1; i <= n; i++) {
                    SegmentIndex si;
                    double s;
                    infile >> si;
                    infile >> s;
                    segments[si].singedge_right = s;
                }
            }

            if (strcmp(str, "singular_face_inside") == 0) {
                infile >> n;
                for (i = 1; i <= n; i++) {
                    SurfaceElementIndex sei;
                    double s;
                    infile >> sei;
                    infile >> s;
                    GetFaceDescriptor(surfelements[sei].GetIndex()).domin_singular = s;
                }
            }

            if (strcmp(str, "singular_face_outside") == 0) {
                infile >> n;
                for (i = 1; i <= n; i++) {
                    SurfaceElementIndex sei;
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

        topology -> Update();
        clusters -> Update();

        SetNextMajorTimeStamp();
    }

    void Mesh::BuildBoundaryEdges(void)
    {
        delete boundaryedges;

        boundaryedges = new INDEX_2_CLOSED_HASHTABLE<int>
                (3 * (GetNSE() + GetNOpenElements()) + GetNSeg() + 1);


        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
            const Element2d & sel = surfelements[sei];
            if (sel.IsDeleted()) continue;

            // int si = sel.GetIndex();

            if (sel.GetNP() <= 4)
                for (int j = 0; j < sel.GetNP(); j++) {
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
            const Element2d & sel = openelements[i];
            for (int j = 0; j < sel.GetNP(); j++) {
                INDEX_2 i2;
                i2.I1() = sel.PNumMod(j + 1);
                i2.I2() = sel.PNumMod(j + 2);
                i2.Sort();
                boundaryedges->Set(i2, 1);

                points[sel[j]].SetType(FIXEDPOINT);
            }
        }

        for (int i = 0; i < GetNSeg(); i++) {
            const Segment & seg = segments[i];
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();

            boundaryedges -> Set(i2, 2);
            //segmentht -> Set (i2, i);
        }


    }

    void Mesh::CalcSurfacesOfNode()
    {
        surfacesonnode.SetSize(GetNP());

        delete boundaryedges;
        boundaryedges = NULL;

        delete surfelementht;
        delete segmentht;

        surfelementht = new INDEX_3_CLOSED_HASHTABLE<int> (3 * GetNSE() + 1);
        segmentht = new INDEX_2_CLOSED_HASHTABLE<int> (3 * GetNSeg() + 1);

        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
            const Element2d & sel = surfelements[sei];
            if (sel.IsDeleted()) continue;

            int si = sel.GetIndex();

            for (int j = 0; j < sel.GetNP(); j++) {
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

        for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
            const Element2d & sel = surfelements[sei];
            if (sel.IsDeleted()) continue;

            INDEX_3 i3;
            i3.I1() = sel.PNum(1);
            i3.I2() = sel.PNum(2);
            i3.I3() = sel.PNum(3);
            i3.Sort();
            surfelementht -> Set(i3, sei); // war das wichtig ???    sel.GetIndex());
        }

        // int np = GetNP();

        if (dimension == 3) {
            for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
                points[pi].SetType(INNERPOINT);
            }
            if (GetNFD() == 0) {
                for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
                    const Element2d & sel = surfelements[sei];
                    if (sel.IsDeleted()) continue;
                    for (int j = 0; j < sel.GetNP(); j++) {
                        PointIndex pi = SurfaceElement(sei)[j];
                        points[pi].SetType(FIXEDPOINT);
                    }
                }
            }
            else {
                for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++) {
                    const Element2d & sel = surfelements[sei];
                    if (sel.IsDeleted()) continue;
                    for (int j = 0; j < sel.GetNP(); j++) {
                        PointIndex pi = sel[j];
                        int ns = surfacesonnode[pi].size();
                        if (ns == 1)
                            points[pi].SetType(SURFACEPOINT);
                        if (ns == 2)
                            points[pi].SetType(EDGEPOINT);
                        if (ns >= 3)
                            points[pi].SetType(FIXEDPOINT);
                    }
                }
            }
        }

        for (int i = 0; i < segments.size(); i++) {
            const Segment & seg = segments[i];
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
            const Segment & seg = segments[i];
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();

            segmentht -> Set(i2, i);
        }
    }

    void Mesh::FixPoints(const BitArray & fixpoints)
    {
        if (fixpoints.Size() != GetNP()) {
            std::cerr << "Mesh::FixPoints: sizes don't fit" << std::endl;
            return;
        }
        int np = GetNP();
        for (int i = 1; i <= np; i++) {
            if (fixpoints.Test(i)) {
                points.Elem(i).SetType(FIXEDPOINT);
            }
        }
    }

    void Mesh::FindOpenElements(int dom)
    {
        int np = GetNP();
        int ne = GetNE();
        int nse = GetNSE();

        Array<int, PointIndex::BASE> numonpoint(np);

        numonpoint = 0;

        for (ElementIndex ei = 0; ei < ne; ei++) {
            const Element & el = volelements[ei];
            if (dom == 0 || dom == el.GetIndex()) {
                if (el.GetNP() == 4) {
                    INDEX_4 i4(el[0], el[1], el[2], el[3]);
                    i4.Sort();
                    numonpoint[i4.I1()]++;
                    numonpoint[i4.I2()]++;
                }
                else {
                    for (int j = 0; j < el.GetNP(); j++) {
                        numonpoint[el[j]]++;
                    }
                }
            }
        }

        TABLE<ElementIndex, PointIndex::BASE> elsonpoint(numonpoint);
        for (ElementIndex ei = 0; ei < ne; ei++) {
            const Element & el = volelements[ei];
            if (dom == 0 || dom == el.GetIndex()) {
                if (el.GetNP() == 4) {
                    INDEX_4 i4(el[0], el[1], el[2], el[3]);
                    i4.Sort();
                    elsonpoint.Add(i4.I1(), ei);
                    elsonpoint.Add(i4.I2(), ei);
                }
                else {
                    for (int j = 0; j < el.GetNP(); j++) {
                        elsonpoint.Add(el[j], ei);
                    }
                }
            }
        }

        Array<char, 1> hasface(GetNFD());

        int i;
        for (i = 1; i <= GetNFD(); i++) {
            int domin = GetFaceDescriptor(i).DomainIn();
            int domout = GetFaceDescriptor(i).DomainOut();
            hasface[i] =
                    (dom == 0 && (domin != 0 || domout != 0)) ||
                    (dom != 0 && (domin == dom || domout == dom));
        }

        numonpoint = 0;
        for (SurfaceElementIndex sii = 0; sii < nse; sii++) {
            int ind = surfelements[sii].GetIndex();

            if (hasface[ind]) {
                const Element2d & hel = surfelements[sii];
                int mini = 0;
                for (int j = 1; j < hel.GetNP(); j++) {
                    if (hel[j] < hel[mini])
                        mini = j;
                }
                numonpoint[hel[mini]]++;
            }
        }

        TABLE<SurfaceElementIndex, PointIndex::BASE> selsonpoint(numonpoint);
        for (SurfaceElementIndex sii = 0; sii < nse; sii++) {
            int ind = surfelements[sii].GetIndex();

            if (hasface[ind]) {
                const Element2d & hel = surfelements[sii];
                int mini = 0;
                for (int j = 1; j < hel.GetNP(); j++) {
                    if (hel[j] < hel[mini])
                        mini = j;
                }
                selsonpoint.Add(hel[mini], sii);
            }
        }

        int ii;
        PointIndex pi;
        SurfaceElementIndex sei;
        Element2d hel;

        INDEX_3_CLOSED_HASHTABLE<INDEX_2> faceht(100);
        openelements.resize(0);

        for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
            if (selsonpoint[pi].size() + elsonpoint[pi].size()) {

                faceht.SetSize(2 * selsonpoint[pi].size() + 4 * elsonpoint[pi].size());

                FlatArray<SurfaceElementIndex> row = selsonpoint[pi];
                for (ii = 0; ii < row.size(); ii++) {
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
                                    ? PointIndex(PointIndex::BASE - 1)
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
                                    ? PointIndex(PointIndex::BASE - 1)
                                    : hel.PNum(4));
                            faceht.Set(i3, i2);
                        }
                    }
                }

                FlatArray<ElementIndex> rowel = elsonpoint[pi];
                for (ii = 0; ii < rowel.size(); ii++) {
                    const Element & el = VolumeElement(rowel[ii]);

                    if (dom == 0 || el.GetIndex() == dom) {
                        for (int j = 1; j <= el.GetNFaces(); j++) {
                            el.GetFace(j, hel);
                            hel.Invert();
                            hel.NormalizeNumbering();

                            if (hel[0] == pi) {
                                INDEX_3 i3(hel[0], hel[1], hel[2]);

                                if (faceht.Used(i3)) {
                                    INDEX_2 i2 = faceht.Get(i3);
                                    if (i2.I1() == el.GetIndex()) {
                                        i2.I1() = PointIndex::BASE - 1;
                                        faceht.Set(i3, i2);
                                    }
                                    else {
                                        if (i2.I1() == 0) {
                                            LOG_ERROR("more elements on face !");
                                            LOG_ERROR("   el = " << el);
                                            LOG_ERROR("  hel = " << hel);
                                            LOG_ERROR(" face = " << i3);
                                        }
                                    }
                                }
                                else {
                                    hel.Invert();
                                    hel.NormalizeNumbering();
                                    INDEX_3 i3(hel[0], hel[1], hel[2]);
                                    INDEX_2 i2(el.GetIndex(),
                                            (hel.GetNP() == 3)
                                            ? PointIndex(PointIndex::BASE - 1)
                                            : hel[3]);
                                    faceht.Set(i3, i2);
                                }
                            }
                        }
                    }
                }
                for (int i = 0; i < faceht.Size(); i++) {
                    if (faceht.UsedPos(i)) {
                        INDEX_3 i3;
                        INDEX_2 i2;
                        faceht.GetData(i, i3, i2);
                        if (i2.I1() != PointIndex::BASE - 1) {
                            Element2d tri;
                            tri.SetType((i2.I2() == PointIndex::BASE - 1) ? TRIG : QUAD);
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
        for (i = 0; i < openelements.size(); i++) {
            if (openelements[i].GetNP() == 3)
                cnt3++;
        }
        int cnt4 = openelements.size() - cnt3;
        if (openelements.size() > 0)
            LOG_WARNING(openelements.size() << " (" << cnt3 << " + " << cnt4 << ")" << " open elements");

        BuildBoundaryEdges();

        for (int i = 1; i <= openelements.size(); i++) {
            const Element2d & sel = openelements.Get(i);

            if (boundaryedges)
                for (int j = 1; j <= sel.GetNP(); j++) {
                    INDEX_2 i2;
                    i2.I1() = sel.PNumMod(j);
                    i2.I2() = sel.PNumMod(j + 1);
                    i2.Sort();
                    boundaryedges->Set(i2, 1);
                }

            for (int j = 1; j <= 3; j++) {
                PointIndex pi = sel.PNum(j);
                if (pi < points.End())
                    points[pi].SetType(FIXEDPOINT);
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

        LOG_DEBUG("Test Opensegments");
        for (int i = 1; i <= GetNSeg(); i++) {
            const Segment & seg = LineSegment(i);

            if (surfnr == 0 || seg.si == surfnr) {
                INDEX_2 key(seg[0], seg[1]);
                INDEX_2 data(seg.si, -i);

                if (faceht.Used(key)) {
                    std::cerr << "ERROR: Segment " << seg << " already used" << std::endl;
                }

                faceht.Set(key, data);
            }
        }

        for (int i = 1; i <= GetNSeg(); i++) {
            const Segment & seg = LineSegment(i);

            if (surfnr == 0 || seg.si == surfnr) {
                INDEX_2 key(seg[1], seg[0]);
                if (!faceht.Used(key)) {
                    std::cerr << "ERROR: Segment " << seg << " brother not used" << std::endl;
                }
            }
        }

        for (int i = 1; i <= GetNSE(); i++) {
            const Element2d & el = SurfaceElement(i);
            if (el.IsDeleted()) continue;

            if (surfnr == 0 || el.GetIndex() == surfnr) {
                for (int j = 1; j <= el.GetNP(); j++) {
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
                            LOG_WARNING("hash table si not fitting for segment: " <<
                                    seg.I1() << "-" << seg.I2() << " other = " << data.I2());
                        }
                    }
                    else {
                        std::swap(seg.I1(), seg.I2());
                        data.I1() = el.GetIndex();
                        data.I2() = i;

                        faceht.Set(seg, data);
                    }
                }
            }
        }

        std::cerr << "open segments: " << std::endl;
        opensegments.resize(0);
        for (int i = 1; i <= faceht.GetNBags(); i++) {
            for (int j = 1; j <= faceht.GetBagSize(i); j++) {
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
                        const Element2d & el = SurfaceElement(data.I2());
                        for (int k = 1; k <= el.GetNP(); k++) {
                            if (seg[0] == el.PNum(k))
                                seg.geominfo[0] = el.GeomInfoPi(k);
                            if (seg[1] == el.PNum(k))
                                seg.geominfo[1] = el.GeomInfoPi(k);
                        }
                        std::cerr << "trig seg: ";
                    }
                    else {
                        // segment due to line
                        const Segment & lseg = LineSegment(-data.I2());
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
        LOG_DEBUG(opensegments.size() << " open segments found");

        for (int i = 1; i <= points.size(); i++) {
            points.Elem(i).SetType(SURFACEPOINT);
        }

        for (int i = 1; i <= GetNSeg(); i++) {
            const Segment & seg = LineSegment(i);
            points[seg[0]].SetType(EDGEPOINT);
            points[seg[1]].SetType(EDGEPOINT);
        }
        for (int i = 1; i <= GetNOpenSegments(); i++) {
            const Segment & seg = GetOpenSegment(i);
            points[seg[0]].SetType(EDGEPOINT);
            points[seg[1]].SetType(EDGEPOINT);
        }
    }

    void Mesh::RemoveOneLayerSurfaceElements()
    {
        int i, j;
        int np = GetNP();

        FindOpenSegments();
        BitArray frontpoints(np);

        frontpoints.Clear();
        for (i = 1; i <= GetNOpenSegments(); i++) {
            const Segment & seg = GetOpenSegment(i);
            frontpoints.Set(seg[0]);
            frontpoints.Set(seg[1]);
        }

        for (i = 1; i <= GetNSE(); i++) {
            Element2d & sel = surfelements.Elem(i);
            int remove = 0;
            for (j = 1; j <= sel.GetNP(); j++) {
                if (frontpoints.Test(sel.PNum(j)))
                    remove = 1;
            }
            if (remove)
                sel.PNum(1) = 0;
        }

        for (i = surfelements.size(); i >= 1; i--) {
            if (surfelements.Elem(i).PNum(1) == 0) {
                surfelements.Elem(i) = surfelements.Last();
                surfelements.DeleteLast();
            }
        }

        RebuildSurfaceElementLists();

        timestamp = NextTimeStamp();
        //  Compress();
    }

    void Mesh::FreeOpenElementsEnvironment(int layers)
    {
        int i, j, k;
        PointIndex pi;
        const int large = 9999;
        Array<int, PointIndex::BASE> dist(GetNP());

        dist = large;

        for (int i = 1; i <= GetNOpenElements(); i++) {
            const Element2d & face = OpenElement(i);
            for (j = 0; j < face.GetNP(); j++) {
                dist[face[j]] = 1;
            }
        }

        for (k = 1; k <= layers; k++) {
            for (i = 1; i <= GetNE(); i++) {
                const Element & el = VolumeElement(i);
                if (el[0] == -1 || el.IsDeleted()) continue;

                int elmin = large;
                for (j = 0; j < el.GetNP(); j++) {
                    if (dist[el[j]] < elmin)
                        elmin = dist[el[j]];
                }
                if (elmin < large) {
                    for (j = 0; j < el.GetNP(); j++) {
                        if (dist[el[j]] > elmin + 1)
                            dist[el[j]] = elmin + 1;
                    }
                }
            }
        }
        int cntfree = 0;
        for (i = 1; i <= GetNE(); i++) {
            Element & el = VolumeElement(i);
            if (el[0] == -1 || el.IsDeleted()) continue;

            int elmin = large;
            for (j = 0; j < el.GetNP(); j++) {
                if (dist[el[j]] < elmin)
                    elmin = dist[el[j]];
            }
            el.flags.fixed = elmin > layers;

            if (elmin <= layers)
                cntfree++;
        }

        LOG_DEBUG("free: " << cntfree << ", fixed: " << GetNE() - cntfree);

        for (pi = PointIndex::BASE;
                pi < GetNP() + PointIndex::BASE; pi++) {
            if (dist[pi] > layers + 1)
                points[pi].SetType(FIXEDPOINT);
        }
    }

    void Mesh::SetLocalH(const Point3d & pmin, const Point3d & pmax, double grading)
    {
        Point3d c = Center(pmin, pmax);
        double d = max3(pmax.X() - pmin.X(),
                pmax.Y() - pmin.Y(),
                pmax.Z() - pmin.Z());
        d /= 2;
        Point3d pmin2 = c - Vec3d(d, d, d);
        Point3d pmax2 = c + Vec3d(d, d, d);


        delete lochfunc;
        lochfunc = new LocalH(pmin2, pmax2, grading);
    }

    void Mesh::RestrictLocalH(const Point3d & p, double hloc)
    {
        if (hloc < hmin)
            hloc = hmin;

        //std::cout << "restrict h in " << p << " to " << hloc <<std::endl;
        if (!lochfunc) {
            LOG_WARNING("RestrictLocalH called, creating mesh-size tree");

            Point3d boxmin, boxmax;
            GetBox(boxmin, boxmax);
            SetLocalH(boxmin, boxmax, 0.8);
        }

        lochfunc -> SetH(p, hloc);
    }

    void Mesh::RestrictLocalHLine(const Point3d & p1,
            const Point3d & p2,
            double hloc)
    {
        if (hloc < hmin)
            hloc = hmin;

        // std::cout << "restrict h along " << p1 << " - " << p2 << " to " << hloc <<std::endl;
        int i;
        int steps = int (Dist(p1, p2) / hloc) + 2;
        Vec3d v(p1, p2);

        for (i = 0; i <= steps; i++) {
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
            return maxhdomain.Get(dom);
        else
            return 1e10;
    }

    void Mesh::SetMaxHDomain(const Array<double> & mhd)
    {
        maxhdomain.resize(mhd.size());
        for (int i = 1; i <= mhd.size(); i++) {
            maxhdomain.Elem(i) = mhd.Get(i);
        }
    }

    double Mesh::GetH(const Point3d & p) const
    {
        double hmin = hglob;
        if (lochfunc) {
            double hl = lochfunc->GetH(p);
            if (hl < hglob)
                hmin = hl;
        }
        return hmin;
    }

    double Mesh::GetMinH(const Point3d & pmin, const Point3d & pmax)
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
            const Element2d & el = SurfaceElement(i);
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

        LOG_DEBUG("minh = " << minh << " avh = " << (hsum / n) << " maxh = " << maxh);
        return (hsum / n);
    }

    void Mesh::CalcLocalH(double grading)
    {
        if (!lochfunc) {
            Point3d pmin, pmax;
            GetBox(pmin, pmax);
            SetLocalH(pmin, pmax, grading);
        }

        LOG_DEBUG("CalcLocalH: " << GetNP() << " points, "
                << GetNE() << " elements, "
                << GetNSE() << " surface elements.");

        for (int i = 0; i < GetNSE(); i++) {
            const Element2d & el = surfelements[i];
            int j;

            if (el.GetNP() == 3) {
                double hel = -1;
                for (j = 1; j <= 3; j++) {
                    const Point3d & p1 = points[el.PNumMod(j)];
                    const Point3d & p2 = points[el.PNumMod(j + 1)];

                    if (!ident -> UsedSymmetric(el.PNumMod(j),
                            el.PNumMod(j + 1))) {
                        double hedge = Dist(p1, p2);
                        if (hedge > hel)
                            hel = hedge;
                    }
                }

                if (hel > 0) {
                    const Point3d & p1 = points[el.PNum(1)];
                    const Point3d & p2 = points[el.PNum(2)];
                    const Point3d & p3 = points[el.PNum(3)];
                    lochfunc->SetH(Center(p1, p2, p3), hel);
                }
            }
            else {
                {
                    const Point3d & p1 = points[el.PNum(1)];
                    const Point3d & p2 = points[el.PNum(2)];
                    lochfunc->SetH(Center(p1, p2), 2 * Dist(p1, p2));
                }
                {
                    const Point3d & p1 = points[el.PNum(3)];
                    const Point3d & p2 = points[el.PNum(4)];
                    lochfunc->SetH(Center(p1, p2), 2 * Dist(p1, p2));
                }
            }
        }

        for (int i = 0; i < GetNSeg(); i++) {
            const Segment & seg = segments[i];
            const Point3d & p1 = points[seg[0]];
            const Point3d & p2 = points[seg[1]];

            if (!ident -> UsedSymmetric(seg[0], seg[1])) {
                lochfunc->SetH(Center(p1, p2), Dist(p1, p2));
            }
        }
    }

    void Mesh::CalcLocalHFromPointDistances(double grading)
    {
        LOG_DEBUG("Calculating local h from point distances");

        if (!lochfunc) {
            Point3d pmin, pmax;
            GetBox(pmin, pmax);

            SetLocalH(pmin, pmax, grading);
        }

        PointIndex i, j;
        double hl;


        for (i = PointIndex::BASE;
                i < GetNP() + PointIndex::BASE; i++) {
            for (j = i + 1; j < GetNP() + PointIndex::BASE; j++) {
                const Point3d & p1 = points[i];
                const Point3d & p2 = points[j];
                hl = Dist(p1, p2);
                RestrictLocalH(p1, hl);
                RestrictLocalH(p2, hl);
            }
        }


    }

    void Mesh::CalcLocalHFromSurfaceCurvature(double grading, double elperr)
    {
        LOG_DEBUG("Calculating local h from surface curvature");

        if (!lochfunc) {
            Point3d pmin, pmax;
            GetBox(pmin, pmax);

            SetLocalH(pmin, pmax, grading);
        }

        INDEX_2_HASHTABLE<int> edges(3 * GetNP() + 2);
        INDEX_2_HASHTABLE<int> bedges(GetNSeg() + 2);
        int i, j;

        for (i = 1; i <= GetNSeg(); i++) {
            const Segment & seg = LineSegment(i);
            INDEX_2 i2(seg[0], seg[1]);
            i2.Sort();
            bedges.Set(i2, 1);
        }
        for (i = 1; i <= GetNSE(); i++) {
            const Element2d & sel = SurfaceElement(i);
            if (!sel.PNum(1))
                continue;
            for (j = 1; j <= 3; j++) {
                INDEX_2 i2(sel.PNumMod(j), sel.PNumMod(j + 1));
                i2.Sort();
                if (bedges.Used(i2)) continue;

                if (edges.Used(i2)) {
                    int other = edges.Get(i2);

                    const Element2d & elother = SurfaceElement(other);

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

                    double rad = ComputeCylinderRadius(Point(i2.I1()),
                            Point(i2.I2()),
                            Point(pi3),
                            Point(pi4));

                    RestrictLocalHLine(Point(i2.I1()), Point(i2.I2()), rad / elperr);
                }
                else
                    edges.Set(i2, i);
            }
        }


        // Restrict h due to line segments

        for (i = 1; i <= GetNSeg(); i++) {
            const Segment & seg = LineSegment(i);
            const Point3d & p1 = Point(seg[0]);
            const Point3d & p2 = Point(seg[1]);
            RestrictLocalH(Center(p1, p2), Dist(p1, p2));
        }
    }

    void Mesh::RestrictLocalH(resthtype rht, int nr, double loch)
    {
        int i;
        switch (rht) {
            case RESTRICTH_FACE:
            {
                for (i = 1; i <= GetNSE(); i++) {
                    const Element2d & sel = SurfaceElement(i);
                    if (sel.GetIndex() == nr)
                        RestrictLocalH(RESTRICTH_SURFACEELEMENT, i, loch);
                }
                break;
            }
            case RESTRICTH_EDGE:
            {
                for (i = 1; i <= GetNSeg(); i++) {
                    const Segment & seg = LineSegment(i);
                    if (seg.edgenr == nr)
                        RestrictLocalH(RESTRICTH_SEGMENT, i, loch);
                }
                break;
            }
            case RESTRICTH_POINT:
            {
                RestrictLocalH(Point(nr), loch);
                break;
            }

            case RESTRICTH_SURFACEELEMENT:
            {
                const Element2d & sel = SurfaceElement(nr);
                Point3d p = Center(Point(sel.PNum(1)),
                        Point(sel.PNum(2)),
                        Point(sel.PNum(3)));
                RestrictLocalH(p, loch);
                break;
            }
            case RESTRICTH_SEGMENT:
            {
                const Segment & seg = LineSegment(nr);
                RestrictLocalHLine(Point(seg[0]), Point(seg[1]), loch);
                break;
            }
        }
    }

    void Mesh::LoadLocalMeshSize(const char * meshsizefilename)
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
            LOG_ERROR("Error loading mesh size file: " << meshsizefilename << "....  Skipping!");
            return;
        }

        LOG_DEBUG("Load local mesh-size file: " << meshsizefilename);

        int nmsp = 0;
        int nmsl = 0;

        msf >> nmsp;
        if (!msf.good())
            throw std::runtime_error("Mesh-size file error: No points found\n");

        if (nmsp > 0)
            LOG_DEBUG("Number of mesh-size restriction points: " << nmsp);

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
            LOG_DEBUG("Number of mesh-size restriction lines: " << nmsl);

        for (int i = 0; i < nmsl; i++) {
            Point3d p1, p2;
            double hi;
            msf >> p1.X() >> p1.Y() >> p1.Z();
            msf >> p2.X() >> p2.Y() >> p2.Z();
            msf >> hi;
            if (!msf.good())
                throw std::runtime_error("Mesh-size file error: Number of line definitions don't match specified list size\n");
            RestrictLocalHLine(p1, p2, hi);
        }

        msf.close();
    }

    void Mesh::GetBox(Point3d & pmin, Point3d & pmax, int dom) const
    {
        if (points.size() == 0) {
            pmin = pmax = Point3d(0, 0, 0);
            return;
        }

        if (dom <= 0) {
            pmin = Point3d(1e10, 1e10, 1e10);
            pmax = Point3d(-1e10, -1e10, -1e10);

            for (PointIndex pi = points.Begin(); pi < points.End(); pi++) {
                pmin.SetToMin(points[pi]);
                pmax.SetToMax(points[pi]);
            }
        }
        else {
            int j, nse = GetNSE();
            SurfaceElementIndex sei;

            pmin = Point3d(1e10, 1e10, 1e10);
            pmax = Point3d(-1e10, -1e10, -1e10);
            for (sei = 0; sei < nse; sei++) {
                const Element2d & el = surfelements[sei];
                if (el.IsDeleted()) continue;

                if (dom == -1 || el.GetIndex() == dom) {
                    for (j = 0; j < 3; j++) {
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

    void Mesh::GetBox(Point3d & pmin, Point3d & pmax, POINTTYPE ptyp) const
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

    double Mesh::ElementError(int eli, const MeshingParameters & mp) const
    {
        const Element & el = volelements.Get(eli);
        return CalcTetBadness(points.Get(el[0]), points.Get(el[1]),
                points.Get(el[2]), points.Get(el[3]), -1, mp);
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
        Array<PointIndex, PointIndex::BASE, PointIndex> op2np(GetNP());
        Array<MeshPoint> hpoints;
        BitArrayChar<PointIndex::BASE> pused(GetNP());

        for (int i = 0; i < volelements.size(); i++) {
            if (volelements[i][0] <= PointIndex::BASE - 1 ||
                    volelements[i].IsDeleted()) {
                volelements.Delete(i);
                i--;
            }
        }

        for (int i = 0; i < surfelements.size(); i++) {
            if (surfelements[i].IsDeleted()) {
                surfelements.Delete(i);
                i--;
            }
        }
        for (int i = 0; i < segments.size(); i++) {
            if (segments[i][0] <= PointIndex::BASE - 1) {
                segments.Delete(i);
                i--;
            }
        }
        pused.Clear();
        for (int i = 0; i < volelements.size(); i++) {
            const Element & el = volelements[i];
            for (int j = 0; j < el.GetNP(); j++) {
                pused.Set(el[j]);
            }
        }

        for (int i = 0; i < surfelements.size(); i++) {
            const Element2d & el = surfelements[i];
            for (int j = 0; j < el.GetNP(); j++) {
                pused.Set(el[j]);
            }
        }

        for (int i = 0; i < segments.size(); i++) {
            const Segment & seg = segments[i];
            pused.Set(seg[0]);
            pused.Set(seg[1]);
        }

        for (int i = 0; i < openelements.size(); i++) {
            const Element2d & el = openelements[i];
            for (int j = 0; j < el.GetNP(); j++) {
                pused.Set(el[j]);
            }
        }

        for (int i = 0; i < lockedpoints.size(); i++) {
            pused.Set(lockedpoints[i]);
        }

        int npi = PointIndex::BASE - 1;

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

        for (int i = 1; i <= volelements.size(); i++) {
            Element & el = VolumeElement(i);
            for (int j = 0; j < el.GetNP(); j++) {
                el[j] = op2np[el[j]];
            }
        }

        for (int i = 1; i <= surfelements.size(); i++) {
            Element2d & el = SurfaceElement(i);
            for (int j = 0; j < el.GetNP(); j++) {
                el[j] = op2np[el[j]];
            }
        }

        for (int i = 0; i < segments.size(); i++) {
            Segment & seg = segments[i];
            seg[0] = op2np[seg[0]];
            seg[1] = op2np[seg[1]];
        }

        for (int i = 1; i <= openelements.size(); i++) {
            Element2d & el = openelements.Elem(i);
            for (int j = 0; j < el.GetNP(); j++) {
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
            const Element2d & sel = OpenElement(i);

            for (int j = 1; j <= sel.GetNP(); j++) {
                i2.I1() = sel.PNumMod(j);
                i2.I2() = sel.PNumMod(j + 1);

                int sign = (i2.I2() > i2.I1()) ? 1 : -1;
                i2.Sort();
                if (!edges.Used(i2))
                    edges.Set(i2, 0);
                edges.Set(i2, edges.Get(i2) + sign);
            }
        }

        for (int i = 1; i <= edges.GetNBags(); i++) {
            for (int j = 1; j <= edges.GetBagSize(i); j++) {
                int cnt = 0;
                edges.GetData(i, j, i2, cnt);
                if (cnt) {
                    LOG_ERROR("Edge " << i2.I1() << " - " << i2.I2() << " multiple times in surface mesh");
                    i2s = i2;
                    i2s.Sort();
                    for (int k = 1; k <= nf; k++) {
                        const Element2d & sel = OpenElement(k);
                        for (int l = 1; l <= sel.GetNP(); l++) {
                            edge.I1() = sel.PNumMod(l);
                            edge.I2() = sel.PNumMod(l + 1);
                            edge.Sort();

                            if (edge == i2s)
                                LOG_ERROR("  edge of element " << sel);
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
        int i, j, k;

        Point3d pmin, pmax;
        GetBox(pmin, pmax);
        Box3dTree setree(pmin, pmax);
        Array<int> inters;

        bool overlap = 0;
        bool incons_layers = 0;


        for (i = 1; i <= GetNSE(); i++) {
            SurfaceElement(i).badel = 0;
        }

        for (i = 1; i <= GetNSE(); i++) {
            const Element2d & tri = SurfaceElement(i);

            Point3d tpmin(Point(tri[0]));
            Point3d tpmax(tpmin);

            for (k = 1; k < tri.GetNP(); k++) {
                tpmin.SetToMin(Point(tri[k]));
                tpmax.SetToMax(Point(tri[k]));
            }
            Vec3d diag(tpmin, tpmax);

            tpmax = tpmax + 0.1 * diag;
            tpmin = tpmin - 0.1 * diag;

            setree.Insert(tpmin, tpmax, i);
        }

        for (i = 1; i <= GetNSE(); i++) {
            const Element2d & tri = SurfaceElement(i);

            Point3d tpmin(Point(tri[0]));
            Point3d tpmax(tpmin);

            for (k = 1; k < tri.GetNP(); k++) {
                tpmin.SetToMin(Point(tri[k]));
                tpmax.SetToMax(Point(tri[k]));
            }

            setree.GetIntersecting(tpmin, tpmax, inters);

            for (j = 1; j <= inters.size(); j++) {
                const Element2d & tri2 = SurfaceElement(inters.Get(j));

                if (points[tri[0]].GetLayer() != points[tri2[0]].GetLayer())
                    continue;

                if (points[tri[0]].GetLayer() != points[tri[1]].GetLayer() ||
                        points[tri[0]].GetLayer() != points[tri[2]].GetLayer()) {
                    incons_layers = 1;
                    LOG_WARNING("inconsistent layers in triangle");
                }


                const meshit::Point<3> *trip1[3], *trip2[3];
                for (k = 1; k <= 3; k++) {
                    trip1[k - 1] = &Point(tri.PNum(k));
                    trip2[k - 1] = &Point(tri2.PNum(k));
                }

                if (IntersectTriangleTriangle(&trip1[0], &trip2[0])) {
                    overlap = 1;
                    LOG_WARNING("Intersecting elements " << i << " and " << inters.Get(j));
                    LOG_DEBUG(" el1 = " << tri);
                    LOG_DEBUG(" el2 = " << tri2);

                    for (k = 1; k <= 3; k++)
                        LOG_DEBUG_CONT(tri.PNum(k) << "  ");
                    LOG_DEBUG("");
                    for (k = 1; k <= 3; k++)
                        LOG_DEBUG_CONT(tri2.PNum(k) << "  ");
                    LOG_DEBUG("");

                    for (k = 0; k <= 2; k++)
                        LOG_DEBUG_CONT(*trip1[k] << "   ");
                    LOG_DEBUG("");
                    for (k = 0; k <= 2; k++)
                        LOG_DEBUG_CONT(*trip2[k] << "   ");
                    LOG_DEBUG("");

                    LOG_DEBUG("Face1 = " << GetFaceDescriptor(tri.GetIndex()));
                    LOG_DEBUG("Face1 = " << GetFaceDescriptor(tri2.GetIndex()));

                    SurfaceElement(i).badel = 1;
                    SurfaceElement(inters.Get(j)).badel = 1;
                }
            }
        }

        // bug 'fix'
        if (incons_layers) overlap = 0;

        return overlap;
    }

    int Mesh::CheckVolumeMesh() const
    {
        LOG_DEBUG("Checking volume mesh");

        int ne = GetNE();
        DenseMatrix dtrans(3, 3);
        int i, j;

        LOG_DEBUG(ne << " elements");
        for (i = 1; i <= ne; i++) {
            Element & el = (Element&) VolumeElement(i);
            el.flags.badel = 0;
            int nip = el.GetNIP();
            for (j = 1; j <= nip; j++) {
                el.GetTransformation(j, Points(), dtrans);
                double det = dtrans.Det();
                if (det > 0) {
                    LOG_ERROR("Element " << i << " has wrong orientation");
                    el.flags.badel = 1;
                }
            }
        }
        return 0;
    }

    bool Mesh::LegalTrig(const Element2d & el)
    {
        return 1;
        if (/* hp */ 1) // needed for old, simple hp-refinement
        {
            // trigs with 2 or more segments are illegal
            int i;
            int nseg = 0;

            if (!segmentht) {
                std::cerr << "no segmentht allocated" << std::endl;
                return 0;
            }

            //      Point3d cp(0.5, 0.5, 0.5);
            for (i = 1; i <= 3; i++) {
                INDEX_2 i2(el.PNumMod(i), el.PNumMod(i + 1));
                i2.Sort();
                if (segmentht -> Used(i2))
                    nseg++;
            }
            if (nseg >= 2)
                return 0;
        }
        return 1;
    }

    bool Mesh::LegalTet2(Element & el)
    {
        // Test, whether 4 points have a common surface plus
        // at least 4 edges at the boundary

        if (!boundaryedges) BuildBoundaryEdges();

        // non-tets are always legal
        if (el.GetType() != TET) {
            el.SetLegal(1);
            return 1;
        }

        POINTTYPE pointtype[4];
        for (int i = 0; i < 4; i++) {
            pointtype[i] = points[el[i]].Type();
        }

        // element has at least 2 inner points ---> legal
        int cnti = 0;
        for (int j = 0; j < 4; j++) {
            if (pointtype[j] == INNERPOINT) {
                cnti++;
                if (cnti >= 2) {
                    el.SetLegal(1);
                    return 1;
                }
            }
        }
        // which faces are boundary faces ?
        int bface[4];
        for (int i = 0; i < 4; i++) {
            bface[i] = surfelementht->Used(INDEX_3::Sort(el[gftetfacesa[i][0]],
                    el[gftetfacesa[i][1]],
                    el[gftetfacesa[i][2]]));
        }

        int bedge[4][4];
        int segedge[4][4];
        static const int pi3map[4][4] = {
            { -1, 2, 1, 1},
            { 2, -1, 0, 0},
            { 1, 0, -1, 0},
            { 1, 0, 0, -1}
        };

        static const int pi4map[4][4] = {
            { -1, 3, 3, 2},
            { 3, -1, 3, 2},
            { 3, 3, -1, 1},
            { 2, 2, 1, -1}
        };

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < i; j++) {
                bool sege = false, be = false;

                int pos = boundaryedges -> Position(INDEX_2::Sort(el[i], el[j]));
                if (pos) {
                    be = true;
                    if (boundaryedges -> GetData(pos) == 2)
                        sege = true;
                }

                segedge[j][i] = segedge[i][j] = sege;
                bedge[j][i] = bedge[i][j] = be;
            }
        }
        // two boundary faces and no edge is illegal
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 4; j++) {
                if (bface[i] && bface[j]) {
                    if (!segedge[pi3map[i][j]][pi4map[i][j]]) {
                        // 2 boundary faces withoud edge in between
                        el.SetLegal(0);
                        return 0;
                    }
                }
            }
        }
        // three boundary edges meeting in a Surface point
        for (int i = 0; i < 4; i++) {
            if (pointtype[i] == SURFACEPOINT) {
                bool alledges = 1;
                for (int j = 0; j < 4; j++) {
                    if (j != i && !bedge[i][j]) {
                        alledges = 0;
                        break;
                    }
                }
                if (alledges) {
                    el.SetLegal(0);
                    return 0;
                }
            }
        }

        for (int fnr = 0; fnr < 4; fnr++) {
            if (!bface[fnr]) {
                for (int i = 0; i < 4; i++) {
                    if (i != fnr) {
                        int pi1 = pi3map[i][fnr];
                        int pi2 = pi4map[i][fnr];

                        if (pointtype[i] == SURFACEPOINT) {
                            // two connected edges on surface, but no face
                            if (bedge[i][pi1] && bedge[i][pi2]) {
                                el.SetLegal(0);
                                return 0;
                            }
                        }
                        if (pointtype[i] == EDGEPOINT) {
                            // connected surface edge and edge edge, but no face
                            if ((bedge[i][pi1] && segedge[i][pi2]) ||
                                    (bedge[i][pi2] && segedge[i][pi1])) {
                                el.SetLegal(0);
                                return 0;
                            }
                        }

                    }
                }
            }
        }
        el.SetLegal(1);
        return 1;
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

    void Mesh::SurfaceMeshOrientation()
    {
        int i, j;
        int nse = GetNSE();

        BitArray used(nse);
        used.Clear();
        INDEX_2_HASHTABLE<int> edges(nse + 1);

        bool haschanged = 0;

        const Element2d & tri = SurfaceElement(1);
        for (j = 1; j <= 3; j++) {
            INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j + 1));
            edges.Set(i2, 1);
        }
        used.Set(1);

        bool unused;
        do {
            bool changed;
            do {
                changed = 0;
                for (i = 1; i <= nse; i++) {
                    if (!used.Test(i)) {
                        Element2d & el = surfelements.Elem(i);
                        int found = 0, foundrev = 0;
                        for (j = 1; j <= 3; j++) {
                            INDEX_2 i2(el.PNumMod(j), el.PNumMod(j + 1));
                            if (edges.Used(i2))
                                foundrev = 1;
                            std::swap(i2.I1(), i2.I2());
                            if (edges.Used(i2))
                                found = 1;
                        }

                        if (found || foundrev) {
                            if (foundrev) {
                                std::swap(el.PNum(2), el.PNum(3));
                            }
                            changed = 1;
                            for (j = 1; j <= 3; j++) {
                                INDEX_2 i2(el.PNumMod(j), el.PNumMod(j + 1));
                                edges.Set(i2, 1);
                            }
                            used.Set(i);
                        }
                    }
                }
                if (changed)
                    haschanged = 1;
            } while (changed);


            unused = 0;
            for (i = 1; i <= nse; i++) {
                if (!used.Test(i)) {
                    unused = 1;
                    const Element2d & tri = SurfaceElement(i);
                    for (j = 1; j <= 3; j++) {
                        INDEX_2 i2(tri.PNumMod(j), tri.PNumMod(j + 1));
                        edges.Set(i2, 1);
                    }
                    used.Set(i);
                    break;
                }
            }
        } while (unused);

        if (haschanged)
            timestamp = NextTimeStamp();
    }

    void Mesh::Split2Tets()
    {
        LOG_DEBUG("Split To Tets");
        bool has_prisms = 0;

        int oldne = GetNE();
        for (int i = 1; i <= oldne; i++) {
            Element el = VolumeElement(i);

            if (el.GetType() == PRISM) {
                // prism, to 3 tets

                // make minimal node to node 1
                int minpi = 0;
                PointIndex minpnum;
                minpnum = GetNP() + 1;

                for (int j = 1; j <= 6; j++) {
                    if (el.PNum(j) < minpnum) {
                        minpnum = el.PNum(j);
                        minpi = j;
                    }
                }

                if (minpi >= 4) {
                    for (int j = 1; j <= 3; j++) {
                        std::swap(el.PNum(j), el.PNum(j + 3));
                    }
                    minpi -= 3;
                }

                while (minpi > 1) {
                    int hi = 0;
                    for (int j = 0; j <= 3; j += 3) {
                        hi = el.PNum(1 + j);
                        el.PNum(1 + j) = el.PNum(2 + j);
                        el.PNum(2 + j) = el.PNum(3 + j);
                        el.PNum(3 + j) = hi;
                    }
                    minpi--;
                }

                /*
                  version 1: edge from pi2 to pi6,
                  version 2: edge from pi3 to pi5,
                 */

                static const int ntets[2][12] = {
                    { 1, 4, 5, 6, 1, 2, 3, 6, 1, 2, 5, 6},
                    { 1, 4, 5, 6, 1, 2, 3, 5, 3, 1, 5, 6}
                };

                const int * min2pi;

                if (min2(el.PNum(2), el.PNum(6)) < min2(el.PNum(3), el.PNum(5))) {
                    min2pi = &ntets[0][0];
                }
                else {
                    min2pi = &ntets[1][0];
                }

                int firsttet = 1;
                for (int j = 1; j <= 3; j++) {
                    Element nel(TET);
                    for (int k = 1; k <= 4; k++) {
                        nel.PNum(k) = el.PNum(min2pi[4 * j + k - 5]);
                    }
                    nel.SetIndex(el.GetIndex());

                    int legal = 1;
                    for (int k = 1; k <= 3; k++) {
                        for (int l = k + 1; l <= 4; l++) {
                            if (nel.PNum(k) == nel.PNum(l))
                                legal = 0;
                        }
                    }
                    if (legal) {
                        if (firsttet) {
                            VolumeElement(i) = nel;
                            firsttet = 0;
                        }
                        else {
                            AddVolumeElement(nel);
                        }
                    }
                }
            }

            else if (el.GetType() == HEX) {
                // hex to A) 2 prisms or B) to 5 tets

                // make minimal node to node 1
                int minpi = 0;
                PointIndex minpnum;
                minpnum = GetNP() + 1;

                for (int j = 1; j <= 8; j++) {
                    if (el.PNum(j) < minpnum) {
                        minpnum = el.PNum(j);
                        minpi = j;
                    }
                }

                if (minpi >= 5) {
                    for (int j = 1; j <= 4; j++) {
                        std::swap(el.PNum(j), el.PNum(j + 4));
                    }
                    minpi -= 4;
                }

                while (minpi > 1) {
                    int hi = 0;
                    for (int j = 0; j <= 4; j += 4) {
                        hi = el.PNum(1 + j);
                        el.PNum(1 + j) = el.PNum(2 + j);
                        el.PNum(2 + j) = el.PNum(3 + j);
                        el.PNum(3 + j) = el.PNum(4 + j);
                        el.PNum(4 + j) = hi;
                    }
                    minpi--;
                }

                static const int to_prisms[3][12] = {
                    { 0, 1, 2, 4, 5, 6, 0, 2, 3, 4, 6, 7},
                    { 0, 1, 5, 3, 2, 6, 0, 5, 4, 3, 6, 7},
                    { 0, 7, 4, 1, 6, 5, 0, 3, 7, 1, 2, 6},
                };

                const int * min2pi = 0;
                if (min2(el[4], el[6]) < min2(el[5], el[7])) {
                    min2pi = &to_prisms[0][0];
                }
                else if (min2(el[3], el[6]) < min2(el[2], el[7])) {
                    min2pi = &to_prisms[1][0];
                }
                else if (min2(el[1], el[6]) < min2(el[2], el[5])) {
                    min2pi = &to_prisms[2][0];
                }

                if (min2pi) {
                    has_prisms = 1;
                    for (int j = 0; j < 2; j++) {
                        Element nel(PRISM);
                        for (int k = 0; k < 6; k++) {
                            nel[k] = el[min2pi[6 * j + k]];
                        }
                        nel.SetIndex(el.GetIndex());

                        if (j == 0)
                            VolumeElement(i) = nel;
                        else
                            AddVolumeElement(nel);
                    }
                }
                else {
                    // split to 5 tets

                    static const int to_tets[20] = {
                        1, 2, 0, 5,
                        3, 0, 2, 7,
                        4, 5, 7, 0,
                        6, 7, 5, 2,
                        0, 2, 7, 5
                    };

                    for (int j = 0; j < 5; j++) {
                        Element nel(TET);
                        for (int k = 0; k < 4; k++) {
                            nel[k] = el[to_tets[4 * j + k]];
                        }
                        nel.SetIndex(el.GetIndex());

                        if (j == 0)
                            VolumeElement(i) = nel;
                        else
                            AddVolumeElement(nel);
                    }

                }
            }
            else if (el.GetType() == PYRAMID) {
                // pyramid, to 2 tets
                static const int ntets[2][8] = {
                    { 1, 2, 3, 5, 1, 3, 4, 5},
                    { 1, 2, 4, 5, 4, 2, 3, 5}
                };

                const int * min2pi;

                if (min2(el[0], el[2]) < min2(el[1], el[3]))
                    min2pi = &ntets[0][0];
                else
                    min2pi = &ntets[1][0];

                bool firsttet = 1;
                for (int j = 0; j < 2; j++) {
                    Element nel(TET);
                    for (int k = 0; k < 4; k++) {
                        nel[k] = el[min2pi[4 * j + k] - 1];
                    }
                    nel.SetIndex(el.GetIndex());

                    bool legal = 1;
                    for (int k = 0; k < 3; k++) {
                        for (int l = k + 1; l < 4; l++) {
                            if (nel[k] == nel[l])
                                legal = 0;
                        }
                    }
                    if (legal) {
                        if (firsttet)
                            VolumeElement(i) = nel;
                        else
                            AddVolumeElement(nel);

                        firsttet = 0;
                    }
                }
            }
        }

        int oldnse = GetNSE();
        for (int i = 1; i <= oldnse; i++) {
            Element2d el = SurfaceElement(i);
            if (el.GetNP() == 4) {

                static const int ntris[2][6] = {
                    { 1, 2, 3, 1, 3, 4},
                    { 1, 2, 4, 4, 2, 3}
                };

                const int * min2pi;

                if (min2(el.PNum(1), el.PNum(3)) <
                        min2(el.PNum(2), el.PNum(4)))
                    min2pi = &ntris[0][0];
                else
                    min2pi = &ntris[1][0];

                int firsttri = 1;
                for (int j = 1; j <= 2; j++) {
                    Element2d nel(3);
                    for (int k = 1; k <= 3; k++) {
                        nel.PNum(k) = el.PNum(min2pi[3 * j + k - 4]);
                    }
                    nel.SetIndex(el.GetIndex());

                    int legal = 1;
                    for (int k = 1; k <= 2; k++) {
                        for (int l = k + 1; l <= 3; l++) {
                            if (nel.PNum(k) == nel.PNum(l))
                                legal = 0;
                        }
                    }
                    if (legal) {
                        if (firsttri) {
                            SurfaceElement(i) = nel;
                            firsttri = 0;
                        }
                        else {
                            AddSurfaceElement(nel);
                        }
                    }
                }
            }
        }

        if (has_prisms)
            Split2Tets();

        else {
            for (int i = 1; i <= GetNE(); i++) {
                Element & el = VolumeElement(i);
                const Point3d & p1 = Point(el.PNum(1));
                const Point3d & p2 = Point(el.PNum(2));
                const Point3d & p3 = Point(el.PNum(3));
                const Point3d & p4 = Point(el.PNum(4));

                double vol = (Vec3d(p1, p2) *
                        Cross(Vec3d(p1, p3), Vec3d(p1, p4)));
                if (vol > 0)
                    std::swap(el.PNum(3), el.PNum(4));
            }



            UpdateTopology();
            timestamp = NextTimeStamp();
        }
    }

    void Mesh::BuildElementSearchTree()
    {
        if (elementsearchtreets == GetTimeStamp()) return;

        //#pragma omp critical (buildsearchtree)
        {
            if (elementsearchtreets != GetTimeStamp()) {

                LOG_DEBUG("Rebuild element searchtree");

                delete elementsearchtree;
                elementsearchtree = NULL;

                int ne = (dimension == 2) ? GetNSE() : GetNE();

                if (ne) {
                    if (dimension == 2) {
                        Box<3> box(Box<3>::EMPTY_BOX);
                        for (SurfaceElementIndex sei = 0; sei < ne; sei++) {
                            box.Add(points[surfelements[sei].PNums()]);
                        }
                        box.Increase(1.01 * box.Diam());
                        elementsearchtree = new Box3dTree(box);

                        for (SurfaceElementIndex sei = 0; sei < ne; sei++) {
                            box.Set(points[surfelements[sei].PNums()]);
                            elementsearchtree -> Insert(box, sei + 1);
                        }
                    }
                    else {
                        Box<3> box(Box<3>::EMPTY_BOX);
                        for (ElementIndex ei = 0; ei < ne; ei++) {
                            box.Add(points[volelements[ei].PNums()]);
                        }
                        box.Increase(1.01 * box.Diam());
                        elementsearchtree = new Box3dTree(box);

                        for (ElementIndex ei = 0; ei < ne; ei++) {
                            box.Set(points[volelements[ei].PNums()]);
                            elementsearchtree -> Insert(box, ei + 1);
                        }
                    }

                    elementsearchtreets = GetTimeStamp();
                }
            }
        }
    }

    bool Mesh::PointContainedIn2DElement(const Point3d & p,
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
            const Element2d & el = SurfaceElement(element);

            const Point3d & p1 = Point(el.PNum(1));
            const Point3d & p2 = Point(el.PNum(2));
            const Point3d & p3 = Point(el.PNum(3));
            const Point3d & p4 = Point(el.PNum(4));

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

            double c0, c1, c2; // ,rt; 
            lami[2] = 0.;
            double eps = 1.E-12;

            if (fabs(d.X()) <= eps && fabs(d.Y()) <= eps) {
                //Solve Linear System
                lami[0] = (c.Y()*(p.X() - a.X()) - c.X()*(p.Y() - a.Y())) /
                        (b.X() * c.Y() - b.Y() * c.X());
                lami[1] = (-b.Y()*(p.X() - a.X()) + b.X()*(p.Y() - a.Y())) /
                        (b.X() * c.Y() - b.Y() * c.X());
            }
            else
                if (fabs(dxb) <= eps) {
                lami[1] = (dxp - dxa) / dxc;
                if (fabs(b.X() - d.X() * lami[1]) >= eps)
                    lami[0] = (p.X() - a.X() - c.X() * lami[1]) / (b.X() + d.X() * lami[1]);
                else
                    lami[0] = (p.Y() - a.Y() - c.Y() * lami[1]) / (b.Y() + d.Y() * lami[1]);
            }
            else
                if (fabs(dxc) <= eps) {
                lami[0] = (dxp - dxa) / dxb;
                if (fabs(c.X() - d.X() * lami[0]) >= eps)
                    lami[1] = (p.X() - a.X() - b.X() * lami[0]) / (c.X() + d.X() * lami[0]);
                else
                    lami[1] = (p.Y() - a.Y() - b.Y() * lami[0]) / (c.Y() + d.Y() * lami[0]);
            }
            else //Solve quadratic equation
            {
                if (fabs(d.X()) >= eps) {
                    c2 = d.X() * dxc;
                    c1 = d.X() * dxc - c.X() * dxb - d.X()*(dxp - dxa);
                    c0 = -b.X()*(dxp - dxa) - (a.X() - p.X()) * dxb;
                }
                else {
                    c2 = d.Y() * dxc;
                    c1 = d.Y() * dxc - c.Y() * dxb - d.Y()*(dxp - dxa);
                    c0 = -b.Y()*(dxp - dxa) - (a.Y() - p.Y()) * dxb;
                }

                double rt = c1 * c1 - 4 * c2*c0;
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

        }
        else {
            //	  SurfaceElement(element).GetTets (loctets);
            loctrigs.resize(1);
            loctrigs.Elem(1) = SurfaceElement(element);

            for (int j = 1; j <= loctrigs.size(); j++) {
                const Element2d & el = loctrigs.Get(j);


                const Point3d & p1 = Point(el.PNum(1));
                const Point3d & p2 = Point(el.PNum(2));
                const Point3d & p3 = Point(el.PNum(3));

                col1 = p2 - p1;
                col2 = p3 - p1;
                col3 = Cross(col1, col2);
                //col3 = Vec3d(0, 0, 1);
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

    bool Mesh::PointContainedIn3DElement(const Point3d & p,
            double lami[3],
            const int element) const
    {
        const double eps = 1.e-4;
        const Element & el = VolumeElement(element);

        meshit::Point<3> lam = 0.0;

        if (el.GetType() == TET || el.GetType() == TET10) {
            lam = 0.25;
        }
        else if (el.GetType() == PRISM) {
            lam(0) = 0.33;
            lam(1) = 0.33;
            lam(2) = 0.5;
        }
        else if (el.GetType() == PYRAMID) {
            lam(0) = 0.4;
            lam(1) = 0.4;
            lam(2) = 0.2;
        }
        else if (el.GetType() == HEX) {
            lam = 0.5;
        }


        Vec<3> deltalam, rhs;
        meshit::Point<3> x;
        Mat<3, 3> Jac, Jact;

        double delta = 1;

        bool retval;

        int i = 0;

        const int maxits = 30;
        while (delta > 1e-16 && i < maxits) {
            curvedelems->CalcElementTransformation(lam, element - 1, x, Jac);
            rhs = p - x;
            Jac.Solve(rhs, deltalam);

            lam += deltalam;
            delta = deltalam.Length2();
            i++;
        }

        if (i == maxits)
            return false;

        for (i = 0; i < 3; i++) {
            lami[i] = lam(i);
        }
        if (el.GetType() == TET || el.GetType() == TET10) {
            retval = (lam(0) > -eps &&
                    lam(1) > -eps &&
                    lam(2) > -eps &&
                    lam(0) + lam(1) + lam(2) < 1 + eps);
        }
        else if (el.GetType() == PRISM) {
            retval = (lam(0) > -eps &&
                    lam(1) > -eps &&
                    lam(2) > -eps &&
                    lam(2) < 1 + eps &&
                    lam(0) + lam(1) < 1 + eps);
        }
        else if (el.GetType() == PYRAMID) {
            retval = (lam(0) > -eps &&
                    lam(1) > -eps &&
                    lam(2) > -eps &&
                    lam(0) + lam(2) < 1 + eps &&
                    lam(1) + lam(2) < 1 + eps);
        }
        else if (el.GetType() == HEX) {
            retval = (lam(0) > -eps && lam(0) < 1 + eps &&
                    lam(1) > -eps && lam(1) < 1 + eps &&
                    lam(2) > -eps && lam(2) < 1 + eps);
        }
        else
            throw std::runtime_error("Da haun i wos vagessn");

        return retval;
    }

    bool Mesh::PointContainedIn3DElementOld(const Point3d & p,
            double lami[3],
            const int element) const
    {
        Vec3d col1, col2, col3;
        Vec3d rhs, sol;
        const double eps = 1.e-4;

        Array<Element> loctets;

        VolumeElement(element).GetTets(loctets);

        for (int j = 1; j <= loctets.size(); j++) {
            const Element & el = loctets.Get(j);

            const Point3d & p1 = Point(el.PNum(1));
            const Point3d & p2 = Point(el.PNum(2));
            const Point3d & p3 = Point(el.PNum(3));
            const Point3d & p4 = Point(el.PNum(4));

            Box3d box;
            box.SetPoint(p1);
            box.AddPoint(p2);
            box.AddPoint(p3);
            box.AddPoint(p4);
            if (!box.IsIn(p))
                continue;

            col1 = p2 - p1;
            col2 = p3 - p1;
            col3 = p4 - p1;
            rhs = p - p1;

            SolveLinearSystem(col1, col2, col3, rhs, sol);

            if (sol.X() >= -eps && sol.Y() >= -eps && sol.Z() >= -eps &&
                    sol.X() + sol.Y() + sol.Z() <= 1 + eps) {
                Array<Element> loctetsloc;
                Array<meshit::Point<3> > pointsloc;

                VolumeElement(element).GetTetsLocal(loctetsloc);
                VolumeElement(element).GetNodesLocalNew(pointsloc);

                const Element & le = loctetsloc.Get(j);


                Point3d pp =
                        pointsloc.Get(le.PNum(1))
                        + sol.X() * Vec3d(pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(2)))
                        + sol.Y() * Vec3d(pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(3)))
                        + sol.Z() * Vec3d(pointsloc.Get(le.PNum(1)), pointsloc.Get(le.PNum(4)));

                lami[0] = pp.X();
                lami[1] = pp.Y();
                lami[2] = pp.Z();
                return true;
            }
        }
        return false;
    }

    int Mesh::GetElementOfPoint(
            const meshit::Point<3> & p,
            double lami[3],
            bool build_searchtree,
            const int index,
            const bool allowindex)
    {
        if (index != -1) {
            Array<int> dummy(1);
            dummy[0] = index;
            return GetElementOfPoint(p, lami, &dummy, build_searchtree, allowindex);
        }
        else
            return GetElementOfPoint(p, lami, NULL, build_searchtree, allowindex);
    }

    int Mesh::GetElementOfPoint(
            const meshit::Point<3> & p,
            double lami[3],
            const Array<int> * const indices,
            bool build_searchtree,
            const bool allowindex)
    {
        if (dimension == 2) {
            int ne;


            if (ps_startelement != 0 && ps_startelement <= GetNSE() && PointContainedIn2DElement(p, lami, ps_startelement))
                return ps_startelement;

            Array<int> locels;
            if (elementsearchtree || build_searchtree) {
                // update if necessary:
                BuildElementSearchTree();
                elementsearchtree->GetIntersecting(p, p, locels);
                ne = locels.size();
            }
            else
                ne = GetNSE();

            for (int i = 1; i <= ne; i++) {
                int ii;

                if (elementsearchtree)
                    ii = locels.Get(i);
                else
                    ii = i;

                if (ii == ps_startelement) continue;

                if (indices != NULL && indices->size() > 0) {
                    bool contained = indices->Contains(SurfaceElement(ii).GetIndex());
                    if ((allowindex && !contained) || (!allowindex && contained)) continue;
                }

                if (PointContainedIn2DElement(p, lami, ii)) return ii;

            }
            return 0;
        }
        else {
            // int i, j;
            int ne;

            if (ps_startelement != 0 && PointContainedIn3DElement(p, lami, ps_startelement))
                return ps_startelement;

            Array<int> locels;
            if (elementsearchtree || build_searchtree) {
                // update if necessary:
                BuildElementSearchTree();
                elementsearchtree->GetIntersecting(p, p, locels);
                ne = locels.size();
            }
            else
                ne = GetNE();

            for (int i = 1; i <= ne; i++) {
                int ii;

                if (elementsearchtree)
                    ii = locels.Get(i);
                else
                    ii = i;
                if (ii == ps_startelement) continue;

                if (indices != NULL && indices->size() > 0) {
                    bool contained = indices->Contains(VolumeElement(ii).GetIndex());
                    if ((allowindex && !contained) || (!allowindex && contained)) continue;
                }

                if (PointContainedIn3DElement(p, lami, ii)) {
                    ps_startelement = ii;
                    return ii;
                }
            }

            // Not found, try uncurved variant:
            for (int i = 1; i <= ne; i++) {
                int ii;

                if (elementsearchtree)
                    ii = locels.Get(i);
                else
                    ii = i;

                if (indices != NULL && indices->size() > 0) {
                    bool contained = indices->Contains(VolumeElement(ii).GetIndex());
                    if ((allowindex && !contained) || (!allowindex && contained)) continue;
                }


                if (PointContainedIn3DElementOld(p, lami, ii)) {
                    ps_startelement = ii;
                    LOG_WARNING("found element of point " << p << " only for uncurved mesh");
                    return ii;
                }
            }

            return 0;
        }
    }

    int Mesh::GetSurfaceElementOfPoint(
            const meshit::Point<3> & p,
            double lami[3],
            bool build_searchtree,
            const int index,
            const bool allowindex)
    {
        if (index != -1) {
            Array<int> dummy(1);
            dummy[0] = index;
            return GetSurfaceElementOfPoint(p, lami, &dummy, build_searchtree, allowindex);
        }
        else
            return GetSurfaceElementOfPoint(p, lami, NULL, build_searchtree, allowindex);
    }

    int Mesh::GetSurfaceElementOfPoint(const meshit::Point<3> & p,
            double lami[3],
            const Array<int> * const indices,
            bool build_searchtree,
            const bool allowindex)
    {
        if (dimension == 2) {
            throw std::runtime_error("GetSurfaceElementOfPoint not yet implemented for 2D meshes");
        }
        else {
            double vlam[3];
            int velement = GetElementOfPoint(p, vlam, NULL, build_searchtree, allowindex);

            Array<int> faces;
            topology->GetElementFaces(velement, faces);

            for (int i = 0; i < faces.size(); i++) {
                faces[i] = topology->GetFace2SurfaceElement(faces[i]);
            }

            for (int i = 0; i < faces.size(); i++) {
                if (faces[i] == 0)
                    continue;

                if (indices && indices->size() != 0) {
                    if (indices->Contains(SurfaceElement(faces[i]).GetIndex()) &&
                            PointContainedIn2DElement(p, lami, faces[i], true))
                        return faces[i];
                }
                else {
                    if (PointContainedIn2DElement(p, lami, faces[i], true)) {
                        return faces[i];
                    }
                }
            }

        }

        return 0;
    }

    void Mesh::GetIntersectingVolEls(const Point3d& p1, const Point3d& p2,
            Array<int> & locels) const
    {
        elementsearchtree->GetIntersecting(p1, p2, locels);
    }

    void Mesh::SplitIntoParts()
    {
        int i, j, dom;
        int ne = GetNE();
        int np = GetNP();
        int nse = GetNSE();

        BitArray surfused(nse);
        BitArray pused(np);

        surfused.Clear();

        dom = 0;

        while (1) {
            int cntd = 1;

            dom++;

            pused.Clear();

            int found = 0;
            for (i = 1; i <= nse; i++) {
                if (!surfused.Test(i)) {
                    SurfaceElement(i).SetIndex(dom);
                    for (j = 1; j <= 3; j++) {
                        pused.Set(SurfaceElement(i).PNum(j));
                    }
                    found = 1;
                    cntd = 1;
                    surfused.Set(i);
                    break;
                }
            }
            if (!found)
                break;

            int change;
            do {
                change = 0;
                for (i = 1; i <= nse; i++) {
                    int is = 0, isnot = 0;
                    for (j = 1; j <= 3; j++) {
                        if (pused.Test(SurfaceElement(i).PNum(j)))
                            is = 1;
                        else
                            isnot = 1;
                    }
                    if (is && isnot) {
                        change = 1;
                        for (j = 1; j <= 3; j++) {
                            pused.Set(SurfaceElement(i).PNum(j));
                        }
                    }

                    if (is) {
                        if (!surfused.Test(i)) {
                            surfused.Set(i);
                            SurfaceElement(i).SetIndex(dom);
                            cntd++;
                        }
                    }
                }


                for (i = 1; i <= ne; i++) {
                    int is = 0, isnot = 0;
                    for (j = 1; j <= 4; j++) {
                        if (pused.Test(VolumeElement(i).PNum(j)))
                            is = 1;
                        else
                            isnot = 1;
                    }
                    if (is && isnot) {
                        change = 1;
                        for (j = 1; j <= 4; j++) {
                            pused.Set(VolumeElement(i).PNum(j));
                        }
                    }

                    if (is) {
                        VolumeElement(i).SetIndex(dom);
                    }
                }
            } while (change);

            LOG_DEBUG("domain " << dom << " has " << cntd << " surfaceelements");
        }

        ClearFaceDescriptors();
        for (i = 1; i <= dom; i++) {
            AddFaceDescriptor(FaceDescriptor(0, i, 0, 0));
        }
        CalcSurfacesOfNode();
        timestamp = NextTimeStamp();
    }

    void Mesh::SplitSeparatedFaces()
    {
        LOG_DEBUG("SplitSeparateFaces");
        int fdi;
        int np = GetNP();

        BitArray usedp(np);
        Array<SurfaceElementIndex> els_of_face;

        fdi = 1;
        while (fdi <= GetNFD()) {
            GetSurfaceElementsOfFace(fdi, els_of_face);

            if (els_of_face.size() == 0) continue;

            SurfaceElementIndex firstel = els_of_face[0];

            usedp.Clear();
            for (int j = 1; j <= SurfaceElement(firstel).GetNP(); j++) {
                usedp.Set(SurfaceElement(firstel).PNum(j));
            }

            bool changed;
            do {
                changed = false;

                for (int i = 0; i < els_of_face.size(); i++) {
                    const Element2d & el = SurfaceElement(els_of_face[i]);

                    bool has = 0;
                    bool hasno = 0;
                    for (int j = 0; j < el.GetNP(); j++) {
                        if (usedp.Test(el[j]))
                            has = true;
                        else
                            hasno = true;
                    }

                    if (has && hasno)
                        changed = true;

                    if (has)
                        for (int j = 0; j < el.GetNP(); j++) {
                            usedp.Set(el[j]);
                        }
                }
            } while (changed);

            int nface = 0;
            for (int i = 0; i < els_of_face.size(); i++) {
                Element2d & el = SurfaceElement(els_of_face[i]);

                int hasno = 0;
                for (int j = 1; j <= el.GetNP(); j++) {
                    if (!usedp.Test(el.PNum(j)))
                        hasno = 1;
                }

                if (hasno) {
                    if (!nface) {
                        FaceDescriptor nfd = GetFaceDescriptor(fdi);
                        nface = AddFaceDescriptor(nfd);
                    }

                    el.SetIndex(nface);
                }
            }

            // reconnect list
            if (nface) {
                facedecoding[nface - 1].firstelement = -1;
                facedecoding[fdi - 1].firstelement = -1;

                for (int i = 0; i < els_of_face.size(); i++) {
                    int ind = SurfaceElement(els_of_face[i]).GetIndex();
                    SurfaceElement(els_of_face[i]).next = facedecoding[ind - 1].firstelement;
                    facedecoding[ind - 1].firstelement = els_of_face[i];
                }
            }

            fdi++;
        }
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

    void Mesh::GetSurfaceElementsOfFace(int facenr, Array<SurfaceElementIndex> & sei) const
    {
        /* Philippose - 01/10/2009
        Commented out the following lines, and activated the originally 
        commented out lines above because of a bug which causes corruption 
        of the variable "facedecoding" when a mesh is converted to second order
         */
        sei.resize(0);

        SurfaceElementIndex si = facedecoding[facenr - 1].firstelement;
        while (si != -1) {
            const Element2d & se = SurfaceElement(si);
            if (se.GetIndex() == facenr
                    && se[0] >= PointIndex::BASE
                    && !se.IsDeleted()) {
                sei.push_back(si);
            }
            si = se.next;
        }
    }

    void Mesh::CalcMinMaxAngle(double badellimit, double * retvalues)
    {
        int i, j;
        int lpi1, lpi2, lpi3, lpi4;
        double phimax = 0, phimin = 10;
        double facephimax = 0, facephimin = 10;
        int illegaltets = 0, negativetets = 0, badtets = 0;

        for (i = 1; i <= GetNE(); i++) {
            bool badel = false;

            Element & el = VolumeElement(i);

            if (el.GetType() != TET) {
                VolumeElement(i).flags.badel = 0;
                continue;
            }

            if (el.Volume(Points()) < 0) {
                badel = true;
                negativetets++;
            }

            if (!LegalTet(el)) {
                badel = true;
                illegaltets++;
                std::cerr << "illegal tet: " << i << " ";
                for (j = 1; j <= el.GetNP(); j++) {
                    std::cerr << el.PNum(j) << " ";
                }
                std::cerr << std::endl;
            }

            // angles between faces
            for (lpi1 = 1; lpi1 <= 3; lpi1++) {
                for (lpi2 = lpi1 + 1; lpi2 <= 4; lpi2++) {
                    lpi3 = 1;
                    while (lpi3 == lpi1 || lpi3 == lpi2)
                        lpi3++;
                    lpi4 = 10 - lpi1 - lpi2 - lpi3;

                    const Point3d & p1 = Point(el.PNum(lpi1));
                    const Point3d & p2 = Point(el.PNum(lpi2));
                    const Point3d & p3 = Point(el.PNum(lpi3));
                    const Point3d & p4 = Point(el.PNum(lpi4));

                    Vec3d n(p1, p2);
                    n /= n.Length();
                    Vec3d v1(p1, p3);
                    Vec3d v2(p1, p4);

                    v1 -= (n * v1) * n;
                    v2 -= (n * v2) * n;

                    double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
                    double phi = acos(cosphi);
                    if (phi > phimax) phimax = phi;
                    if (phi < phimin) phimin = phi;

                    if ((180 / M_PI) * phi > badellimit)
                        badel = true;
                }
            }

            // angles in faces
            for (j = 1; j <= 4; j++) {
                Element2d face;
                el.GetFace(j, face);
                for (lpi1 = 1; lpi1 <= 3; lpi1++) {
                    lpi2 = lpi1 % 3 + 1;
                    lpi3 = lpi2 % 3 + 1;

                    const Point3d & p1 = Point(el.PNum(lpi1));
                    const Point3d & p2 = Point(el.PNum(lpi2));
                    const Point3d & p3 = Point(el.PNum(lpi3));

                    Vec3d v1(p1, p2);
                    Vec3d v2(p1, p3);
                    double cosphi = (v1 * v2) / (v1.Length() * v2.Length());
                    double phi = acos(cosphi);
                    if (phi > facephimax) facephimax = phi;
                    if (phi < facephimin) facephimin = phi;

                    if ((180 / M_PI) * phi > badellimit)
                        badel = true;

                }
            }


            VolumeElement(i).flags.badel = badel;
            if (badel) badtets++;
        }

        if (!GetNE()) {
            phimin = phimax = facephimin = facephimax = 0;
        }

        if (!retvalues) {
            LOG_WARNING("between planes:  phimin = " << (180 / M_PI) * phimin <<
                    " phimax = " << (180 / M_PI) * phimax);
            LOG_WARNING("inside planes:   phimin = " << (180 / M_PI) * facephimin <<
                    " phimax = " << (180 / M_PI) * facephimax);
        }
        else {
            retvalues[0] = (180 / M_PI) * facephimin;
            retvalues[1] = (180 / M_PI) * facephimax;
            retvalues[2] = (180 / M_PI) * phimin;
            retvalues[3] = (180 / M_PI) * phimax;
        }
        LOG_DEBUG("negative tets: " << negativetets);
        LOG_DEBUG("illegal tets:  " << illegaltets);
        LOG_DEBUG("bad tets:      " << badtets);
    }

    int Mesh::MarkIllegalElements()
    {
        int cnt = 0;
        int i;

        for (i = 1; i <= GetNE(); i++) {
            LegalTet(VolumeElement(i));
            cnt += VolumeElement(i).Illegal();
        }
        return cnt;
    }

    void Mesh::InitPointCurve(double red, double green, double blue) const
    {
        pointcurves_startpoint.push_back(pointcurves.size());
        pointcurves_red.push_back(red);
        pointcurves_green.push_back(green);
        pointcurves_blue.push_back(blue);
    }

    void Mesh::AddPointCurvePoint(const Point3d & pt) const
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

    Point3d & Mesh::GetPointCurvePoint(int curve, int n) const
    {
        return pointcurves[pointcurves_startpoint[curve] + n];
    }

    void Mesh::GetPointCurveColor(int curve, double & red, double & green, double & blue) const
    {
        red = pointcurves_red[curve];
        green = pointcurves_green[curve];
        blue = pointcurves_blue[curve];
    }

    void Mesh::ComputeNVertices()
    {
        int i, j, nv;
        int ne = GetNE();
        int nse = GetNSE();

        numvertices = 0;
        for (i = 1; i <= ne; i++) {
            const Element & el = VolumeElement(i);
            nv = el.GetNV();
            for (j = 0; j < nv; j++) {
                if (el[j] > numvertices)
                    numvertices = el[j];
            }
        }
        for (i = 1; i <= nse; i++) {
            const Element2d & el = SurfaceElement(i);
            nv = el.GetNV();
            for (j = 1; j <= nv; j++) {
                if (el.PNum(j) > numvertices)
                    numvertices = el.PNum(j);
            }
        }

        numvertices += 1 - PointIndex::BASE;
    }

    int Mesh::GetNV() const
    {
        if (numvertices < 0)
            return GetNP();
        else
            return numvertices;
    }

    void Mesh::SetNP(int np)
    {
        points.resize(np);
        //  ptyps.SetSize(np);

        int mlold = mlbetweennodes.size();
        mlbetweennodes.resize(np);
        if (np > mlold)
            for (int i = mlold + PointIndex::BASE;
                    i < np + PointIndex::BASE; i++) {
                mlbetweennodes[i].I1() = PointIndex::BASE - 1;
                mlbetweennodes[i].I2() = PointIndex::BASE - 1;
            }

        GetIdentifications().SetMaxPointNr(np + PointIndex::BASE - 1);
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

    bool Mesh::PureTetMesh() const
    {
        for (ElementIndex ei = 0; ei < GetNE(); ei++) {
            if (VolumeElement(ei).GetNP() != 4)
                return 0;
        }
        return 1;
    }

    void Mesh::UpdateTopology()
    {
        topology->Update();
        clusters->Update();
    }

    void Mesh::SetMaterial(int domnr, const char * mat)
    {
        if (domnr > materials.size()) {
            int olds = materials.size();
            materials.resize(domnr);
            for (int i = olds; i < domnr; i++) {
                materials[i] = 0;
            }
        }
        materials.Elem(domnr) = new char[strlen(mat) + 1];
        strcpy(materials.Elem(domnr), mat);
    }

    const char * Mesh::GetMaterial(int domnr) const
    {
        if (domnr <= materials.size())
            return materials.Get(domnr);
        return 0;
    }

    void Mesh::SetNBCNames(int nbcn)
    {
        if (bcnames.size())
            for (int i = 0; i < bcnames.size(); i++) {
                if (bcnames[i]) delete bcnames[i];
            }
        bcnames.resize(nbcn);
        bcnames = 0;
    }

    void Mesh::SetBCName(int bcnr, const std::string & abcname)
    {
        if (bcnames[bcnr]) delete bcnames[bcnr];
        if (abcname != "default")
            bcnames[bcnr] = new std::string(abcname);
        else
            bcnames[bcnr] = 0;
    }

    const std::string & Mesh::GetBCName(int bcnr) const
    {
        static std::string defaultstring = "default";

        if (!bcnames.size())
            return defaultstring;
        if (bcnames[bcnr])
            return *bcnames[bcnr];
        else
            return defaultstring;
    }

    void Mesh::SetUserData(const char * id, Array<int> & data)
    {
        if (userdata_int.Used(id))
            delete userdata_int.Get(id);

        Array<int> * newdata = new Array<int>(data);

        userdata_int.Set(id, newdata);
    }

    bool Mesh::GetUserData(const char * id, Array<int> & data, int shift) const
    {
        if (userdata_int.Used(id)) {
            if (data.size() < (*userdata_int.Get(id)).size() + shift)
                data.resize((*userdata_int.Get(id)).size() + shift);
            for (int i = 0; i < (*userdata_int.Get(id)).size(); i++) {
                data[i + shift] = (*userdata_int.Get(id))[i];
            }
            return true;
        }
        else {
            data.resize(0);
            return false;
        }
    }

    void Mesh::SetUserData(const char * id, Array<double> & data)
    {
        if (userdata_double.Used(id))
            delete userdata_double.Get(id);

        Array<double> * newdata = new Array<double>(data);

        userdata_double.Set(id, newdata);
    }

    bool Mesh::GetUserData(const char * id, Array<double> & data, int shift) const
    {
        if (userdata_double.Used(id)) {
            if (data.size() < (*userdata_double.Get(id)).size() + shift)
                data.resize((*userdata_double.Get(id)).size() + shift);
            for (int i = 0; i < (*userdata_double.Get(id)).size(); i++) {
                data[i + shift] = (*userdata_double.Get(id))[i];
            }
            return true;
        }
        else {
            data.resize(0);
            return false;
        }
    }

    void Mesh::PrintMemInfo(std::ostream & ost) const
    {
        ost << "Mesh Mem:" << std::endl;

        ost << GetNP() << " Points, of size "
                << sizeof (Point3d) << " + " << sizeof (POINTTYPE) << " = "
                << GetNP() * (sizeof (Point3d) + sizeof (POINTTYPE)) << std::endl;

        ost << GetNSE() << " Surface elements, of size "
                << sizeof (Element2d) << " = "
                << GetNSE() * sizeof (Element2d) << std::endl;

        ost << GetNE() << " Volume elements, of size "
                << sizeof (Element) << " = "
                << GetNE() * sizeof (Element) << std::endl;

        ost << "surfs on node:";
        surfacesonnode.PrintMemInfo(ost);

        ost << "boundaryedges: ";
        if (boundaryedges)
            boundaryedges->PrintMemInfo(ost);

        ost << "surfelementht: ";
        if (surfelementht)
            surfelementht->PrintMemInfo(ost);
    }
}
