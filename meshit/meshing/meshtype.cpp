#include <stdexcept>

#include "meshtype.hpp"
#include "meshclass.hpp"

namespace meshit {

    Segment::Segment()
    {
        pnums[0] = -1;
        pnums[1] = -1;
        edgenr = -1;
        seginfo = 0;

        si = -1;

        domin = -1;
        domout = -1;
        tlosurf = -1;

        surfnr1 = -1;
        surfnr2 = -1;
        pnums[2] = -1;
        meshdocval = 0;
    }

    Segment::Segment(const Segment& other) :
            edgenr(other.edgenr),
            seginfo(other.seginfo),
            si(other.si),
            domin(other.domin),
            domout(other.domout),
            tlosurf(other.tlosurf),
            surfnr1(other.surfnr1),
            surfnr2(other.surfnr2),
            epgeominfo(),
            meshdocval(other.meshdocval),
            hp_elnr(other.hp_elnr)
    {
        for (int j = 0; j < 3; j++) {
            pnums[j] = other.pnums[j];
        }
        epgeominfo[0] = other.epgeominfo[0];
        epgeominfo[1] = other.epgeominfo[1];
    }

    Segment& Segment::operator=(const Segment& other)
    {
        if (&other != this) {
            pnums[0] = other[0];
            pnums[1] = other[1];
            edgenr = other.edgenr;
            seginfo = other.seginfo;
            si = other.si;
            domin = other.domin;
            domout = other.domout;
            tlosurf = other.tlosurf;
            surfnr1 = other.surfnr1;
            surfnr2 = other.surfnr2;
            epgeominfo[0] = other.epgeominfo[0];
            epgeominfo[1] = other.epgeominfo[1];
            pnums[2] = other.pnums[2];
            meshdocval = other.meshdocval;
            hp_elnr = other.hp_elnr;
        }

        return *this;
    }

    std::ostream& operator<<(std::ostream& s, const Segment& seg)
    {
        s << seg[0] << " - " << seg[1] << " domin = " << seg.domin << ", domout = " << seg.domout
        << " si = " << seg.si << ", edgenr = " << seg.edgenr;
        return s;
    }

    Element2d::Element2d()
    {
        for (int i = 0; i < 3; i++) {
            pnum[i] = 0;
        }
        index = 0;
        deleted = 0;
    }

    void Element2d::GetBox(const Array<MeshPoint>& points, Box3d& box) const
    {
        box.SetPoint(points[pnum[0] - 1]);
        for (unsigned i = 1; i < 3; i++) {
            box.AddPoint(points[pnum[i] - 1]);
        }
    }

    bool Element2d::operator==(const Element2d& el2) const
    {
        bool retval = true;
        for (int i = 0; retval && i < 3; i++) {
            retval = (el2[i] == (*this)[i]);
        }

        return retval;
    }

    Array<IntegrationPointData*> ipdtrig;

    void Element2d::GetTransformation(class DenseMatrix& pmat, class DenseMatrix& trans) const
    {
        ComputeIntegrationPointData();
        DenseMatrix* dshapep = &ipdtrig[0]->dshape;

        CalcABt(pmat, *dshapep, trans);
    }

    double Element2d::CalcJacobianBadness(const Array<MeshPoint>& points) const
    {
        DenseMatrix trans(2, 2);
        DenseMatrix pmat;

        pmat.SetSize(2, 3);

        for (size_t i = 0; i < 3; i++) {
            const Point3d& p = points[PNum(i + 1)];
            pmat.Elem(0, i) = p.Y();
            pmat.Elem(1, i) = -p.X();
        }

        GetTransformation(pmat, trans);

        // Frobenius norm
        double frob = 0;
        for (size_t j = 1; j <= 4; j++) {
            double d = trans.Get(j);
            frob += d * d;
        }
        frob = 0.5 * sqrt(frob);

        double err;
        double det = trans.Det();
        if (det <= 0)
            err = 1e12;
        else
            err = frob * frob / det;

        return err;
    }

    void Element2d::ComputeIntegrationPointData() const
    {
        if (ipdtrig.size()) return;

        IntegrationPointData* ipd = new IntegrationPointData;
        ipd->p.X() = 1.0 / 3.0;
        ipd->p.Y() = 1.0 / 3.0;
        ipd->p.Z() = 0;
        ipd->weight = 0.5;

        ipd->shape.SetSize(3);
        ipd->dshape.SetSize(2, 3);

        ipd->shape[0] = 1.0 / 3.0;
        ipd->shape[1] = 1.0 / 3.0;
        ipd->shape[2] = 1.0 / 3.0;

        ipd->dshape.Elem(0, 0) = -1;
        ipd->dshape.Elem(0, 1) = 1;
        ipd->dshape.Elem(0, 2) = 0;
        ipd->dshape.Elem(1, 0) = -1;
        ipd->dshape.Elem(1, 1) = 0;
        ipd->dshape.Elem(1, 2) = 1;

        ipdtrig.push_back(ipd);
    }

    std::ostream& operator<<(std::ostream& s, const Element2d& el)
    {
        for (size_t j = 0; j < 3; j++) {
            s << " " << el.PNum(j + 1);
        }
        return s;
    }

    FaceDescriptor::FaceDescriptor()
    {
        surfnr = domin = domout = bcprop = 0;
        tlosurf = -1;
        firstelement = -1;
    }

    FaceDescriptor::FaceDescriptor(const FaceDescriptor& other)
            : surfnr(other.surfnr), domin(other.domin), domout(other.domout),
              tlosurf(other.tlosurf), bcprop(other.bcprop)
    {
        firstelement = -1;
    }

    FaceDescriptor::FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi)
    {
        surfnr = surfnri;
        domin = domini;
        domout = domouti;
        tlosurf = tlosurfi;
        bcprop = surfnri;
        firstelement = -1;
    }

    FaceDescriptor::FaceDescriptor(const Segment& seg)
    {
        surfnr = seg.si;
        domin = seg.domin + 1;
        domout = seg.domout + 1;
        tlosurf = seg.tlosurf + 1;
        bcprop = 0;
        firstelement = -1;
    }

    std::ostream& operator<<(std::ostream& s, const FaceDescriptor& fd)
    {
        s << "surfnr = " << fd.SurfNr()
        << ", domin = " << fd.DomainIn()
        << ", domout = " << fd.DomainOut()
        << ", tlosurf = " << fd.TLOSurface()
        << ", bcprop = " << fd.BCProperty();
        return s;
    }

    Identifications::Identifications(Mesh& amesh)
            : mesh(amesh)
    {
        identifiedpoints = new INDEX_2_HASHTABLE<int>(100);
        identifiedpoints_nr = new INDEX_3_HASHTABLE<int>(100);
        maxidentnr = 0;
    }

    Identifications::~Identifications()
    {
        delete identifiedpoints;
        delete identifiedpoints_nr;
    }

    void Identifications::Add(PointIndex pi1, PointIndex pi2, int identnr)
    {
        INDEX_2 pair(pi1, pi2);
        identifiedpoints->Set(pair, identnr);

        INDEX_3 tripl(pi1, pi2, identnr);
        identifiedpoints_nr->Set(tripl, 1);

        if (identnr > maxidentnr) maxidentnr = identnr;

        if (identnr + 1 > idpoints_table.Size())
            idpoints_table.ChangeSize(identnr + 1);
        idpoints_table.Add(identnr, pair);

        //  timestamp = NextTimeStamp();
    }

    int Identifications::Get(PointIndex pi1, PointIndex pi2) const
    {
        INDEX_2 pair(pi1, pi2);
        if (identifiedpoints->Used(pair))
            return identifiedpoints->Get(pair);
        else
            return 0;
    }

    bool Identifications::Get(PointIndex pi1, PointIndex pi2, int nr) const
    {
        INDEX_3 tripl(pi1, pi2, nr);
        if (identifiedpoints_nr->Used(tripl))
            return 1;
        else
            return 0;
    }

    void Identifications::GetMap(size_t identnr, Array<int>& identmap, bool symmetric) const
    {
        identmap.resize(mesh.GetNP());
        identmap = 0;

        if (identnr > 0) {
            for (size_t i = 0; i < idpoints_table[identnr].size(); i++) {
                INDEX_2 pair = idpoints_table[identnr][i];
                identmap[pair.I1()] = pair.I2();
                if (symmetric)
                    identmap[pair.I2()] = pair.I1();
            }
        } else {
            MESHIT_LOG_DEBUG("getmap, identnr = " << identnr);
            for (size_t i = 0; i < identifiedpoints_nr->GetNBags(); i++) {
                for (size_t j = 0; j < identifiedpoints_nr->GetBagSize(i); j++) {
                    INDEX_3 i3;
                    int dummy;
                    identifiedpoints_nr->GetData(i, j, i3, dummy);

                    if (i3.I3() == identnr || !identnr) {
                        identmap[i3.I1() - 1] = i3.I2();
                        if (symmetric)
                            identmap[i3.I2() - 1] = i3.I1();
                    }
                }
            }
        }

    }

    void Identifications::GetPairs(int identnr, Array<INDEX_2>& identpairs) const
    {
        identpairs.resize(0);

        if (identnr == 0) {
            for (size_t i = 0; i < identifiedpoints->GetNBags(); i++) {
                for (size_t j = 0; j < identifiedpoints->GetBagSize(i); j++) {
                    INDEX_2 i2;
                    int nr;
                    identifiedpoints->GetData(i, j, i2, nr);
                    identpairs.push_back(i2);
                }
            }
        }
        else {
            for (size_t i = 0; i < identifiedpoints_nr->GetNBags(); i++) {
                for (size_t j = 0; j < identifiedpoints_nr->GetBagSize(i); j++) {
                    INDEX_3 i3;
                    int dummy;
                    identifiedpoints_nr->GetData(i, j, i3, dummy);

                    if (i3.I3() == identnr)
                        identpairs.push_back(INDEX_2(i3.I1(), i3.I2()));
                }
            }
        }
    }

    void Identifications::SetMaxPointNr(int maxpnum)
    {
        for (size_t i = 0; i < identifiedpoints->GetNBags(); i++) {
            for (size_t j = 0; j < identifiedpoints->GetBagSize(i); j++) {
                INDEX_2 i2;
                int nr;
                identifiedpoints->GetData(i, j, i2, nr);

                if (i2.I1() > maxpnum || i2.I2() > maxpnum) {
                    i2.I1() = i2.I2() = -1;
                    identifiedpoints->SetData(i + 1, j + 1, i2, -1);
                }
            }
        }
    }

    MeshingParameters::MeshingParameters()
    {
        optimize2d = "smsmsmSmSmSm";
        optsteps2d = 3;
        opterrpow = 2;
        blockfill = 1;
        filldist = 0.1;
        safety = 5;
        relinnersafety = 3;
        uselocalh = 1;
        grading = -1.0;
        delaunay = 1;
        maxh = 1e10;
        minh = 0;
        meshsizefilename = NULL;
        startinsurface = 0;
        checkoverlap = 1;
        checkoverlappingboundary = 1;
        checkchartboundary = 1;
        curvaturesafety = 2;
        segmentsperedge = 1;
        parthread = 0;

        elsizeweight = 0.2;
        giveuptol2d = 200;
        giveuptol = 10;
        maxoutersteps = 10;
        starshapeclass = 5;
        baseelnp = 0;
        sloppy = 1;
        n_steps = 0;

        badellimit = 175;
        secondorder = 0;
    }

    void MeshingParameters::Print(std::ostream& ost) const
    {
        ost << "Meshing parameters: " << std::endl
        << " optimize2d = " << optimize2d << std::endl
        << " optsteps2d = " << optsteps2d << std::endl
        << " opterrpow = " << opterrpow << std::endl
        << " blockfill = " << blockfill << std::endl
        << " filldist = " << filldist << std::endl
        << " safety = " << safety << std::endl
        << " relinnersafety = " << relinnersafety << std::endl
        << " uselocalh = " << uselocalh << std::endl
        << " grading = " << grading << std::endl
        << " delaunay = " << delaunay << std::endl
        << " maxh = " << maxh << std::endl;
        if (meshsizefilename)
            ost << " meshsizefilename = " << meshsizefilename << std::endl;
        else
            ost << " meshsizefilename = NULL" << std::endl;
        ost << " startinsurface = " << startinsurface << std::endl
        << " checkoverlap = " << checkoverlap << std::endl
        << " checkchartboundary = " << checkchartboundary << std::endl
        << " curvaturesafety = " << curvaturesafety << std::endl
        << " segmentsperedge = " << segmentsperedge << std::endl
        << " parthread = " << parthread << std::endl
        << " elsizeweight = " << elsizeweight << std::endl
        << " giveuptol2d = " << giveuptol2d << std::endl
        << " giveuptol = " << giveuptol << std::endl
        << " maxoutersteps = " << maxoutersteps << std::endl
        << " starshapeclass = " << starshapeclass << std::endl
        << " baseelnp        = " << baseelnp << std::endl
        << " sloppy = " << sloppy << std::endl
        << " badellimit = " << badellimit << std::endl
        << " secondorder = " << secondorder << std::endl
        << " elementorder = " << elementorder << std::endl;
    }

    void MeshingParameters::CopyFrom(const MeshingParameters& other)
    {
        optimize2d = other.optimize2d;
        optsteps2d = other.optsteps2d;
        opterrpow = other.opterrpow;
        blockfill = other.blockfill;
        filldist = other.filldist;
        safety = other.safety;
        relinnersafety = other.relinnersafety;
        uselocalh = other.uselocalh;
        grading = other.grading;
        delaunay = other.delaunay;
        maxh = other.maxh;
        // strcpy(const_cast<char*>(meshsizefilename), other.meshsizefilename);
        // const_cast<char*>(meshsizefilename) = other.meshsizefilename; //???
        startinsurface = other.startinsurface;
        checkoverlap = other.checkoverlap;
        checkoverlappingboundary = other.checkoverlappingboundary;
        checkchartboundary = other.checkchartboundary;
        curvaturesafety = other.curvaturesafety;
        segmentsperedge = other.segmentsperedge;
        parthread = other.parthread;
        elsizeweight = other.elsizeweight;
        giveuptol2d = other.giveuptol2d;
        giveuptol = other.giveuptol;
        maxoutersteps = other.maxoutersteps;
        starshapeclass = other.starshapeclass;
        baseelnp = other.baseelnp;
        sloppy = other.sloppy;
        badellimit = other.badellimit;
        secondorder = other.secondorder;
        elementorder = other.elementorder;
    }

    DebugParameters::DebugParameters()
    {
        haltsuccess = 0;
        haltnosuccess = 0;
        haltlargequalclass = 0;
    };
}

