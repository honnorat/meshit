#include <stdexcept>

#include "meshtype.hpp"
#include "meshclass.hpp"

namespace meshit
{
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
    }

    Segment::Segment(const Segment& other)
        : edgenr(other.edgenr),
          seginfo(other.seginfo),
          si(other.si),
          domin(other.domin),
          domout(other.domout),
          tlosurf(other.tlosurf),
          surfnr1(other.surfnr1),
          surfnr2(other.surfnr2),
          epgeominfo(),
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

    bool Element2d::operator==(const Element2d& el2) const
    {
        bool retval = true;
        for (int i = 0; retval && i < 3; i++) {
            retval = (el2[i] == (*this)[i]);
        }

        return retval;
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

    FaceDescriptor::FaceDescriptor(size_t surfnri, size_t domini, size_t domouti, int tlosurfi)
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
    }

    Identifications::~Identifications()
    {
        delete identifiedpoints;
        delete identifiedpoints_nr;
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

    void Identifications::SetMaxPointNr(int maxpnum)
    {
        for (size_t i = 0; i < identifiedpoints->GetNBags(); i++) {
            for (size_t j = 0; j < identifiedpoints->GetBagSize(i); j++) {
                INDEX_2 i2;
                int nr;
                identifiedpoints->GetData(i, j, i2, nr);

                if (i2.I1() > maxpnum || i2.I2() > maxpnum) {
                    i2.I1() = i2.I2() = -1;
                    identifiedpoints->SetData(i, j, i2, -1);
                }
            }
        }
    }

    MeshingParameters::MeshingParameters()
    {
        optimize2d = "smsmsmSmSmSm";
        optsteps2d = 3;
        grading = -1.0;
        maxh = 1e10;
        minh = 0;
        check_overlap = 1;
        check_chart_boundary = 1;
        curvature_safety = 2;
        segments_per_edge = 1;

        giveup_tol2d = 200;
        n_steps = 0;
    }

    MeshingParameters::MeshingParameters(const MeshingParameters& other)
    {
        optimize2d = other.optimize2d;
        optsteps2d = other.optsteps2d;
        grading = other.grading;
        minh = other.minh;
        maxh = other.maxh;
        check_overlap = other.check_overlap;
        check_chart_boundary = other.check_chart_boundary;
        curvature_safety = other.curvature_safety;
        segments_per_edge = other.segments_per_edge;
        giveup_tol2d = other.giveup_tol2d;
        n_steps = other.n_steps;
    }

    void MeshingParameters::Print(std::ostream& ost) const
    {
        ost << "Meshing parameters: " << std::endl
        << "  optimize2d = " << optimize2d << std::endl
        << "  optsteps2d = " << optsteps2d << std::endl
        << "  grading = " << grading << std::endl
        << "  maxh = " << maxh << std::endl
        << "  check_overlap = " << check_overlap << std::endl
        << "  check_chart_boundary = " << check_chart_boundary << std::endl
        << "  curvature_safety = " << curvature_safety << std::endl
        << "  segments_per_edge = " << segments_per_edge << std::endl
        << "  giveup_tol2d = " << giveup_tol2d << std::endl;
    }

    DebugParameters::DebugParameters()
    {
        haltsuccess = 0;
        haltnosuccess = 0;
        haltlargequalclass = 0;
    };
}  // namespace meshit

