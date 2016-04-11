#include <stdexcept>

#include "meshtype.hpp"
#include "meshclass.hpp"

namespace meshit
{
    Segment::Segment()
    {
        pnums[0] = -1;
        pnums[1] = -1;
        pnums[2] = -1;

        edgenr = -1;
        seginfo = 0;

        si = -1;

        domin = static_cast<size_t>(-1);
        domout = static_cast<size_t>(-1);
    }

    Segment::Segment(const Segment& other)
        : edgenr(other.edgenr),
          seginfo(other.seginfo),
          si(other.si),
          domin(other.domin),
          domout(other.domout),
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
            s << " " << el.PointID(j);
        }
        return s;
    }

    FaceDescriptor::FaceDescriptor()
    {
        surfnr = bcprop = 0;
        firstelement = -1;
    }

    FaceDescriptor::FaceDescriptor(const FaceDescriptor& other)
        : surfnr(other.surfnr), bcprop(other.bcprop)
    {
        firstelement = -1;
    }

    FaceDescriptor::FaceDescriptor(size_t surfnri)
    {
        surfnr = surfnri;
        bcprop = surfnri;
        firstelement = -1;
    }

    std::ostream& operator<<(std::ostream& s, const FaceDescriptor& fd)
    {
        s << "surfnr = " << fd.SurfNr() << ", bcprop = " << fd.BCProperty();
        return s;
    }

    MeshingParameters::MeshingParameters()
    {
        optimize2d = "smsmsmSmSmSm";
        optsteps2d = 3;
        grading = -1.0;
        maxh = 1e10;
        minh = 0;
        restrict_segment = false;
        curvature_safety = 2;
        segments_per_edge = 1;

        giveup_tol2d = 200;
        n_steps = 0;
    }

    MeshingParameters::MeshingParameters(const MeshingParameters& other)
    {
        optimize2d = other.optimize2d;
        optsteps2d = other.optsteps2d;
        restrict_segment = other.restrict_segment;
        grading = other.grading;
        minh = other.minh;
        maxh = other.maxh;
        curvature_safety = other.curvature_safety;
        segments_per_edge = other.segments_per_edge;
        giveup_tol2d = other.giveup_tol2d;
        n_steps = other.n_steps;
    }

    DebugParameters::DebugParameters()
    {
        haltnosuccess = 0;
        haltlargequalclass = 0;
    };
}  // namespace meshit

