#include <stdexcept>

#include "meshtype.hpp"
#include "meshclass.hpp"

namespace meshit
{
    Segment::Segment()
    {
        pnums[0] = -1;
        pnums[1] = -1;
        edge_id = -1;
        dom_left = static_cast<size_t>(-1);
        dom_right = static_cast<size_t>(-1);
    }

    Segment::Segment(const Segment& other)
        : edge_id(other.edge_id),
          dom_left(other.dom_left),
          dom_right(other.dom_right)
    {
        pnums[0] = other.pnums[0];
        pnums[1] = other.pnums[1];
    }

    Segment& Segment::operator=(const Segment& other)
    {
        if (&other != this) {
            pnums[0] = other.pnums[0];
            pnums[1] = other.pnums[1];
            edge_id = other.edge_id;
            dom_left = other.dom_left;
            dom_right = other.dom_right;
        }
        return *this;
    }

    std::ostream& operator<<(std::ostream& s, const Segment& seg)
    {
        s << seg[0] << " - " << seg[1] << " domin = " << seg.dom_left << ", domout = " << seg.dom_right
        << " si = " << seg.edge_id;
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

    std::ostream& operator<<(std::ostream& os, const FaceDescriptor& fd)
    {
        os << "surf_id = " << fd.face_id();
        return os;
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

}  // namespace meshit

