#include <stdexcept>

#include "mesh_class.hpp"
#include "mesh_types.hpp"

namespace meshit {

Segment::Segment()
{
    pnums[0] = MeshPoint::undefined;
    pnums[1] = MeshPoint::undefined;
    edge_id = -1;
    face_left = static_cast<size_t>(-1);
    face_right = static_cast<size_t>(-1);
}

Segment::Segment(const Segment& other)
    : edge_id(other.edge_id),
      face_left(other.face_left),
      face_right(other.face_right)
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
        face_left = other.face_left;
        face_right = other.face_right;
    }
    return *this;
}

std::ostream& operator<<(std::ostream& s, const Segment& seg)
{
    s << seg[0] << " - " << seg[1] << " domin = " << seg.face_left << ", domout = " << seg.face_right
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
    curvature_safety = other.curvature_safety;
    segments_per_edge = other.segments_per_edge;
    giveup_tol2d = other.giveup_tol2d;
    n_steps = other.n_steps;
}

}  // namespace meshit
