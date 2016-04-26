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

#include <stdexcept>

#include "mesh_types.hpp"
#include "mesh_class.hpp"

namespace meshit {

template<typename T> constexpr T CONST<T>::undefined;

Segment::Segment()
{
    pnums[0] = CONST<PointIndex>::undefined;
    pnums[1] = CONST<PointIndex>::undefined;
    edge_id = CONST<EdgeIndex>::undefined;
    face_left = CONST<DomainIndex>::undefined;
    face_right = CONST<DomainIndex>::undefined;
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
    s << seg[0] << " - " << seg[1];
    s << " domin = " << seg.face_left << ",";
    s << " domout = " << seg.face_right;
    s << " si = " << seg.edge_id;
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
