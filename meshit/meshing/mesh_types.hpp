#ifndef MESHIT_MESH_TYPES_HPP
#define MESHIT_MESH_TYPES_HPP
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

#include <iostream>

#include "../general/index.hpp"
#include "../general/logging.hpp"
#include "../gprim/geom3d.hpp"
#include "../linalg/densemat.hpp"

namespace meshit {

enum PointType { EDGE_POINT = 2, INNER_POINT = 3 };

typedef PointType point_type_t;

typedef GenericIndex PointIndex;
typedef GenericIndex EdgeIndex;
typedef GenericIndex ElementIndex;
typedef GenericIndex DomainIndex;

template<typename T>
struct CONST
{
    static constexpr T undefined = static_cast<T>(-1);
};

/**
   Point in the mesh.
 */
class MeshPoint : public Point2d
{
 public:
    MeshPoint() { }

    explicit MeshPoint(const Point2d& p, point_type_t type = INNER_POINT)
        : Point2d{p}, type_{type} { }

    explicit MeshPoint(double px, double py, point_type_t type = INNER_POINT)
        : Point2d{px, py}, type_{type} { }

    point_type_t Type() const { return type_; }
    void SetType(point_type_t at) { type_ = at; }

 protected:
    point_type_t type_;
};

/**
   Triangle element for surface mesh generation.
 */
class Element2d
{
 public:
    Element2d()
        : point_ids_{undef_pid, undef_pid, undef_pid}, face_id_{0}, deleted_{false} { }

    PointIndex& operator[](size_t i) { return point_ids_[i]; }
    const PointIndex& operator[](size_t i) const { return point_ids_[i]; }

    PointIndex& PointID(size_t i) { return point_ids_[i]; }
    const PointIndex& PointID(size_t i) const { return point_ids_[i]; }

    void SetFaceID(DomainIndex fid) { face_id_ = fid; }
    DomainIndex FaceID() const { return face_id_; }

    void Invert() { std::swap(point_ids_[1], point_ids_[2]); }  // invert orientation
    void NormalizeNumbering();

    friend class Mesh;

    void Delete()
    {
        deleted_ = true;
        point_ids_[0] = point_ids_[1] = point_ids_[2] = undef_pid;
    }

    bool IsDeleted() const { return deleted_; }

    bool operator==(const Element2d& el2) const;

 protected:
    static constexpr PointIndex undef_pid = CONST<PointIndex>::undefined;

    PointIndex point_ids_[3];
    DomainIndex face_id_;
    ElementIndex next_;  // a linked list for all segments in the same face
    bool deleted_;       // element is deleted
};

std::ostream& operator<<(std::ostream& s, const Element2d& el);

/**
   Edge segment.
 */
class Segment
{
 public:
    Segment();
    Segment(const Segment& other);

    ~Segment() { }

    PointIndex pnums[2];  // p1, p2

    EdgeIndex edge_id;
    DomainIndex face_left;   // domain number left side
    DomainIndex face_right;  // domain number right side

 public:
    Segment& operator=(const Segment& other);

    PointIndex& operator[](size_t i) { return pnums[i]; }
    const PointIndex& operator[](size_t i) const { return pnums[i]; }
};

std::ostream& operator<<(std::ostream& s, const Segment& seg);

class FaceDescriptor
{
 public:
    FaceDescriptor()
        : index_{CONST<DomainIndex>::undefined},
          first_element{CONST<ElementIndex>::undefined} { }

    explicit FaceDescriptor(DomainIndex face_index)
        : index_{face_index},
          first_element{CONST<ElementIndex>::undefined} { }

    DomainIndex face_id() const { return index_; }

 protected:
    DomainIndex index_;
    ElementIndex first_element;  // root of linked list

    friend class Mesh;
};

std::ostream& operator<<(std::ostream& os, const FaceDescriptor& fd);


class MeshingParameters
{
 public:
    MeshingParameters();
    MeshingParameters(const MeshingParameters& other);

 public:
    /**
       2d optimization strategy:
       // s .. Swap, opt 6 lines/node
       // S .. Swap, optimal elements
       // m .. move nodes
       // c .. combine
       // p .. print for debug
     **/
    const char* optimize2d;
    int optsteps2d;           // number of 2d optimization steps
    double grading;           // grading for local h
    double maxh;              // maximal mesh size
    double minh;              // minimal mesh size
    double curvature_safety;  // safty factor for curvatures (elemetns per radius)
    int segments_per_edge;    // minimal number of segments per edge
    uint32_t giveup_tol2d;    // give up quality class

    int n_steps;
};

inline void Element2d::NormalizeNumbering()
{
    if (point_ids_[0] < point_ids_[1] && point_ids_[0] < point_ids_[2]) {
        return;
    } else {
        if (point_ids_[1] < point_ids_[2]) {
            PointIndex pi1 = point_ids_[1];
            point_ids_[1] = point_ids_[2];
            point_ids_[2] = point_ids_[0];
            point_ids_[0] = pi1;
        } else {
            PointIndex pi1 = point_ids_[2];
            point_ids_[2] = point_ids_[1];
            point_ids_[1] = point_ids_[0];
            point_ids_[0] = pi1;
        }
    }
}

}  // namespace meshit

#endif
