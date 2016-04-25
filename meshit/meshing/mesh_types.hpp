#ifndef MESHTYPE_HPP
#define MESHTYPE_HPP

/**************************************************************************/
/* File:   meshtype.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include <iostream>

#include "../general/logging.hpp"
#include "../gprim/geom3d.hpp"
#include "../linalg/densemat.hpp"

namespace meshit {

enum PointType
{
    EDGE_POINT = 2, INNER_POINT = 3
};
typedef PointType point_type_t;

typedef uint32_t PointIndex;
typedef uint32_t ElementIndex;
typedef uint32_t DomainIndex;

/**
   Point in the mesh.
 */
class MeshPoint : public Point2d
{
 public:
    static constexpr PointIndex undefined = static_cast<PointIndex>(-1);

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
    static constexpr ElementIndex undefined = static_cast<ElementIndex>(-1);

 public:
    Element2d()
        : pnum{undefined, undefined, undefined}, face_id{0}, deleted{false} { }

    PointIndex& operator[](size_t i) { return pnum[i]; }
    const PointIndex& operator[](size_t i) const { return pnum[i]; }

    PointIndex& PointID(size_t i) { return pnum[i]; }
    const PointIndex& PointID(size_t i) const { return pnum[i]; }

    void SetFaceID(DomainIndex fid) { face_id = fid; }
    DomainIndex FaceID() const { return face_id; }

    void Invert() { std::swap(pnum[1], pnum[2]); }  // invert orientation
    void NormalizeNumbering();

    friend class Mesh;

    void Delete()
    {
        deleted = true;
        pnum[0] = pnum[1] = pnum[2] = Element2d::undefined;
    }

    bool IsDeleted() const { return deleted; }

    bool operator==(const Element2d& el2) const;

 protected:
    PointIndex pnum[3];
    DomainIndex face_id;
    ElementIndex next;  // a linked list for all segments in the same face
    bool deleted;              // element is deleted
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

    int edge_id;
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
    static constexpr DomainIndex undefined = static_cast<DomainIndex>(-1);

 public:
    FaceDescriptor()
        : index_{undefined}, first_element{Element2d::undefined} { }

    explicit FaceDescriptor(DomainIndex face_index)
        : index_{face_index}, first_element{Element2d::undefined} { }

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
    int giveup_tol2d;         // give up quality class

    int n_steps;
};

inline void Element2d::NormalizeNumbering()
{
    if (pnum[0] < pnum[1] && pnum[0] < pnum[2]) {
        return;
    } else {
        if (pnum[1] < pnum[2]) {
            PointIndex pi1 = pnum[1];
            pnum[1] = pnum[2];
            pnum[2] = pnum[0];
            pnum[0] = pi1;
        } else {
            PointIndex pi1 = pnum[2];
            pnum[2] = pnum[1];
            pnum[1] = pnum[0];
            pnum[0] = pi1;
        }
    }
}

}  // namespace meshit

#endif
