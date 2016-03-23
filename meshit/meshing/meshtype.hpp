#ifndef MESHTYPE_HPP
#define MESHTYPE_HPP

/**************************************************************************/
/* File:   meshtype.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include <iostream>

#include "../meshit.hpp"
#include "msghandler.hpp"
#include "../linalg/densemat.hpp"
#include "../gprim/geom3d.hpp"
#include "../general/hashtabl.hpp"

namespace meshit
{
    typedef int ELEMENT_EDGE[2];  // initial point, end point

    enum POINTTYPE
    {
        FIXEDPOINT = 1, EDGEPOINT = 2, SURFACEPOINT = 3, INNERPOINT = 4
    };

    class EdgePointGeomInfo
    {
     public:
        EdgePointGeomInfo()
            : edgenr(0), body(0), dist(0.0), u(0.0), v(0.0) { }

        EdgePointGeomInfo& operator=(const EdgePointGeomInfo& gi2)
        {
            edgenr = gi2.edgenr;
            body = gi2.body;
            dist = gi2.dist;
            u = gi2.u;
            v = gi2.v;
            return *this;
        }

     public:
        int edgenr;
        int body;     // for ACIS
        double dist;  // for 2d meshing
        double u, v;  // for OCC Meshing
    };

    typedef int PointIndex;
    typedef int SurfaceElementIndex;

    /**
       Point in the mesh.
       Contains layer (a new feature in 4.3 for overlapping meshes.
     */
    class MeshPoint : public Point3d
    {
     public:
        MeshPoint() { }

        explicit MeshPoint(const Point3d& ap, int alayer = 1, POINTTYPE apt = INNERPOINT)
            : Point3d(ap), layer(alayer), type(apt) { }

        int GetLayer() const
        {
            return layer;
        }

        POINTTYPE Type() const
        {
            return type;
        }

        void SetType(POINTTYPE at)
        {
            type = at;
        }

     protected:
        int layer;
        POINTTYPE type;
    };

    /**
       Triangle element for surface mesh generation.
     */
    class Element2d
    {
        /// surface nr
        size_t index;

        /// point numbers
        PointIndex pnum[3];

        // marked for refinement
        bool deleted;  // element is deleted

        // element visible

        /// a linked list for all segments in the same face
        SurfaceElementIndex next;

     public:
        Element2d()
            : index{0}, pnum{0, 0, 0}, deleted{false} { }

        PointIndex& operator[](size_t i)
        {
            return pnum[i];
        }

        const PointIndex& operator[](size_t i) const
        {
            return pnum[i];
        }

        PointIndex& PNum(size_t i)
        {
            return pnum[i - 1];
        }

        const PointIndex& PNum(size_t i) const
        {
            return pnum[i - 1];
        }

        PointIndex& PNumMod(size_t i)
        {
            return pnum[(i - 1) % 3];
        }

        const PointIndex& PNumMod(size_t i) const
        {
            return pnum[(i - 1) % 3];
        }

        void SetIndex(int si)
        {
            index = si;
        }

        size_t GetIndex() const
        {
            return index;
        };

        /// invert orientation
        inline void Invert();
        /// first point number is smallest
        inline void NormalizeNumbering();

        // friend ostream & operator<<(ostream  & s, const Element2d & el);
        friend class Mesh;

        void Delete()
        {
            deleted = 1;
            pnum[0] = pnum[1] = pnum[2] = -1;
        }

        bool IsDeleted() const
        {
            return deleted;
        }

        // Philippose - 08 August 2010
        // Access functions for the new property: visible

        bool operator==(const Element2d& el2) const;
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

        PointIndex pnums[3];  // p1, p2, pmid

        int edgenr;

        /// 0.. not first segment of segs, 1..first of class, 2..first of class, inverse
        unsigned int seginfo : 2;

        /// surface decoding index
        int si;
        /// domain number inner side
        size_t domin;
        /// domain number outer side
        size_t domout;
        /// top-level object number of surface
        int tlosurf;

        /// surfaces describing edge
        int surfnr1, surfnr2;
        EdgePointGeomInfo epgeominfo[2];

     public:
        Segment& operator=(const Segment& other);

        int hp_elnr;

        PointIndex& operator[](int i)
        {
            return pnums[i];
        }

        const PointIndex& operator[](int i) const
        {
            return pnums[i];
        }
    };

    std::ostream& operator<<(std::ostream& s, const Segment& seg);

    class FaceDescriptor
    {
        /// which surface, 0 if not available
        size_t surfnr;
        /// domain nr inside
        size_t domin;
        /// domain nr outside
        size_t domout;
        /// top level object number of surface
        int tlosurf;
        /// boundary condition property
        size_t bcprop;

        /// root of linked list
        SurfaceElementIndex firstelement;

     public:
        FaceDescriptor();
        FaceDescriptor(const FaceDescriptor& other);
        FaceDescriptor(size_t surfnri, size_t domini, size_t domouti, int tlosurfi);
        explicit FaceDescriptor(const Segment& seg);

        ~FaceDescriptor() { }

        size_t SurfNr() const
        {
            return surfnr;
        }

        size_t DomainIn() const
        {
            return domin;
        }

        size_t DomainOut() const
        {
            return domout;
        }

        int TLOSurface() const
        {
            return tlosurf;
        }

        size_t BCProperty() const
        {
            return bcprop;
        }

        void SetBCProperty(int bc)
        {
            bcprop = bc;
        }

        friend class Mesh;
    };

    std::ostream& operator<<(std::ostream& s, const FaceDescriptor& fd);

    class MeshingParameters
    {
     public:
        /**
           2d optimization strategy:
           // s .. Swap, opt 6 lines/node
           // S .. Swap, optimal elements
           // m .. move nodes
           // c .. combine
         **/
        const char* optimize2d;
        /// number of 2d optimization steps
        int optsteps2d;
        /// use local h ?
        int uselocalh;
        /// grading for local h
        double grading;
        /// maximal mesh size
        double maxh;
        /// minimal mesh size
        double minh;
        /// check overlapping surfaces (debug)
        int check_overlap;
        /// check chart boundary (sometimes too restrictive)
        int check_chart_boundary;
        /// safty factor for curvatures (elemetns per radius)
        double curvature_safety;
        /// minimal number of segments per edge
        double segments_per_edge;
        /// give up quality class, 2d meshing
        int giveup_tol2d;

        int n_steps;

        MeshingParameters();
        void Print(std::ostream& ost) const;

        void CopyFrom(const MeshingParameters& other);
    };

    class DebugParameters
    {
     public:
        int haltsuccess;
        int haltnosuccess;
        int haltlargequalclass;
        int haltface;
        int haltfacenr;
        DebugParameters();
    };

    inline void Element2d::Invert()
    {
        std::swap(pnum[1], pnum[2]);
    }

    inline void Element2d::NormalizeNumbering()
    {
        if (PNum(1) < PNum(2) && PNum(1) < PNum(3)) {
            return;
        } else {
            if (PNum(2) < PNum(3)) {
                PointIndex pi1 = PNum(2);
                PNum(2) = PNum(3);
                PNum(3) = PNum(1);
                PNum(1) = pi1;
            } else {
                PointIndex pi1 = PNum(3);
                PNum(3) = PNum(2);
                PNum(2) = PNum(1);
                PNum(1) = pi1;
            }
        }
    }

    /**
       Identification of periodic surfaces, close surfaces, etc. 
     */

    class Mesh;

    class Identifications
    {
     public:
        explicit Identifications(Mesh& amesh);
        ~Identifications();

        int Get(PointIndex pi1, PointIndex pi2) const;
        bool Get(PointIndex pi1, PointIndex pi2, int identnr) const;

        bool UsedSymmetric(PointIndex pi1, PointIndex pi2)
        {
            return identifiedpoints->Used(INDEX_2(pi1, pi2)) ||
                   identifiedpoints->Used(INDEX_2(pi2, pi1));
        }

        /// remove secondorder
        void SetMaxPointNr(int maxpnum);

     private:
        Mesh& mesh;

        /// identify points (thin layers, periodic b.c.)
        INDEX_2_HASHTABLE<PointIndex>* identifiedpoints;

        /// the same, with info about the id-nr
        INDEX_3_HASHTABLE<PointIndex>* identifiedpoints_nr;
    };

}  // namespace meshit

#endif

