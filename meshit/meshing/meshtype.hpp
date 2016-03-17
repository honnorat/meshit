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

namespace meshit {

    typedef int ELEMENT_EDGE[2]; // initial point, end point

#define ELEMENT2D_MAXPOINTS 3

    enum POINTTYPE
    {
        FIXEDPOINT = 1, EDGEPOINT = 2, SURFACEPOINT = 3, INNERPOINT = 4
    };

    extern int NextTimeStamp();

    class EdgePointGeomInfo
    {
     public:
        int edgenr;
        int body;     // for ACIS
        double dist;  // for 2d meshing
        double u, v;  // for OCC Meshing

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
    };

    typedef int PointIndex;
    typedef int SurfaceElementIndex;

    /**
       Point in the mesh.
       Contains layer (a new feature in 4.3 for overlapping meshes.
     */
    class MeshPoint : public Point3d
    {
        int layer;
        double singular; // singular factor for hp-refinement
        POINTTYPE type;

     public:

        MeshPoint() { }

        MeshPoint(const Point3d& ap, int alayer = 1, POINTTYPE apt = INNERPOINT)
                : Point3d(ap), layer(alayer), singular(0.), type(apt) { }

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

        double Singularity() const
        {
            return singular;
        }

        void Singularity(double s)
        {
            singular = s;
        }

    };

    /**
       Triangle element for surface mesh generation.
     */
    class Element2d
    {
        /// point numbers
        PointIndex pnum[3];

        /// surface nr
        int index : 16;

        // marked for refinement
        bool deleted : 1; // element is deleted

        // element visible

        /// a linked list for all segments in the same face
        SurfaceElementIndex next;

     public:
        Element2d();

        PointIndex& operator[](size_t i)
        {
            return pnum[i];
        }

        const PointIndex& operator[](size_t i) const
        {
            return pnum[i];
        }

        FlatArray<const PointIndex> PNums() const
        {
            return FlatArray<const PointIndex>(3, &pnum[0]);
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

        int GetIndex() const
        {
            return index;
        };

        void GetBox(const Array<MeshPoint>& points, Box3d& box) const;
        /// invert orientation
        inline void Invert();
        /// first point number is smallest
        inline void NormalizeNumbering();

        // friend ostream & operator<<(ostream  & s, const Element2d & el);
        friend class Mesh;

        void GetTransformation(class DenseMatrix& pmat, class DenseMatrix& trans) const;

        void ComputeIntegrationPointData() const;

        double CalcJacobianBadness(const Array<MeshPoint>& points) const;

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

    class IntegrationPointData
    {
     public:
        Point3d p;
        double weight;
        Vector shape;
        DenseMatrix dshape;
    };

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
        double singedge_left;
        double singedge_right;

        /// 0.. not first segment of segs, 1..first of class, 2..first of class, inverse
        unsigned int seginfo : 2;

        /// surface decoding index
        int si;
        /// domain number inner side
        int domin;
        /// domain number outer side
        int domout;
        /// top-level object number of surface
        int tlosurf;

        /// surfaces describing edge
        int surfnr1, surfnr2;
        EdgePointGeomInfo epgeominfo[2];
        // int pmid; // for second order
        int meshdocval;

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
        int surfnr;
        /// domain nr inside
        int domin;
        /// domain nr outside
        int domout;
        /// top level object number of surface
        int tlosurf;
        /// boundary condition property
        int bcprop;

        /// root of linked list
        SurfaceElementIndex firstelement;

        double domin_singular;
        double domout_singular;

     public:
        FaceDescriptor();
        FaceDescriptor(int surfnri, int domini, int domouti, int tlosurfi);
        FaceDescriptor(const Segment& seg);
        FaceDescriptor(const FaceDescriptor& other);

        ~FaceDescriptor() { }

        int SurfNr() const
        {
            return surfnr;
        }

        int DomainIn() const
        {
            return domin;
        }

        int DomainOut() const
        {
            return domout;
        }

        int TLOSurface() const
        {
            return tlosurf;
        }

        int BCProperty() const
        {
            return bcprop;
        }

        double DomainInSingular() const
        {
            return domin_singular;
        }

        double DomainOutSingular() const
        {
            return domout_singular;
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
        /// power of error (to approximate max err optimization)
        double opterrpow;
        /// do block filling ?  
        int blockfill;
        /// block filling up to distance
        double filldist;
        /// radius of local environment (times h)
        double safety;
        /// radius of active environment (times h)
        double relinnersafety;
        /// use local h ?
        int uselocalh;
        /// grading for local h
        double grading;
        /// use delaunay meshing
        int delaunay;
        /// maximal mesh size
        double maxh;
        /// minimal mesh size
        double minh;
        /// file for meshsize
        const char* meshsizefilename;
        /// start surfacemeshing from everywhere in surface
        int startinsurface;
        /// check overlapping surfaces (debug)
        int checkoverlap;
        /// check overlapping surface mesh before volume meshing
        int checkoverlappingboundary;
        /// check chart boundary (sometimes too restrictive)
        int checkchartboundary;
        /// safty factor for curvatures (elemetns per radius)
        double curvaturesafety;
        /// minimal number of segments per edge
        double segmentsperedge;
        /// use parallel threads
        int parthread;
        /// weight of element size w.r.t element shape
        double elsizeweight;
        /// init with default values

        /// from mp3:
        /// give up quality class, 2d meshing
        int giveuptol2d;
        /// give up quality class, 3d meshing
        int giveuptol;
        /// maximal outer steps
        int maxoutersteps;
        /// class starting star-shape filling
        int starshapeclass;
        /// if non-zero, baseelement must have baseelnp points
        int baseelnp;
        /// quality tolerances are handled less careful
        int sloppy;

        /// limit for max element angle (150-180)
        double badellimit;

        int secondorder;
        /// high order element curvature
        int elementorder;

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
            }
            else {
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

        enum ID_TYPE
        {
            UNDEFINED = 1, PERIODIC = 2
        };

     private:
        Mesh& mesh;

        /// identify points (thin layers, periodic b.c.)  
        INDEX_2_HASHTABLE<int>* identifiedpoints;

        /// the same, with info about the id-nr
        INDEX_3_HASHTABLE<int>* identifiedpoints_nr;

        /// sorted by identification nr
        TABLE<INDEX_2> idpoints_table;

        Array<ID_TYPE> type;

        /// number of identifications (or, actually used identifications ?)
        int maxidentnr;

     public:
        Identifications(Mesh& amesh);
        ~Identifications();

        /*
          Identify points pi1 and pi2, due to
          identification nr identnr
         */
        void Add(PointIndex pi1, PointIndex pi2, int identnr);

        int Get(PointIndex pi1, PointIndex pi2) const;
        bool Get(PointIndex pi1, PointIndex pi2, int identnr) const;

        bool UsedSymmetric(PointIndex pi1, PointIndex pi2)
        {
            return identifiedpoints->Used(INDEX_2(pi1, pi2)) ||
                   identifiedpoints->Used(INDEX_2(pi2, pi1));
        }

        void GetMap(size_t identnr, Array<int>& identmap, bool symmetric = false) const;

        ID_TYPE GetType(size_t identnr) const
        {
            if (identnr < type.size()) {
                return type[identnr];
            } else {
                return UNDEFINED;
            }
        }

        void SetType(size_t identnr, ID_TYPE t)
        {
            while (type.size() <= identnr) {
                type.push_back(UNDEFINED);
            }
            type[identnr] = t;
        }

        void GetPairs(int identnr, Array<INDEX_2>& identpairs) const;

        int GetMaxNr() const
        {
            return maxidentnr;
        }

        /// remove secondorder
        void SetMaxPointNr(int maxpnum);

    };

}  // namelist meshit

#endif

