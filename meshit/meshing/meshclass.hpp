#ifndef MESHCLASS_HPP
#define MESHCLASS_HPP

/**************************************************************************/
/* File:   meshclass.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Nov. 99                                                    */
/**************************************************************************/

#include <iostream>
#include <string>
#include "meshtype.hpp"
#include "localh.hpp"
#include "../gprim/adtree.hpp"
#include "../general/bitarray.hpp"
#include "../general/symbolta.hpp"
#include "../gprim/geomops.hpp"

#include "topology.hpp"

/*
  The mesh class
 */

namespace meshit
{
    enum resthtype
    {
        RESTRICTH_FACE, RESTRICTH_EDGE,
        RESTRICTH_SURFACEELEMENT, RESTRICTH_POINT, RESTRICTH_SEGMENT
    };

    class SplineGeometry;

    class Mesh
    {
     private:
        /// point coordinates
        std::vector<MeshPoint> points;

        /// line-segments at edges
        std::vector<Segment> segments;
        /// surface elements, 2d-inner elements
        std::vector<Element2d> surfelements;
        /// points will be fixed forever
        std::vector<PointIndex> lockedpoints;

        /// surface indices at boundary nodes
        TABLE<int> surfacesonnode;
        /// boundary edges  (1..normal bedge, 2..segment)
        INDEX_2_CLOSED_HASHTABLE<int>* boundaryedges;
        INDEX_2_CLOSED_HASHTABLE<int>* segmentht;
        INDEX_3_CLOSED_HASHTABLE<int>* surfelementht;

        /// faces of rest-solid
        std::vector<Element2d> openelements;
        /// open segmenets for surface meshing
        std::vector<Segment> opensegments;

        /**
           Representation of local mesh-size h
         */
        LocalH* lochfunc;
        double hglob_;
        double hmin_;

        /**
           the face-index of the surface element maps into
           this table.
         */
        std::vector<FaceDescriptor> facedecoding;

        /// sub-domain materials
        std::vector<char*> materials;

        /// Periodic surface, close surface, etc. identifications
        Identifications* ident;

        /// number of vertices
        size_t numvertices;

        /// element -> face, element -> edge etc ...
        class MeshTopology* topology;

     public:
        /// refinement hierarchy
        std::vector<INDEX_2> mlbetweennodes;

     private:
        void BuildBoundaryEdges(void);

     public:
        bool PointContainedIn2DElement(
            const Point3d& p,
            double lami[3],
            const int element,
            bool consider3D = false) const;

     public:
        Mesh();
        ~Mesh();

        Mesh& operator=(const Mesh& mesh2);

        void BuildFromSplineGeometry(SplineGeometry& geometry, MeshingParameters& mp);

        size_t AddPoint(const Point3d& p, POINTTYPE type = INNERPOINT);

        size_t GetNP() const
        {
            return points.size();
        }

        MeshPoint& Point(size_t pi)
        {
            return points[pi];
        }

        const MeshPoint& Point(size_t pi) const
        {
            return points[pi];
        }

        MeshPoint& operator[](size_t pi)
        {
            return points[pi];
        }

        const MeshPoint& operator[](size_t pi) const
        {
            return points[pi];
        }

        void AddSegment(const Segment& s);

        size_t GetNSeg() const
        {
            return segments.size();
        }

        Segment& LineSegment(size_t i)
        {
            return segments[i];
        }

        const Segment& LineSegment(size_t i) const
        {
            return segments[i];
        }

        void AddSurfaceElement(const Element2d& el);

        size_t GetNSE() const
        {
            return surfelements.size();
        }

        Element2d& SurfaceElement(size_t i)
        {
            return surfelements[i];
        }

        const Element2d& SurfaceElement(size_t i) const
        {
            return surfelements[i];
        }

        void RebuildSurfaceElementLists();
        void GetSurfaceElementsOfFace(size_t facenr, std::vector<SurfaceElementIndex>& sei) const;

        void AddLockedPoint(PointIndex pi);
        void ClearLockedPoints();

        const std::vector<PointIndex>& LockedPoints() const
        {
            return lockedpoints;
        }

        /// Returns number of domains
        int GetNDomains() const;

        /// sets internal tables
        void CalcSurfacesOfNode();

        /**
           finds elements without neighbour and
           boundary elements without inner element.
           Results are stored in openelements.
           if dom == 0, all sub-domains, else subdomain dom */
        void FindOpenElements(size_t dom = 0);

        /**
           finds segments without surface element,
           and surface elements without neighbours.
           store in opensegmentsy
         */
        void FindOpenSegments(size_t surfnr = 0);

        size_t GetNOpenSegments()
        {
            return opensegments.size();
        }

        const Segment& GetOpenSegment(int nr)
        {
            return opensegments[nr];
        }

        /**
           Checks overlap of boundary
           return == 1, iff overlap
         */
        int CheckOverlappingBoundary();
        /**
           Checks consistent boundary
           return == 0, everything ok
         */
        int CheckConsistentBoundary() const;

        /**
           finds average h of surface surfnr if surfnr > 0,
           else of all surfaces.
         */
        double AverageH(size_t surfnr = 0) const;
        /// Calculates localh
        void CalcLocalH();
        void SetLocalH(const Point3d& pmin, const Point3d& pmax, double grading);
        void RestrictLocalH(const Point3d& p, double hloc);
        void RestrictLocalHLine(const Point3d& p1, const Point3d& p2,
                                double hloc);
        void RestrictLocalH(resthtype rht, size_t nr, double loc_h);
        double GetH(const Point3d& p) const;
        double GetMinH(const Point3d& pmin, const Point3d& pmax);

        /// Find bounding box
        void GetBox(Point3d& pmin, Point3d& pmax) const;

        size_t GetNOpenElements() const
        {
            return openelements.size();
        }

        const Element2d& OpenElement(size_t i) const
        {
            return openelements[i];
        }

        /// Refines mesh and projects points to true surface
        // void Refine (int levels, const CSGeometry * geom);

        bool IsSegment(PointIndex pi1, PointIndex pi2) const
        {
            INDEX_2 i2(pi1, pi2);
            i2.Sort();
            return segmentht->Used(i2);
        }

        /**
           Remove unused points. etc.
         */
        void Compress();

        void Export(const std::string& filetype, const std::string& filename) const;
        void Export(const std::string& filetype, std::ostream& os) const;
        void Export(const std::string& filename) const;
        void Export(std::ostream& os) const;
        void Save(std::ostream& outfile) const;
        void Load(std::istream& infile);
        void Save(const std::string& filename) const;
        void Load(const std::string& filename);

        void SetMaterial(size_t domnr, const char* mat);

        size_t GetNFD() const
        {
            return facedecoding.size();
        }

        const FaceDescriptor& GetFaceDescriptor(int i) const
        {
            return facedecoding[i - 1];
        }

        FaceDescriptor& GetFaceDescriptor(int i)
        {
            return facedecoding[i - 1];
        }

        /// return periodic, close surface etc. identifications

        Identifications& GetIdentifications()
        {
            return *ident;
        }
        /// return periodic, close surface etc. identifications

        /// find number of vertices
        void ComputeNVertices();

        /// number of vertices (no edge-midpoints)
        size_t GetNV() const {
            return numvertices;
        }

        /// remove edge points
        void SetNP(size_t np);

        void UpdateTopology();

        class CSurfaceArea
        {
            const Mesh& mesh;
            bool valid;
            double area;

         public:
            CSurfaceArea(const Mesh& amesh)
                : mesh(amesh), valid(false) { }

            void Add(const Element2d& sel)
            {
                area += Cross(mesh.Point(sel[1]) - mesh.Point(sel[0]),
                              mesh.Point(sel[2]) - mesh.Point(sel[0])).Length() / 2;
            }

            double Area() const { return area; }

            bool Valid() const
            {
                return valid;
            }
        };

        CSurfaceArea surfarea;

        double SurfaceArea()
        {
            return surfarea.Area();
        }

        void PrintMemInfo(std::ostream& ost) const;

        friend class Meshing3;
    };

    inline std::ostream& operator<<(std::ostream& ost, const Mesh& mesh)
    {
        ost << "mesh: " << std::endl;
        mesh.Save(ost);
        return ost;
    }

}  // namespace meshit

#endif

