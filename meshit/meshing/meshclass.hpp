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

namespace meshit {

    enum resthtype
    {
        RESTRICTH_FACE, RESTRICTH_EDGE,
        RESTRICTH_SURFACEELEMENT, RESTRICTH_POINT, RESTRICTH_SEGMENT
    };

    class SplineGeometry2d;

    /// 2d/3d mesh

    class Mesh
    {
     public:
        typedef ::meshit::T_POINTS T_POINTS;
        typedef Array<Element2d> T_SURFELEMENTS;

     private:
        /// point coordinates
        T_POINTS points;

        /// line-segments at edges
        Array<Segment> segments;
        /// surface elements, 2d-inner elements
        T_SURFELEMENTS surfelements;
        /// points will be fixed forever
        Array<PointIndex> lockedpoints;

        /// surface indices at boundary nodes
        TABLE<int, PointIndex::BASE> surfacesonnode;
        /// boundary edges  (1..normal bedge, 2..segment)
        INDEX_2_CLOSED_HASHTABLE<int>* boundaryedges;
        INDEX_2_CLOSED_HASHTABLE<int>* segmentht;
        INDEX_3_CLOSED_HASHTABLE<int>* surfelementht;

        /// faces of rest-solid
        Array<Element2d> openelements;
        /// open segmenets for surface meshing  
        Array<Segment> opensegments;

        /**
           Representation of local mesh-size h
         */
        LocalH* lochfunc;
        double hglob;
        double hmin;
        Array<double> maxhdomain;

        /**
           the face-index of the surface element maps into
           this table.
         */
        Array<FaceDescriptor> facedecoding;

        /// sub-domain materials
        Array<char*> materials;

        /// labels for boundary conditions
        Array<std::string*> bcnames;

        /// Periodic surface, close surface, etc. identifications
        Identifications* ident;

        /// number of vertices (if < 0, use np)
        int numvertices;

        /// geometric search tree for interval intersection search
        Box3dTree* elementsearchtree;
        /// time stamp for tree
        mutable int elementsearchtreets;

        /// element -> face, element -> edge etc ...
        class MeshTopology* topology;

        /// space dimension (2 or 3)
        int dimension;

        /// changed by every minor modification (addpoint, ...)
        int timestamp;

        SYMBOLTABLE<Array<int>*> userdata_int;
        SYMBOLTABLE<Array<double>*> userdata_double;

        mutable Array<Point3d> pointcurves;
        mutable Array<int> pointcurves_startpoint;
        mutable Array<double> pointcurves_red, pointcurves_green, pointcurves_blue;

     private:
        void BuildBoundaryEdges(void);

     public:
        bool PointContainedIn2DElement(
                const Point3d& p,
                double lami[3],
                const int element,
                bool consider3D = false) const;

     public:

        /// refinement hierarchy
        Array<INDEX_2, PointIndex::BASE> mlbetweennodes;

        Mesh();
        ~Mesh();

        Mesh& operator=(const Mesh& mesh2);

        void BuildFromSpline2D(SplineGeometry2d& geometry, MeshingParameters& mp);

        void DeleteMesh();

        void ClearSurfaceElements();

        PointIndex AddPoint(const Point3d& p, int layer = 1);
        PointIndex AddPoint(const Point3d& p, int layer, POINTTYPE type);

        int GetNP() const
        {
            return points.size();
        }

        MeshPoint& Point(int i)
        {
            return points.Elem(i);
        }

        MeshPoint& Point(PointIndex pi)
        {
            return points[pi];
        }

        const MeshPoint& Point(int i) const
        {
            return points.Get(i);
        }

        const MeshPoint& Point(PointIndex pi) const
        {
            return points[pi];
        }

        const MeshPoint& operator[](PointIndex pi) const
        {
            return points[pi];
        }

        MeshPoint& operator[](PointIndex pi)
        {
            return points[pi];
        }

        T_POINTS& Points()
        {
            return points;
        }

        SegmentIndex AddSegment(const Segment& s);

        int GetNSeg() const
        {
            return segments.size();
        }

        Segment& LineSegment(int i)
        {
            return segments[i - 1];
        }

        const Segment& LineSegment(int i) const
        {
            return segments[i - 1];
        }

        Segment& LineSegment(SegmentIndex si)
        {
            return segments[si];
        }

        const Segment& operator[](SegmentIndex si) const
        {
            return segments[si];
        }

        Segment& operator[](SegmentIndex si)
        {
            return segments[si];
        }

        SurfaceElementIndex AddSurfaceElement(const Element2d& el);

        void DeleteSurfaceElement(int eli)
        {
            surfelements.Elem(eli).Delete();
            surfelements.Elem(eli).PNum(1) = -1;
            surfelements.Elem(eli).PNum(2) = -1;
            surfelements.Elem(eli).PNum(3) = -1;
            timestamp = NextTimeStamp();
        }

        void DeleteSurfaceElement(SurfaceElementIndex eli)
        {
            DeleteSurfaceElement(int(eli) + 1);
        }

        int GetNSE() const
        {
            return surfelements.size();
        }

        Element2d& SurfaceElement(int i)
        {
            return surfelements.Elem(i);
        }

        const Element2d& SurfaceElement(int i) const
        {
            return surfelements.Get(i);
        }

        Element2d& SurfaceElement(SurfaceElementIndex i)
        {
            return surfelements[i];
        }

        const Element2d& SurfaceElement(SurfaceElementIndex i) const
        {
            return surfelements[i];
        }

        void RebuildSurfaceElementLists();
        void GetSurfaceElementsOfFace(int facenr, Array<SurfaceElementIndex>& sei) const;

        void AddLockedPoint(PointIndex pi);
        void ClearLockedPoints();

        const Array<PointIndex>& LockedPoints() const
        {
            return lockedpoints;
        }

        /// Returns number of domains
        int GetNDomains() const;

        int GetDimension() const
        {
            return dimension;
        }

        /// sets internal tables
        void CalcSurfacesOfNode();

        /// additional (temporarily) fix points 
        void FixPoints(const BitArray& fixpoints);

        /**
           finds elements without neighbour and
           boundary elements without inner element.
           Results are stored in openelements.
           if dom == 0, all sub-domains, else subdomain dom */
        void FindOpenElements(int dom = 0);

        /**
           finds segments without surface element,
           and surface elements without neighbours.
           store in opensegmentsy
         */
        void FindOpenSegments(int surfnr = 0);

        int GetNOpenSegments()
        {
            return opensegments.size();
        }

        const Segment& GetOpenSegment(int nr)
        {
            return opensegments.Get(nr);
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
        double AverageH(int surfnr = 0) const;
        /// Calculates localh 
        void CalcLocalH();
        void SetLocalH(const Point3d& pmin, const Point3d& pmax, double grading);
        void RestrictLocalH(const Point3d& p, double hloc);
        void RestrictLocalHLine(const Point3d& p1, const Point3d& p2,
                                double hloc);
        /// number of elements per radius
        void CalcLocalHFromSurfaceCurvature(double elperr);
        void CalcLocalHFromPointDistances();
        void RestrictLocalH(resthtype rht, int nr, double loch);
        void LoadLocalMeshSize(const char* meshsizefilename);
        void SetGlobalH(double h);
        void SetMinimalH(double h);
        double MaxHDomain(int dom) const;
        void SetMaxHDomain(const Array<double>& mhd);
        double GetH(const Point3d& p) const;
        double GetMinH(const Point3d& pmin, const Point3d& pmax);

        /// Find bounding box
        void GetBox(Point3d& pmin, Point3d& pmax, int dom = -1) const;

        /// Find bounding box of points of typ ptyp or less
        void GetBox(Point3d& pmin, Point3d& pmax, POINTTYPE ptyp) const;

        int GetNOpenElements() const
        {
            return openelements.size();
        }

        const Element2d& OpenElement(int i) const
        {
            return openelements.Get(i);
        }

        /// are also quads open elements
        bool HasOpenQuads() const;

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

        /// orient surface mesh, for one sub-domain only
        void SurfaceMeshOrientation();

        /// build box-search tree
        void BuildElementSearchTree();

        int AddFaceDescriptor(const FaceDescriptor& fd)
        {
            return facedecoding.push_back(fd);
        }

        void SetMaterial(int domnr, const char* mat);
        const char* GetMaterial(int domnr) const;
        void SetNBCNames(int nbcn);
        void SetBCName(int bcnr, const std::string& abcname);

        const std::string& GetBCName(int bcnr) const;

        void ClearFaceDescriptors()
        {
            facedecoding.resize(0);
        }

        int GetNFD() const
        {
            return facedecoding.size();
        }

        const FaceDescriptor& GetFaceDescriptor(int i) const
        {
            return facedecoding.Get(i);
        }

        FaceDescriptor& GetFaceDescriptor(int i)
        {
            return facedecoding.Elem(i);
        }

        /// return periodic, close surface etc. identifications

        Identifications& GetIdentifications()
        {
            return *ident;
        }
        /// return periodic, close surface etc. identifications

        void InitPointCurve(double red = 1, double green = 0, double blue = 0) const;
        void AddPointCurvePoint(const Point3d& pt) const;
        int GetNumPointCurves(void) const;
        int GetNumPointsOfPointCurve(int curve) const;
        Point3d& GetPointCurvePoint(int curve, int n) const;
        void GetPointCurveColor(int curve, double& red, double& green, double& blue) const;

        /// find number of vertices
        void ComputeNVertices();
        /// number of vertices (no edge-midpoints)
        int GetNV() const;
        /// remove edge points
        void SetNP(int np);

        bool PureTrigMesh(int faceindex = 0) const;

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
                if (sel.GetNP() == 3) {
                    area += Cross(
                            mesh.Point(sel[1]) - mesh.Point(sel[0]),
                            mesh.Point(sel[2]) - mesh.Point(sel[0])).Length() / 2;
                }
                else {
                    area += Cross(
                            Vec3d(mesh.Point(sel.PNum(1)), mesh.Point(sel.PNum(3))),
                            Vec3d(mesh.Point(sel.PNum(1)), mesh.Point(sel.PNum(4)))
                    ).Length() / 2;
                }
            }

            operator double() const
            {
                return area;
            }

            bool Valid() const
            {
                return valid;
            }
        };

        CSurfaceArea surfarea;

        CSurfaceArea& SurfaceArea()
        {
            return surfarea;
        }

        int GetTimeStamp() const
        {
            return timestamp;
        }

        void SetNextTimeStamp()
        {
            timestamp = NextTimeStamp();
        }

        void SetUserData(const char* id, Array<int>& data);
        bool GetUserData(const char* id, Array<int>& data, int shift = 0) const;
        void SetUserData(const char* id, Array<double>& data);
        bool GetUserData(const char* id, Array<double>& data, int shift = 0) const;

        void PrintMemInfo(std::ostream& ost) const;

        friend class Meshing3;
    };

    inline std::ostream& operator<<(std::ostream& ost, const Mesh& mesh)
    {
        ost << "mesh: " << std::endl;
        mesh.Save(ost);
        return ost;
    }

}

#endif

