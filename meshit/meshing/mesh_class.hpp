#ifndef MESHIT_MESH_CLASS_HPP
#define MESHIT_MESH_CLASS_HPP
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
#include <string>

#include "../general/table.hpp"
#include "../gprim/adtree.hpp"
#include "../gprim/geomobjects.hpp"
#include "localh.hpp"
#include "mesh_types.hpp"

namespace meshit {

class SplineGeometry;

class Mesh
{
 private:
    std::vector<MeshPoint> points;    // point coordinates
    std::vector<Segment> segments;    // line-segments at edges
    std::vector<Element2d> elements;  // 2d-inner elements

    IndexPair_map<int> segment_ht;  // boundary edges  (1..normal bedge, 2..segment)

    // Representation of local mesh-size h
    LocalH* loc_h_func;
    double hglob_;
    double hmin_;

    /**
       the face-index of the surface element maps into
       this table.
     */
    std::vector<FaceDescriptor> faces;

    /// sub-domain materials
    std::vector<char*> materials;

 public:
    Mesh();
    ~Mesh();

    Mesh& operator=(const Mesh& mesh2);

    void BuildFromSplineGeometry(SplineGeometry& geometry, MeshingParameters& mp, bool geometry_only = false);

    size_t AddPoint(const Point2d& p, point_type_t type = INNER_POINT);
    void AddSegment(const Segment& s);
    void AddElement(const Element2d& el);

    MeshPoint& Point(size_t pi) { return points[pi]; }
    Element2d& Element(size_t i) { return elements[i]; }

    const MeshPoint& Point(size_t pi) const { return points[pi]; }
    const Segment& LineSegment(size_t i) const { return segments[i]; }
    const Element2d& Element(size_t i) const { return elements[i]; }

    size_t GetNbPoints() const { return points.size(); }
    size_t GetNbSegments() const { return segments.size(); }
    size_t GetNbElements() const { return elements.size(); }
    size_t GetNbFaces() const { return faces.size(); }

    void RebuildSurfaceElementLists();
    void GetElementsOfFace(DomainIndex facenr, std::vector<ElementIndex>& sei) const;

    /// sets internal tables
    void IndexBoundaryEdges();

    int CheckSurface();
    int CheckOverlappingBoundary();  // Checks overlap of boundary, return == 1, iff overlap
    void PrintQuality();

    // finds average h of surface surfnr if surfnr > 0, else of all surfaces.
    double AverageH(DomainIndex surf_id = 0) const;

    void CalcLocalH();
    void SetLocalH(const Box2d& bbox, double grading);
    void RestrictLocalH(const Point2d& p, double hloc);
    void RestrictLocalHLine(const Point2d& p1, const Point2d& p2, double hloc);
    double GetH(const Point2d& p) const;
    double GetMinH(const Point2d& pmin, const Point2d& pmax);

    /// Find bounding box
    void GetBox(Point2d& pmin, Point2d& pmax) const;

    bool IsSegment(const IndexPair& e) const { return segment_ht.count(e) == 1; }

    void Compress();  // Remove unused points. etc.
    size_t ComputeNVertices();

    void Export(const std::string& filetype, const std::string& filename) const;
    void Export(const std::string& filetype, std::ostream& os) const;
    void Export(const std::string& filename) const;
    void Export(std::ostream& os) const;
    void Save(std::ostream& outfile) const;
    void Load(std::istream& infile);
    void Save(const std::string& filename) const;
    void Load(const std::string& filename);

    void Refine();
    void Optimize2d(MeshingParameters& mp);

    void SetMaterial(size_t domnr, const char* mat);

    int GetDomainNumber(const Element2d& el) const { return GetDomainNumber(el.FaceID()); }
    int GetDomainNumber(DomainIndex face_id) const { return faces[face_id - 1].face_id() + 100; }

    class CSurfaceArea
    {
     public:
        explicit CSurfaceArea(const Mesh& amesh)
            : mesh(amesh), valid(false) { }

        void Add(const Element2d& sel)
        {
            const MeshPoint& p0 = mesh.Point(sel[0]);
            const MeshPoint& p1 = mesh.Point(sel[1]);
            const MeshPoint& p2 = mesh.Point(sel[2]);

            double v1_x = p1.X() - p0.X();  // v1 = p1 - p0
            double v1_y = p1.Y() - p0.Y();
            double v2_x = p2.X() - p0.X();  // v2 = p2 - p0
            double v2_y = p2.Y() - p0.Y();

            // (1/2) * || (p1-p0) x (p2-p0) ||
            area += 0.5 * fabs(v1_x * v2_y - v1_y * v2_x);
        }

        double Area() const { return area; }
        bool Valid() const { return valid; }

     protected:
        const Mesh& mesh;
        bool valid;
        double area;
    };

    CSurfaceArea surfarea;

    double SurfaceArea() { return surfarea.Area(); }

    void PrintMemInfo(std::ostream& ost) const;
};

inline std::ostream& operator<<(std::ostream& ost, const Mesh& mesh)
{
    ost << "mesh: " << std::endl;
    mesh.Save(ost);
    return ost;
}

}  // namespace meshit

#endif
