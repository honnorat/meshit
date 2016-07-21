#ifndef MESHIT_GEOMETRY2D_HPP
#define MESHIT_GEOMETRY2D_HPP
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

#define _USE_MATH_DEFINES 1

#include <string>
#include <vector>

#include "../gprim/geomobjects.hpp"
#include "../gprim/spline.hpp"

namespace meshit {

class SplineGeometry
{
 public:
    SplineGeometry()
        : elto0{0.3} { }

    ~SplineGeometry();

    void GetBoundingBox(Box2d& box) const;
    Box2d GetBoundingBox() const;

    void Load(const std::string& filename);
    void LoadData(std::istream& infile);

    DomainIndex AddFace(const std::string& name, double maxh_f = 1e99);

    size_t AddPoint(const Point2d& point);

    void AddPoints(const std::vector<Point2d>& points);

    void AddSegment(const GeomPoint& p0, const GeomPoint& p1, double hmax, int spline_id,
                    DomainIndex domain_left = 1, DomainIndex domain_right = 0);

    void AddClosedLine(const std::vector<Point2d>& points, double hmax, int spline_id,
                       DomainIndex domain_left = 1, DomainIndex domain_right = 0);

    void AddOpenLine(const std::vector<size_t>& points, double hmax, int spline_id,
                     DomainIndex domain_left, DomainIndex domain_right);

    void AddStructureLine(const std::vector<Point2d>& points, double hmax = 1e99, int bc = 1, DomainIndex domain = 1);

    void AddSpline(const std::vector<Point2d>& points, double hmax, int spline_id = 1,
                   DomainIndex domain_left = 1, DomainIndex domain_right = 0);

    void AddCircle(const Point2d& center, double radius, double hmax, int spline_id = 1,
                   DomainIndex face_left = 1, DomainIndex face_right = 0);

    void AddHole(const std::vector<Point2d>& points, double hmax, int bc, DomainIndex domain = 1);

    void PartitionBoundary(Mesh& mesh2d, MeshingParameters& mp);

    void GetMaterial(DomainIndex domnr, char*& material);
    double GetDomainMaxh(DomainIndex domain_id);
    double GetGrading() { return elto0; }
    void SetGrading(const double grading) { elto0 = grading; }

 protected:
    std::vector<GeomPoint> geompoints;
    std::vector<SplineSeg*> splines;
    std::vector<char*> materials;
    std::vector<double> maxh;
    double elto0;

 private:
    char TestComment(std::istream& infile);
};

class SplineSegmenter
{
 public:
    SplineSegmenter(Mesh& mesh, MeshingParameters& mp, double elto0, const Box2d& bbox)
        : mesh_{mesh}, mp_{mp}, elto0_{elto0}, searchtree_{bbox.PMin(), bbox.PMax()} { }

    void Partition(const std::vector<SplineSeg*>& splines)
    {
        for (size_t i = 0; i < splines.size(); i++) {
            Partition(*splines[i]);
        }
    }

 protected:
    void Partition(const SplineSeg& spline);
    void CalcPartition(const SplineSeg& spline, std::vector<double>& points);

 protected:
    Mesh& mesh_;
    MeshingParameters& mp_;
    double elto0_;
    Point3dTree searchtree_;
};

}  // namespace meshit

#endif
