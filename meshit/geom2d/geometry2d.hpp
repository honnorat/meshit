#ifndef MESHIT_GEOM2D_GEOMETRY2D_HPP
#define MESHIT_GEOM2D_GEOMETRY2D_HPP

/* *************************************************************************/
/* File:   geometry2d.hpp                                                  */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

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

    int AddFace(const std::string& name, double maxh_f = 1e99);
    void AddHole(const std::vector<Point2d>& points, double hmax, int bc, int domain = 1);

    void AddLine(const std::vector<Point2d>& points, double hmax, int spline_id,
                 int domain_left = 1, int domain_right = 0);

    void AddStructureLine(const std::vector<Point2d>& points, double hmax = 1e99, int bc = 1, int domain = 1);

    void AddSpline(const std::vector<Point2d>& points, double hmax, int spline_id = 1,
                   int domain_left = 1, int domain_right = 0);

    void AddCircle(const Point2d& center, double radius, double hmax, int spline_id = 1,
                   int face_left = 1, int face_right = 0);

    void PartitionBoundary(Mesh& mesh2d, MeshingParameters& mp);

    void GetMaterial(size_t domnr, char*& material);
    double GetDomainMaxh(size_t domain_id);
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
        : mesh_{mesh}, mp_{mp}, elto0_{elto0}, searchtree_{bbox.PMin(), bbox.PMax()}
    {
    }

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
