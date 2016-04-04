#ifndef MESHIT_GEOM2D_GEOMETRY2D_HPP
#define MESHIT_GEOM2D_GEOMETRY2D_HPP

/* *************************************************************************/
/* File:   geometry2d.hpp                                                  */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

#define _USE_MATH_DEFINES 1

#include <vector>
#include <string>

#include "../gprim/spline.hpp"
#include "../gprim/geomobjects.hpp"
#include "../meshing/refine.hpp"

namespace meshit
{
    void Optimize2d(Mesh& mesh, MeshingParameters& mp);

    class SplineGeometry
    {
     protected:
        std::vector<GeomPoint<2> > geompoints;
        std::vector<SplineSeg*> splines;
        std::vector<char*> materials;
        std::vector<double> maxh;
        double elto0;

     public:
        SplineGeometry() : elto0{0.3} { }

        ~SplineGeometry();

        void GetBoundingBox(Box<2>& box) const;
        Box<2> GetBoundingBox() const;

        void Load(const std::string& filename);
        void LoadData(std::istream& infile);

        void AddLine(const std::vector<Point2d>& points,
                     double hmax, int bc,
                     int face_left = 1, int face_right = 0);
        void AddHole(const std::vector<Point2d>& points,
                     double hmax, int bc, int face = 1);
        void AddStructureLine(const std::vector<Point2d>& points,
                              double hmax = 1e99,
                              int bc = 1,
                              int face = 1);

        void AddSpline(const std::vector<Point2d>& points,
                       double hmax, int bc = 1,
                       int face_left = 1,
                       int face_right = 0);

        void AddCircle(const Point2d& center, double radius,
                       double hmax, int bc = 1,
                       int face_left = 1,
                       int face_right = 0);

        int AddFace(const std::string& name, double maxh_f = 1e99);

        void FakeData();

        char TestComment(std::istream& infile);

        void PartitionBoundary(MeshingParameters& mp, double h, Mesh& mesh2d);

        void GetMaterial(size_t domnr, char*& material);

        double GetDomainMaxh(size_t domnr);

        void SetGrading(const double grading)
        {
            elto0 = grading;
        }

        double GetGrading()
        {
            return elto0;
        }
    };

}  // namespace meshit

#endif
