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
#include "../gprim/splinegeometry.hpp"
#include "../gprim/geomobjects.hpp"
#include "../meshing/refine.hpp"

namespace meshit
{
    void Optimize2d(Mesh& mesh, MeshingParameters& mp);

    class SplineSegExt : public SplineSeg
    {
     public:
        const SplineSeg& seg;

        size_t leftdom;     // left domain
        size_t rightdom;    // right domain
        double reffak;      // refinement at line
        double hmax;        // maximal h
        int bc;             // boundary condition number

        explicit SplineSegExt(const SplineSeg& hseg)
            : seg(hseg) { }

        ~SplineSegExt()
        {
            delete &seg;
        }

        virtual const GeomPoint<2>& StartPI() const
        {
            return seg.StartPI();
        }

        virtual const GeomPoint<2>& EndPI() const
        {
            return seg.EndPI();
        }

        virtual Point<2> GetPoint(double t) const
        {
            return seg.GetPoint(t);
        }

        virtual void GetDerivatives(const double t,
                                    Point<2>& point,
                                    Vec<2>& first,
                                    Vec<2>& second) const
        {
            seg.GetDerivatives(t, point, first, second);
        }

        virtual void GetPoints(size_t n, std::vector<Point<2> >& points) const
        {
            seg.GetPoints(n, points);
        }

        virtual double CalcCurvature(double t) const
        {
            Point<2> point;
            Vec<2> first, second;
            GetDerivatives(t, point, first, second);
            double fl = first.Length();
            return fabs(first[0] * second[1] - first[1] * second[0]) / (fl * fl * fl);
        }
    };

    class SplineGeometry2d : public SplineGeometry
    {
     protected:
        std::vector<char*> materials;
        std::vector<double> maxh;
        double elto0;

     public:
        SplineGeometry2d() : elto0{0.3} { }

        virtual ~SplineGeometry2d();

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

        int AddFace(const char* name, double maxh_f = 1e99);

        void FakeData();

        void TestComment(std::istream& infile);

        SplineSegExt& GetSpline(const size_t i)
        {
            return dynamic_cast<SplineSegExt&> (*splines[i]);
        }

        int GenerateMesh(Mesh*& mesh, MeshingParameters& mp);

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
