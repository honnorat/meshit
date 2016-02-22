#ifndef MESHIT_GEOM2D_GEOMETRY2D_H
#define MESHIT_GEOM2D_GEOMETRY2D_H

/* *************************************************************************/
/* File:   geometry2d.hpp                                                  */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

#include <vector>
#include <string>

#include "../gprim/spline.hpp"
#include "../gprim/splinegeometry.hpp"
#include "geom2dmesh.hpp"

namespace meshit {

    void Optimize2d(Mesh& mesh, MeshingParameters& mp);

    class SplineSegExt : public SplineSeg<2>
    {
     public:
        const SplineSeg<2>& seg;

        int leftdom;    // left domain
        int rightdom;   // right domain
        double reffak;  // refinement at line
        double hmax;    // maximal h
        int bc;         // boundary condition number
        int copyfrom;   // copy spline mesh from other spline (-1.. do not copy)
        bool hpref_left;    // perfrom anisotropic refinement (hp-refinement) to edge
        bool hpref_right;   // perfrom anisotropic refinement (hp-refinement) to edge
        int layer;

        explicit SplineSegExt(const SplineSeg<2>& hseg)
                : seg(hseg)
        {
            layer = 1;
            copyfrom = -1;
        }

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

        virtual Vec<2> GetTangent(const double t) const
        {
            return seg.GetTangent(t);
        }

        virtual void GetDerivatives(const double t,
                                    Point<2>& point,
                                    Vec<2>& first,
                                    Vec<2>& second) const
        {
            seg.GetDerivatives(t, point, first, second);
        }

        virtual void GetPoints(int n, Array<Point<2> >& points) const
        {
            seg.GetPoints(n, points);
        }

        virtual double CalcCurvature(double t) const
        {
            Point<2> point;
            Vec<2> first, second;
            GetDerivatives(t, point, first, second);
            double fl = first.Length();
            return fabs(first(0) * second(1) - first(1) * second(0)) / (fl * fl * fl);
        }

    };

    class SplineGeometry2d : public SplineGeometry<2>
    {
     protected:
        std::vector<char*> materials;
        std::vector<double> maxh;
        std::vector<bool> quadmeshing;
        std::vector<bool> tensormeshing;
        std::vector<int> layer;
        std::vector<std::string*> bcnames;
        double elto0;

     public:
        SplineGeometry2d() : elto0(1.0) { }

        virtual ~SplineGeometry2d();

        void Load(const std::string& filename);
        void LoadData(std::istream& infile);

        void AddLine(const std::vector<Point2d>& point_list,
                     double hmax = 1e99,
                     bool hole = false,
                     int bc = 1);
        void AddStructureLine(const std::vector<Point2d>& point_list,
                              double hmax = 1e99,
                              int bc = 1);
        void FakeData();

        void TestComment(std::istream& infile);

        SplineSegExt& GetSpline(const size_t i)
        {
            return dynamic_cast<SplineSegExt&> (*splines[i]);
        }

        int GenerateMesh(Mesh*& mesh, MeshingParameters& mp);

        void PartitionBoundary(MeshingParameters& mp, double h, Mesh& mesh2d);

        void CopyEdgeMesh(int from, int to, Mesh& mesh2d, Point3dTree& searchtree);

        void GetMaterial(const int domnr, char*& material);

        double GetDomainMaxh(const int domnr);

        bool GetDomainQuadMeshing(int domnr)
        {
            if (quadmeshing.size()) return quadmeshing[domnr - 1];
            else return false;
        }

        bool GetDomainTensorMeshing(int domnr)
        {
            if (tensormeshing.size()) return tensormeshing[domnr - 1];
            else return false;
        }

        int GetDomainLayer(int domnr)
        {
            if (layer.size()) return layer[domnr - 1];
            else return 1;
        }

        void SetGrading(const double grading)
        {
            elto0 = grading;
        }

        std::string GetBCName(const int bcnr) const;

        std::string* BCNamePtr(const int bcnr);

        Refinement& GetRefinement() const;
    };

}  // namespace meshit

#endif
