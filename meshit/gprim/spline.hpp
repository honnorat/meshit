#ifndef FILE_SPLINE_HPP
#define FILE_SPLINE_HPP

/**************************************************************************/
/* File:   spline.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Jul. 96                                                    */
/**************************************************************************/

#define _USE_MATH_DEFINES 1

#include <iostream>
#include <string>
#include <cmath>

#include "../linalg/vector.hpp"
#include "geomobjects.hpp"
#include "geom2d.hpp"

namespace meshit
{
    /*
      Spline curves for 2D mesh generation
     */

    /// Geometry point
    class GeomPoint : public Point2d
    {
     public:
        /// refinement factor at point
        double refatpoint;
        /// max mesh-size at point
        double hmax;

        GeomPoint() { }

        explicit GeomPoint(const Point2d& ap, double aref = 1, double ahmax = 1e99)
            : Point2d{ap}, refatpoint{aref}, hmax{ahmax} { }
    };


    /// base class for 2d - segment

    class SplineSeg
    {
     public:
        SplineSeg() { }

        virtual ~SplineSeg() { }

        virtual double Length() const;

        /// returns point at curve, 0 <= t <= 1
        virtual Point2d GetPoint(double t) const = 0;

        /// returns a (not necessarily unit-length) tangent vector for 0 <= t <= 1
        virtual void GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const = 0;

        /// returns initial point on curve
        virtual const GeomPoint& StartPI() const = 0;
        /// returns terminal point on curve
        virtual const GeomPoint& EndPI() const = 0;

        virtual void GetPoints(size_t n, std::vector<Point2d>& points) const;

        double CalcCurvature(double t) const
        {
            Point2d point;
            Vec2d first, second;
            GetDerivatives(t, point, first, second);
            double fl = first.Length();
            return fabs(first.X() * second.Y() - first.Y() * second.X()) / (fl * fl * fl);
        }

     public:
        size_t leftdom;     // left domain
        size_t rightdom;    // right domain
        double reffak;      // refinement at line
        double hmax;        // maximal h
        int bc;             // boundary condition number
    };

    /// Straight line form p1 to p2
    class LineSeg : public SplineSeg
    {
     public:
        LineSeg(const GeomPoint& ap1, const GeomPoint& ap2)
            : p1(ap1), p2(ap2) { }

        virtual double Length() const;

        inline virtual Point2d GetPoint(double t) const;

        virtual void GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const;

        virtual const GeomPoint& StartPI() const
        {
            return p1;
        };

        virtual const GeomPoint& EndPI() const
        {
            return p2;
        }

     protected:
        GeomPoint p1, p2;
    };

    /// curve given by a rational, quadratic spline (including ellipses)

    class SplineSeg3 : public SplineSeg
    {
     public:
        SplineSeg3(const GeomPoint& ap1,
                   const GeomPoint& ap2,
                   const GeomPoint& ap3);

        inline virtual Point2d GetPoint(double t) const;

        virtual void GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const;

        virtual const GeomPoint& StartPI() const
        {
            return p1;
        };

        virtual const GeomPoint& EndPI() const
        {
            return p3;
        }

     protected:
        GeomPoint p1, p2, p3;
        double weight;
    };

}  // namespace meshit

#endif
