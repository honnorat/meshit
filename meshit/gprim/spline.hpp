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
#include "geomfuncs.hpp"
#include "geom2d.hpp"

namespace meshit {

    /*
      Spline curves for 2D mesh generation
     */

    /// Geometry point

    template<int D>
    class GeomPoint : public Point<D>
    {
     public:
        /// refinement factor at point
        double refatpoint;
        /// max mesh-size at point
        double hmax;

        GeomPoint() { }

        GeomPoint(const Point<D>& ap, double aref = 1, double ahmax=1e99)
                : Point<D>(ap), refatpoint{aref}, hmax{ahmax} { }
    };


    /// base class for 2d - segment

    class SplineSeg
    {
     public:
        SplineSeg() { }

        virtual ~SplineSeg() { }

        virtual double Length() const;

        /// returns point at curve, 0 <= t <= 1
        virtual Point<2> GetPoint(double t) const = 0;

        /// returns a (not necessarily unit-length) tangent vector for 0 <= t <= 1

        virtual Vec<2> GetTangent(const double t) const
        {
            std::cerr << "GetTangent not implemented for spline base-class" << std::endl;
            Vec<2> dummy;
            return dummy;
        }

        virtual void GetDerivatives(const double t,
                                    Point<2>& point,
                                    Vec<2>& first,
                                    Vec<2>& second) const
        {
            double eps = 1e-6;
            point = GetPoint(t);
            Point<2> pl = GetPoint(t - eps);
            Point<2> pr = GetPoint(t + eps);
            first = 1.0 / (2 * eps) * (pr - pl);
            second = 1.0 / (eps * eps) * ((pr - point) + (pl - point));
        }

        /// returns initial point on curve
        virtual const GeomPoint<2>& StartPI() const = 0;
        /// returns terminal point on curve
        virtual const GeomPoint<2>& EndPI() const = 0;

        virtual void GetPoints(size_t n, std::vector<Point<2> >& points) const;
    };

    /// Straight line form p1 to p2
    class LineSeg : public SplineSeg
    {
     public:
        LineSeg(const GeomPoint<2>& ap1, const GeomPoint<2>& ap2)
                : p1(ap1), p2(ap2) { }

        virtual double Length() const;

        inline virtual Point<2> GetPoint(double t) const;

        virtual Vec<2> GetTangent(const double t) const;

        virtual void GetDerivatives(const double t,
                                    Point<2>& point,
                                    Vec<2>& first,
                                    Vec<2>& second) const;

        virtual const GeomPoint<2>& StartPI() const
        {
            return p1;
        };

        virtual const GeomPoint<2>& EndPI() const
        {
            return p2;
        }

        virtual std::string GetType(void) const
        {
            return "line";
        }

     protected:
        GeomPoint<2> p1, p2;
    };

    /// curve given by a rational, quadratic spline (including ellipses)

    class SplineSeg3 : public SplineSeg
    {
        GeomPoint<2> p1, p2, p3;
        double weight;
        mutable double proj_latest_t;
     public:

        SplineSeg3(const GeomPoint<2>& ap1,
                   const GeomPoint<2>& ap2,
                   const GeomPoint<2>& ap3);

        inline virtual Point<2> GetPoint(double t) const;

        virtual Vec<2> GetTangent(const double t) const;

        virtual void GetDerivatives(const double t,
                                    Point<2>& point,
                                    Vec<2>& first,
                                    Vec<2>& second) const;

        virtual const GeomPoint<2>& StartPI() const
        {
            return p1;
        };

        virtual const GeomPoint<2>& EndPI() const
        {
            return p3;
        }

        virtual std::string GetType(void) const
        {
            return "spline3";
        }

        const GeomPoint<2>& TangentPoint(void) const
        {
            return p2;
        }
    };
}  // namespace meshit

#endif
