#ifndef MESHIT_SPLINE_HPP
#define MESHIT_SPLINE_HPP
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

#include <cmath>
#include <iostream>
#include <string>

#include "../linalg/vector.hpp"
#include "../meshing/mesh_class.hpp"
#include "geom2d.hpp"
#include "geomobjects.hpp"

namespace meshit {

/*
  Spline curves for 2D mesh generation
 */

// Geometry point
class GeomPoint : public Point2d
{
 public:
    GeomPoint() { }

    explicit GeomPoint(const Point2d& ap, double aref = 1, double ahmax = 1e99)
        : Point2d{ap}, refatpoint{aref}, hmax{ahmax} { }

 public:
    double refatpoint;  // refinement factor at point
    double hmax;        // max mesh-size at point
};

// base class for 2d-segment
class SplineSeg
{
 public:
    SplineSeg() { }

    virtual ~SplineSeg() { }

    virtual double Length() const;

    // returns point at curve, 0 <= t <= 1
    virtual Point2d GetPoint(double t) const = 0;

    // returns a (not necessarily unit-length) tangent vector for 0 <= t <= 1
    virtual void GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const = 0;

    virtual const GeomPoint& StartPI() const = 0;  // returns initial point on curve
    virtual const GeomPoint& EndPI() const = 0;    // returns terminal point on curve

    void GetPoints(size_t n, std::vector<Point2d>& points) const;

    double CalcCurvature(double t) const
    {
        Point2d point;
        Vec2d first, second;
        GetDerivatives(t, point, first, second);
        double fl = first.Length();
        return fabs(first.X() * second.Y() - first.Y() * second.X()) / (fl * fl * fl);
    }

    EdgeIndex GetID() const { return id_; }
    void SetID(EdgeIndex id) { id_ = id; }

    void SetDomains(DomainIndex left, DomainIndex right)
    {
        dom_left = left;
        dom_right = right;
    }

    void SetHRef(double max_h, double ref_fac = 1.0)
    {
        hmax_ = max_h;
        ref_fac_ = ref_fac;
    }

 public:
    DomainIndex dom_left;   // left domain
    DomainIndex dom_right;  // right domain
    double ref_fac_;        // refinement at line
    double hmax_;           // maximal h
    EdgeIndex id_;          // spline index number
};

// Straight line form p1 to p2
class LineSeg : public SplineSeg
{
 public:
    LineSeg(const GeomPoint& ap1, const GeomPoint& ap2)
        : p1(ap1), p2(ap2) { }

    virtual double Length() const;

    Point2d GetPoint(double t) const override;

    void GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const override;

    const GeomPoint& StartPI() const override { return p1; }
    const GeomPoint& EndPI() const override { return p2; }

 protected:
    GeomPoint p1, p2;
};

/// curve given by a rational, quadratic spline (including ellipses)

class SplineSeg3 : public SplineSeg
{
 public:
    SplineSeg3(const GeomPoint& ap1, const GeomPoint& ap2, const GeomPoint& ap3);

    Point2d GetPoint(double t) const override;

    void GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const override;

    const GeomPoint& StartPI() const override { return p1; }
    const GeomPoint& EndPI() const override { return p3; }

 protected:
    GeomPoint p1, p2, p3;
    double weight;
};

}  // namespace meshit

#endif
