/*
 * Spline curve for Mesh generator
 */

#include "spline.hpp"

namespace meshit {

// calculates length of spline-curve
double SplineSeg::Length() const
{
    uint32_t n = 1000;
    double dt = 1.0 / n;

    Point2d pold = GetPoint(0);

    double l = 0.0;
    for (uint32_t i = 1; i <= n; i++) {
        Point2d p = GetPoint(i * dt);
        l += Dist(p, pold);
        pold = p;
    }

    return l;
}

void SplineSeg::GetPoints(size_t n, std::vector<Point2d>& points) const
{
    points.resize(n);
    if (n >= 2) {
        double nm1 = static_cast<double>(n - 1);
        for (size_t i = 0; i < n; i++) {
            points[i] = GetPoint(static_cast<double>(i) / nm1);
        }
    }
}

/*
   Implementation of line-segment from p1 to p2
 */

inline Point2d LineSeg::GetPoint(double t) const
{
    double px = p1.X() + t * (p2.X() - p1.X());
    double py = p1.Y() + t * (p2.Y() - p1.Y());
    return Point2d(px, py);
}

void LineSeg::GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const
{
    first.X() = p2.X() - p1.X();
    first.Y() = p2.Y() - p1.Y();

    point.X() = p1.X() + t * first.X();
    point.Y() = p1.Y() + t * first.Y();

    second.X() = 0.0;
    second.Y() = 0.0;
}

inline double LineSeg::Length() const
{
    return Dist(p1, p2);
}

SplineSeg3::SplineSeg3(const GeomPoint& ap1, const GeomPoint& ap2, const GeomPoint& ap3)
    : p1(ap1)
    , p2(ap2)
    , p3(ap3)
{
    weight = Dist(p1, p3) / sqrt(0.5 * (Dist2(p1, p2) + Dist2(p2, p3)));
}

inline Point2d SplineSeg3::GetPoint(double t) const
{
    double b1 = (1 - t) * (1 - t);
    double b2 = weight * t * (1 - t);
    double b3 = t * t;
    double w = 1.0 / (b1 + b2 + b3);
    double px = b1 * p1.X() + b2 * p2.X() + b3 * p3.X();
    double py = b1 * p1.Y() + b2 * p2.Y() + b3 * p3.Y();

    return Point2d(w * px, w * py);
}

void SplineSeg3::GetDerivatives(const double t, Point2d& point, Vec2d& first, Vec2d& second) const
{
    Vec2d v1{p1}, v2{p2}, v3{p3};

    double b1 = (1. - t) * (1. - t);
    double b2 = weight * t * (1. - t);
    double b3 = t * t;
    double w = 1.0 / (b1 + b2 + b3);
    b1 *= w;
    b2 *= w;
    b3 *= w;

    double b1p = 2. * (t - 1.);
    double b2p = weight * (1. - 2. * t);
    double b3p = 2. * t;
    const double wp = b1p + b2p + b3p;
    const double fac1 = wp * w;
    b1p *= w;
    b2p *= w;
    b3p *= w;

    const double b1pp = 2.;
    const double b2pp = -2. * weight;
    const double b3pp = 2.;
    const double wpp = b1pp + b2pp + b3pp;
    const double fac2 = wpp * w - 2.0 * wp * wp * w * w;

    point.X() = b1 * p1.X() + b2 * p2.X() + b3 * p3.X();
    point.Y() = b1 * p1.Y() + b2 * p2.Y() + b3 * p3.Y();

    first = (b1p - b1 * fac1) * v1 + (b2p - b2 * fac1) * v2 + (b3p - b3 * fac1) * v3;

    second = (b1pp * w - 2.0 * b1p * fac1 - b1 * fac2) * v1 + (b2pp * w - 2.0 * b2p * fac1 - b2 * fac2) * v2 +
             (b3pp * w - 2.0 * b3p * fac1 - b3 * fac2) * v3;
}

}  // namespace meshit
