/*
 * Spline curve for Mesh generator
 */

#include "spline.hpp"

namespace meshit
{
    // calculates length of spline-curve

    double SplineSeg::Length() const
    {
        int n = 100;
        double dt = 1.0 / n;

        Point<2> pold = GetPoint(0);

        double l = 0;
        for (int i = 1; i <= n; i++) {
            Point<2> p = GetPoint(i * dt);
            l += Dist(p, pold);
            pold = p;
        }

        return l;
    }

    void SplineSeg::GetPoints(size_t n, std::vector<Point<2> >& points) const
    {
        points.resize(n);
        if (n >= 2) {
            for (size_t i = 0; i < n; i++) {
                points[i] = GetPoint(static_cast<double>(i) / static_cast<double>(n - 1));
            }
        }
    }

    /*
       Implementation of line-segment from p1 to p2
     */

    inline Point<2> LineSeg::GetPoint(double t) const
    {
        return p1 + t * (p2 - p1);
    }

    void LineSeg::GetDerivatives(const double t, Point<2>& point, Vec<2>& first, Vec<2>& second) const
    {
        first = p2 - p1;
        point = p1 + t * first;
        second = 0;
    }

    double LineSeg::Length() const
    {
        return Dist(p1, p2);
    }

    SplineSeg3::SplineSeg3(const GeomPoint<2>& ap1,
                           const GeomPoint<2>& ap2,
                           const GeomPoint<2>& ap3)
        : p1(ap1), p2(ap2), p3(ap3)
    {
        weight = Dist(p1, p3) / sqrt(0.5 * (Dist2(p1, p2) + Dist2(p2, p3)));
    }

    inline Point<2> SplineSeg3::GetPoint(double t) const
    {
        double b1, b2, b3;

        b1 = (1 - t) * (1 - t);
        b2 = weight * t * (1 - t);
        b3 = t * t;

        Vec<2> hp = b1 * Vec<2>(p1) + b2 * Vec<2>(p2) + b3 * Vec<2>(p3);
        double w = b1 + b2 + b3;
        return Point<2>((1.0 / w) * hp);
    }

    void SplineSeg3::GetDerivatives(const double t, Point<2>& point, Vec<2>& first, Vec<2>& second) const
    {
        Vec<2> v1(p1), v2(p2), v3(p3);

        double b1 = (1. - t) * (1. - t);
        double b2 = weight * t * (1. - t);
        double b3 = t * t;
        double w = b1 + b2 + b3;
        b1 *= 1. / w;
        b2 *= 1. / w;
        b3 *= 1. / w;

        double b1p = 2. * (t - 1.);
        double b2p = weight * (1. - 2. * t);
        double b3p = 2. * t;
        const double wp = b1p + b2p + b3p;
        const double fac1 = wp / w;
        b1p *= 1. / w;
        b2p *= 1. / w;
        b3p *= 1. / w;

        const double b1pp = 2.;
        const double b2pp = -2. * weight;
        const double b3pp = 2.;
        const double wpp = b1pp + b2pp + b3pp;
        const double fac2 = (wpp * w - 2. * wp * wp) / (w * w);

        point[0] = b1 * p1[0] + b2 * p2[0] + b3 * p3[0];
        point[1] = b1 * p1[1] + b2 * p2[1] + b3 * p3[1];

        first = (b1p - b1 * fac1) * v1 +
                (b2p - b2 * fac1) * v2 +
                (b3p - b3 * fac1) * v3;

        second = (b1pp / w - 2 * b1p * fac1 - b1 * fac2) * v1 +
                 (b2pp / w - 2 * b2p * fac1 - b2 * fac2) * v2 +
                 (b3pp / w - 2 * b3p * fac1 - b3 * fac2) * v3;
    }

}  // namespace meshit
