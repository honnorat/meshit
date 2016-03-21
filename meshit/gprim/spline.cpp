/*
 * Spline curve for Mesh generator
 */

#include "spline.hpp"
#include "../linalg/densemat.hpp"

namespace meshit {

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

    void SplineSeg::GetPoints(int n, Array<Point<2> >& points) const
    {
        points.resize(n);
        if (n >= 2) {
            for (int i = 0; i < n; i++) {
                points[i] = GetPoint(double(i) / (n - 1));
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

    Vec<2> LineSeg::GetTangent(const double t) const
    {
        return p2 - p1;
    }

    void LineSeg::GetDerivatives(
            const double t,
            Point<2>& point,
            Vec<2>& first,
            Vec<2>& second) const
    {
        first = p2 - p1;
        point = p1 + t * first;
        second = 0;
    }

    double LineSeg::Length() const
    {
        return Dist(p1, p2);
    }

    void LineSeg::LineIntersections(const double a, const double b, const double c,
                                    Array<Point<2> >& points, const double eps) const
    {
        points.resize(0);

        double denom = -a * p2(0) + a * p1(0) - b * p2(1) + b * p1(1);
        if (fabs(denom) < 1e-20)
            return;

        double t = (a * p1(0) + b * p1(1) + c) / denom;
        if ((t > -eps) && (t < 1. + eps))
            points.push_back(GetPoint(t));
    }

    void LineSeg::Project(const Point<2>& point, Point<2>& point_on_curve, double& t) const
    {
        Vec<2> v = p2 - p1;
        double l = v.Length();
        v *= 1. / l;
        t = (point - p1) * v;

        if (t < 0) t = 0;
        if (t > l) t = l;

        point_on_curve = p1 + t * v;

        t *= 1. / l;
    }

    SplineSeg3::SplineSeg3(const GeomPoint<2>& ap1,
                           const GeomPoint<2>& ap2,
                           const GeomPoint<2>& ap3)
            : p1(ap1), p2(ap2), p3(ap3)
    {
        weight = Dist(p1, p3) / sqrt(0.5 * (Dist2(p1, p2) + Dist2(p2, p3)));
        proj_latest_t = 0.5;
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

    Vec<2> SplineSeg3::GetTangent(const double t) const
    {
        const double b1 = (1. - t) * ((weight - 2.) * t - weight);
        const double b2 = weight * (1. - 2. * t);
        const double b3 = t * ((weight - 2) * t + 2.);


        Vec<2> retval;
        for (int i = 0; i < 2; i++)
            retval(i) = b1 * p1(i) + b2 * p2(i) + b3 * p3(i);

        return retval;

    }

    void SplineSeg3::Project(const Point<2>& point, Point<2>& point_on_curve, double& t) const
    {
        double t_old = -1;

        if (proj_latest_t > 0. && proj_latest_t < 1.)
            t = proj_latest_t;
        else
            t = 0.5;

        Point<2> phi;
        Vec<2> phip, phipp, phimp;

        int i = 0;

        while (t > -0.5 && t < 1.5 && i < 20 && fabs(t - t_old) > 1e-15) {
            GetDerivatives(t, phi, phip, phipp);

            t_old = t;

            phimp = phi - point;

            // t = min2(std::max(t-(phip*phimp)/(phipp*phimp + phip*phip),0.),1.);
            t -= (phip * phimp) / (phipp * phimp + phip * phip);

            i++;
        }

        if (i < 20 && t > -0.4 && t < 1.4) {
            t = std::max(0.0, std::min(t, 1.0));

            point_on_curve = SplineSeg3::GetPoint(t);

            double dist = Dist(point, point_on_curve);

            phi = SplineSeg3::GetPoint(0);
            double auxdist = Dist(phi, point);
            if (auxdist < dist) {
                t = 0.;
                point_on_curve = phi;
                dist = auxdist;
            }
            phi = SplineSeg3::GetPoint(1);
            auxdist = Dist(phi, point);
            if (auxdist < dist) {
                t = 1.;
                point_on_curve = phi;
                dist = auxdist;
            }
        } else {
            double t0 = 0.0;
            double t1 = 0.5;
            double t2 = 1.0;
            double d0, d1, d2;

            while (t2 - t0 > 1e-8) {
                phi = SplineSeg3::GetPoint(t0);
                d0 = Dist(phi, point);
                phi = SplineSeg3::GetPoint(t1);
                d1 = Dist(phi, point);
                phi = SplineSeg3::GetPoint(t2);
                d2 = Dist(phi, point);

                double a = (2. * d0 - 4. * d1 + 2. * d2) / pow(t2 - t0, 2);

                if (a <= 0) {
                    if (d0 < d2)
                        t2 -= 0.3 * (t2 - t0);
                    else
                        t0 += 0.3 * (t2 - t0);

                    t1 = 0.5 * (t2 + t0);
                }
                else {
                    double b = (d1 - d0 - a * (t1 * t1 - t0 * t0)) / (t1 - t0);
                    double auxt1 = -0.5 * b / a;

                    if (auxt1 < t0) {
                        t2 -= 0.4 * (t2 - t0);
                        t0 = std::max(0., t0 - 0.1 * (t2 - t0));
                    }
                    else if (auxt1 > t2) {
                        t0 += 0.4 * (t2 - t0);
                        t2 = std::min(1., t2 + 0.1 * (t2 - t0));
                    }
                    else {
                        t1 = auxt1;
                        auxt1 = 0.25 * (t2 - t0);
                        t0 = std::max(0., t1 - auxt1);
                        t2 = std::min(1., t1 + auxt1);
                    }

                    t1 = 0.5 * (t2 + t0);
                }
            }

            phi = SplineSeg3::GetPoint(t0);
            d0 = Dist(phi, point);
            phi = SplineSeg3::GetPoint(t1);
            d1 = Dist(phi, point);
            phi = SplineSeg3::GetPoint(t2);
            d2 = Dist(phi, point);

            double mind = d0;
            t = t0;
            if (d1 < mind) {
                t = t1;
                mind = d1;
            }
            if (d2 < mind) {
                t = t2;
                mind = d2;
            }

            point_on_curve = SplineSeg3::GetPoint(t);
        }
        proj_latest_t = t;
    }

    void SplineSeg3::GetDerivatives(const double t,
                                    Point<2>& point,
                                    Vec<2>& first,
                                    Vec<2>& second) const
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

        for (int i = 0; i < 2; i++)
            point(i) = b1 * p1(i) + b2 * p2(i) + b3 * p3(i);


        first = (b1p - b1 * fac1) * v1 +
                (b2p - b2 * fac1) * v2 +
                (b3p - b3 * fac1) * v3;

        second = (b1pp / w - 2 * b1p * fac1 - b1 * fac2) * v1 +
                 (b2pp / w - 2 * b2p * fac1 - b2 * fac2) * v2 +
                 (b3pp / w - 2 * b3p * fac1 - b3 * fac2) * v3;
    }

    void SplineSeg3::LineIntersections(const double a, const double b, const double c,
                                       Array<Point<2> >& points, const double eps) const
    {
        points.resize(0);

        double t;

        const double c1 = a * p1(0) - weight * a * p2(0) + a * p3(0)
                          + b * p1(1) - weight * b * p2(1) + b * p3(1)
                          + (2. - weight) * c;
        const double c2 =
                -2. * a * p1(0) + weight * a * p2(0) - 2. * b * p1(1) + weight * b * p2(1) + (weight - 2.) * c;
        const double c3 = a * p1(0) + b * p1(1) + c;

        if (fabs(c1) < 1e-20) {
            if (fabs(c2) < 1e-20)
                return;

            t = -c3 / c2;
            if ((t > -eps) && (t < 1. + eps))
                points.push_back(GetPoint(t));
            return;
        }

        const double discr = c2 * c2 - 4. * c1 * c3;

        if (discr < 0)
            return;

        if (fabs(discr / (c1 * c1)) < 1e-14) {
            t = -0.5 * c2 / c1;
            if ((t > -eps) && (t < 1. + eps))
                points.push_back(GetPoint(t));
            return;
        }

        t = (-c2 + sqrt(discr)) / (2. * c1);
        if ((t > -eps) && (t < 1. + eps))
            points.push_back(GetPoint(t));

        t = (-c2 - sqrt(discr)) / (2. * c1);
        if ((t > -eps) && (t < 1. + eps))
            points.push_back(GetPoint(t));
    }
}  // namespace meshit
