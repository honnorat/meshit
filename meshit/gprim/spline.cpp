/*
 * Spline curve for Mesh generator
 */

#include "spline.hpp"
#include "../linalg/densemat.hpp"

namespace meshit {

    template<>
    void CircleSeg<2>::LineIntersections(const double a, const double b, const double c,
                                         Array<Point<2> >& points, const double eps) const
    {
        points.resize(0);

        double px = 0, py = 0;

        if (fabs(b) > 1e-20)
            py = -c / b;
        else
            px = -c / a;

        const double c1 = a * a + b * b;
        const double c2 = 2. * (a * (py - pm(1)) - b * (px - pm(0)));
        const double c3 = pow(px - pm(0), 2) + pow(py - pm(1), 2) - pow(Radius(), 2);

        const double discr = c2 * c2 - 4 * c1 * c3;

        if (discr < 0)
            return;

        Array<double> t;

        if (fabs(discr) < 1e-20) {
            t.push_back(-0.5 * c2 / c1);
        } else {
            t.push_back((-c2 + sqrt(discr)) / (2. * c1));
            t.push_back((-c2 - sqrt(discr)) / (2. * c1));
        }

        for (int i = 0; i < t.size(); i++) {
            Point<2> p(px - t[i] * b, py + t[i] * a);

            double angle = atan2(p(1), p(0)) + M_PI;

            if (angle > StartAngle() - eps && angle < EndAngle() + eps)
                points.push_back(p);
        }
    }

    template<int D>
    SplineSeg3<D>::SplineSeg3(const GeomPoint<D>& ap1,
                              const GeomPoint<D>& ap2,
                              const GeomPoint<D>& ap3)
            : p1(ap1), p2(ap2), p3(ap3)
    {
        weight = Dist(p1, p3) / sqrt(0.5 * (Dist2(p1, p2) + Dist2(p2, p3)));
        proj_latest_t = 0.5;
    }

    template<int D>
    inline Point<D> SplineSeg3<D>::GetPoint(double t) const
    {
        double b1, b2, b3;

        b1 = (1 - t) * (1 - t);
        b2 = weight * t * (1 - t);
        b3 = t * t;

        Vec<D> hp = b1 * Vec<D>(p1) + b2 * Vec<D>(p2) + b3 * Vec<D>(p3);
        double w = b1 + b2 + b3;
        return Point<D>((1.0 / w) * hp);
    }

    template<int D>
    Vec<D> SplineSeg3<D>::GetTangent(const double t) const
    {
        const double b1 = (1. - t) * ((weight - 2.) * t - weight);
        const double b2 = weight * (1. - 2. * t);
        const double b3 = t * ((weight - 2) * t + 2.);


        Vec<D> retval;
        for (int i = 0; i < D; i++)
            retval(i) = b1 * p1(i) + b2 * p2(i) + b3 * p3(i);

        return retval;

    }

    template<int D>
    void SplineSeg3<D>::Project(const Point<D>& point, Point<D>& point_on_curve, double& t) const
    {
        double t_old = -1;

        if (proj_latest_t > 0. && proj_latest_t < 1.)
            t = proj_latest_t;
        else
            t = 0.5;

        Point<D> phi;
        Vec<D> phip, phipp, phimp;

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

            point_on_curve = SplineSeg3<D>::GetPoint(t);

            double dist = Dist(point, point_on_curve);

            phi = SplineSeg3<D>::GetPoint(0);
            double auxdist = Dist(phi, point);
            if (auxdist < dist) {
                t = 0.;
                point_on_curve = phi;
                dist = auxdist;
            }
            phi = SplineSeg3<D>::GetPoint(1);
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
                phi = SplineSeg3<D>::GetPoint(t0);
                d0 = Dist(phi, point);
                phi = SplineSeg3<D>::GetPoint(t1);
                d1 = Dist(phi, point);
                phi = SplineSeg3<D>::GetPoint(t2);
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

            phi = SplineSeg3<D>::GetPoint(t0);
            d0 = Dist(phi, point);
            phi = SplineSeg3<D>::GetPoint(t1);
            d1 = Dist(phi, point);
            phi = SplineSeg3<D>::GetPoint(t2);
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

            point_on_curve = SplineSeg3<D>::GetPoint(t);
        }
        proj_latest_t = t;
    }

    template<int D>
    void SplineSeg3<D>::GetDerivatives(const double t,
                                       Point<D>& point,
                                       Vec<D>& first,
                                       Vec<D>& second) const
    {
        Vec<D> v1(p1), v2(p2), v3(p3);

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

        for (int i = 0; i < D; i++)
            point(i) = b1 * p1(i) + b2 * p2(i) + b3 * p3(i);


        first = (b1p - b1 * fac1) * v1 +
                (b2p - b2 * fac1) * v2 +
                (b3p - b3 * fac1) * v3;

        second = (b1pp / w - 2 * b1p * fac1 - b1 * fac2) * v1 +
                 (b2pp / w - 2 * b2p * fac1 - b2 * fac2) * v2 +
                 (b3pp / w - 2 * b3p * fac1 - b3 * fac2) * v3;
    }

    template<int D>
    void SplineSeg3<D>::LineIntersections(const double a, const double b, const double c,
                                          Array<Point<D> >& points, const double eps) const
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

    template
    class SplineSeg3<2>;

    template
    class SplineSeg3<3>;

}  // namespace meshit
