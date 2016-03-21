#ifndef FILE_GEOM2D
#define FILE_GEOM2D

/* *************************************************************************/
/* File:   geom2d.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   5. Aug. 95                                                      */
/* *************************************************************************/

#include <iostream>
#include <algorithm>

#include "../general/template.hpp"
#include "geomfuncs.hpp"

namespace meshit {

    /* Geometric Algorithms */

#define EPSGEOM 1E-5

    class Point2d;

    class Vec2d;

    class LINE2D;

    class Line2d;

    class PLine2d;

    class TRIANGLE2D;

    class PTRIANGLE2D;

    inline Vec2d operator-(const Point2d& p1, const Point2d& p2);
    inline Point2d operator-(const Point2d& p1, const Vec2d& v);
    inline Point2d operator+(const Point2d& p1, const Vec2d& v);
    inline Point2d Center(const Point2d& p1, const Point2d& p2);

    std::ostream& operator<<(std::ostream& s, const Point2d& p);
    inline Vec2d operator-(const Point2d& p1, const Point2d& p2);
    inline Point2d operator-(const Point2d& p1, const Vec2d& v);
    inline Point2d operator+(const Point2d& p1, const Vec2d& v);
    inline Vec2d operator-(const Vec2d& p1, const Vec2d& v);
    inline Vec2d operator+(const Vec2d& p1, const Vec2d& v);
    inline Vec2d operator*(double scal, const Vec2d& v);
    double Angle(const Vec2d& v);
    double Angle(const Vec2d& v1, const Vec2d& v2);
    std::ostream& operator<<(std::ostream& s, const Vec2d& v);
    double Dist2(const Line2d& g, const Line2d& h); // GH
    int Near(const Point2d& p1, const Point2d& p2, const double eps);

    int Parallel(const Line2d& l1, const Line2d& l2, double peps = EPSGEOM);
    int IsOnLine(const Line2d& l, const Point2d& p, double heps = EPSGEOM);
    int IsOnLongLine(const Line2d& l, const Point2d& p);
    int Hit(const Line2d& l1, const Line2d& l2, double heps = EPSGEOM);
    std::ostream& operator<<(std::ostream& s, const Line2d& l);
    Point2d CrossPoint(const PLine2d& l1, const PLine2d& l2);
    Point2d CrossPoint(const Line2d& l1, const Line2d& l2);
    int Parallel(const PLine2d& l1, const PLine2d& l2, double peps = EPSGEOM);
    int IsOnLine(const PLine2d& l, const Point2d& p, double heps = EPSGEOM);
    int IsOnLongLine(const PLine2d& l, const Point2d& p);
    int Hit(const PLine2d& l1, const Line2d& l2, double heps = EPSGEOM);
    std::ostream& operator<<(std::ostream& s, const Line2d& l);
    std::ostream& operator<<(std::ostream& s, const TRIANGLE2D& t);
    std::ostream& operator<<(std::ostream& s, const PTRIANGLE2D& t);
    double Dist2(const Point2d& p1, const Point2d& p2);

    class Point2d
    {
        friend class Vec2d;

     protected:
        double px, py;

     public:
        Point2d() { /* px = py = 0; */ }

        Point2d(double ax, double ay)
                : px(ax), py(ay) { }

        Point2d(const Point2d& p2)
                : px(p2.px), py(p2.py) { }

        Point2d(const Point<2>& p2)
                : px(p2[0]), py(p2[1]) { }

        Point2d& operator=(const Point2d& p2)
        {
            px = p2.px;
            py = p2.py;
            return *this;
        }

        int operator==(const Point2d& p2) const // GH
        {
            return (px == p2.px && py == p2.py);
        }

        double& X()
        {
            return px;
        }

        double& Y()
        {
            return py;
        }

        double X() const
        {
            return px;
        }

        double Y() const
        {
            return py;
        }

        operator Point<2>() const
        {
            return Point<2>(px, py);
        }

        friend inline Vec2d operator-(const Point2d& p1, const Point2d& p2);
        friend inline Point2d operator-(const Point2d& p1, const Vec2d& v);
        friend inline Point2d operator+(const Point2d& p1, const Vec2d& v);

        friend inline Point2d Center(const Point2d& p1, const Point2d& p2);

        const Point2d& SetToMin(const Point2d& p2)
        {
            if (p2.px < px) px = p2.px;
            if (p2.py < py) py = p2.py;
            return *this;
        }

        const Point2d& SetToMax(const Point2d& p2)
        {
            if (p2.px > px) px = p2.px;
            if (p2.py > py) py = p2.py;
            return *this;
        }

        friend double Dist(const Point2d& p1, const Point2d& p2)
        {
            return sqrt((p1.px - p2.px) * (p1.px - p2.px) +
                        (p1.py - p2.py) * (p1.py - p2.py));
        }
        //    { return sqrt ( sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ); }

        friend double Dist2(const Point2d& p1, const Point2d& p2)
        {
            return ((p1.px - p2.px) * (p1.px - p2.px) +
                    (p1.py - p2.py) * (p1.py - p2.py));
        }
        //    { return sqr (p1.X()-p2.X()) + sqr (p1.Y()-p2.Y()) ; }

        /**
           Points clock-wise ?
           Are the points (p1, p2, p3) clock-wise ?
         */
        friend inline int CW(const Point2d& p1, const Point2d& p2, const Point2d& p3)
        {
            //      return Cross (p2 - p1, p3 - p2) < 0;      
            return
                    (p2.px - p1.px) * (p3.py - p2.py) -
                    (p2.py - p1.py) * (p3.px - p2.px) < 0;
        }

        /**
           Points counter-clock-wise ?
           Are the points (p1, p2, p3) counter-clock-wise ?
         */
        friend inline bool CCW(const Point2d& p1, const Point2d& p2, const Point2d& p3)
        {
            //      return Cross (p2 - p1, p3 - p2) > 0;
            return
                    (p2.px - p1.px) * (p3.py - p2.py) -
                    (p2.py - p1.py) * (p3.px - p2.px) > 0;
        }

        /**
             Points counter-clock-wise ?
             Are the points (p1, p2, p3) counter-clock-wise ?
              */

        friend inline bool CCW(const Point2d& p1, const Point2d& p2, const Point2d& p3, double eps)
        {
            //      return Cross (p2 - p1, p3 - p2) > 0;
            double ax = p2.px - p1.px;
            double ay = p2.py - p1.py;
            double bx = p3.px - p2.px;
            double by = p3.py - p2.py;

            return ax * by - ay * bx > eps * eps * std::max(ax * ax + ay * ay, bx * bx + by * by);
        }

        friend std::ostream& operator<<(std::ostream& s, const Point2d& p);
    };

    inline int Near(const Point2d& p1, const Point2d& p2, const double eps = 1e-4)
    {
        return Dist2(p1, p2) <= eps * eps;
    }

    class Vec2d
    {
     protected:
        double vx, vy;

     public:

        Vec2d() { /* vx = vy = 0; */ }

        Vec2d(double ax, double ay)
                : vx{ax}, vy{ay} { }

        Vec2d(const Vec2d& v2)
                : vx{v2.vx}, vy{v2.vy} { }

        explicit Vec2d(const Vec<2>& v2)
                : vx{v2[0]}, vy{v2[1]} { }

        Vec2d(const Point2d& p1, const Point2d& p2)
        {
            vx = p2.px - p1.px;
            vy = p2.py - p1.py;
        }

        Vec2d& operator=(const Vec2d& p2)
        {
            vx = p2.vx;
            vy = p2.vy;
            return *this;
        }

        double& X()
        {
            return vx;
        }

        double& Y()
        {
            return vy;
        }

        double X() const
        {
            return vx;
        }

        double Y() const
        {
            return vy;
        }

        double Length() const
        {
            return sqrt(vx * vx + vy * vy);
        }

        double Length2() const
        {
            return vx * vx + vy * vy;
        }

        void GetNormal(Vec2d& n) const
        {
            n.vx = -vy;
            n.vy = vx;
        } // GH

        inline Vec2d& operator+=(const Vec2d& v2);
        inline Vec2d& operator-=(const Vec2d& v2);
        inline Vec2d& operator*=(double s);
        inline Vec2d& operator/=(double s);

        friend inline Vec2d operator-(const Point2d& p1, const Point2d& p2);
        friend inline Point2d operator-(const Point2d& p1, const Vec2d& v);
        friend inline Point2d operator+(const Point2d& p1, const Vec2d& v);
        friend inline Vec2d operator-(const Vec2d& p1, const Vec2d& v);
        friend inline Vec2d operator+(const Vec2d& p1, const Vec2d& v);
        friend inline Vec2d operator*(double scal, const Vec2d& v);

        friend double operator*(const Vec2d& v1, const Vec2d& v2)
        {
            return v1.X() * v2.X() + v1.Y() * v2.Y();
        }

        friend double Cross(const Vec2d& v1, const Vec2d& v2)
        {
            return double(v1.X()) * double(v2.Y()) -
                   double(v1.Y()) * double(v2.X());
        }

        friend inline void PpSmV(const Point2d& p1, double s, const Vec2d& v, Point2d& p2);
        friend inline void PmP(const Point2d& p1, const Point2d& p2, Vec2d& v);

        ///						Angle in [0,2*PI)

        friend double Angle(const Vec2d& v);
        friend double FastAngle(const Vec2d& v);
        friend double Angle(const Vec2d& v1, const Vec2d& v2);
        friend double FastAngle(const Vec2d& v1, const Vec2d& v2);

        friend std::ostream& operator<<(std::ostream& s, const Vec2d& v);
    };

    class Line2d
    {
     protected:
        Point2d p1, p2;

     public:
        Line2d()
                : p1(), p2() { }

        Line2d(const Point2d& ap1, const Point2d& ap2)
                : p1(ap1), p2(ap2) { }

        Line2d& operator=(const Line2d& l2)
        {
            p1 = l2.p1;
            p2 = l2.p2;
            return *this;
        }

        Point2d& P1()
        {
            return p1;
        }

        Point2d& P2()
        {
            return p2;
        }

        const Point2d& P1() const
        {
            return p1;
        }

        const Point2d& P2() const
        {
            return p2;
        }

        double XMax() const
        {
            return std::max(p1.X(), p2.X());
        }

        double YMax() const
        {
            return std::max(p1.Y(), p2.Y());
        }

        double XMin() const
        {
            return std::min(p1.X(), p2.X());
        }

        double YMin() const
        {
            return std::min(p1.Y(), p2.Y());
        }

        Vec2d Delta() const
        {
            return Vec2d(p2.X() - p1.X(), p2.Y() - p1.Y());
        }

        double Length() const
        {
            return Delta().Length();
        }

        double Length2() const
        {
            return (p1.X() - p2.X()) * (p1.X() - p2.X()) + (p1.Y() - p2.Y()) * (p1.Y() - p2.Y());
        }

        void GetNormal(Line2d& n) const; // GH
        Vec2d NormalDelta() const; // GH

        /// square of the distance between two 2d-lines.
        friend double Dist2(const Line2d& g, const Line2d& h); // GH

        friend Point2d CrossPoint(const Line2d& l1, const Line2d& l2);
        /// returns 1 iff parallel
        friend int CrossPointBarycentric(const Line2d& l1, const Line2d& l2,
                                         double& lam1, double& lam2);

        friend int Parallel(const Line2d& l1, const Line2d& l2, double peps);
        friend int IsOnLine(const Line2d& l, const Point2d& p, double heps);
        friend int IsOnLongLine(const Line2d& l, const Point2d& p);
        friend int Hit(const Line2d& l1, const Line2d& l2, double heps);

        friend std::ostream& operator<<(std::ostream& s, const Line2d& l);
    };

    inline Vec2d& Vec2d::operator+=(const Vec2d& v2)
    {
        vx += v2.vx;
        vy += v2.vy;
        return *this;
    }

    inline Vec2d& Vec2d::operator-=(const Vec2d& v2)
    {
        vx -= v2.vx;
        vy -= v2.vy;
        return *this;
    }

    inline Vec2d& Vec2d::operator*=(double s)
    {
        vx *= s;
        vy *= s;
        return *this;
    }

    inline Vec2d& Vec2d::operator/=(double s)
    {
        if (s != 0) {
            vx /= s;
            vy /= s;
        }
        else {
            std::cerr << "Vec2d::operator /=: Division by zero\n";
        }
        return *this;
    }

    inline Vec2d operator-(const Point2d& p1, const Point2d& p2)
    {
        return Vec2d(p1.X() - p2.X(), p1.Y() - p2.Y());
    }

    inline Point2d operator-(const Point2d& p1, const Vec2d& v)
    {
        return Point2d(p1.X() - v.X(), p1.Y() - v.Y());
    }

    inline Point2d operator+(const Point2d& p1, const Vec2d& v)
    {
        return Point2d(p1.X() + v.X(), p1.Y() + v.Y());
    }

    inline Point2d Center(const Point2d& p1, const Point2d& p2)
    {
        return Point2d((p1.X() + p2.X()) / 2, (p1.Y() + p2.Y()) / 2);
    }

    inline Vec2d operator-(const Vec2d& v1, const Vec2d& v2)
    {
        return Vec2d(v1.X() - v2.X(), v1.Y() - v2.Y());
    }

    inline Vec2d operator+(const Vec2d& v1, const Vec2d& v2)
    {
        return Vec2d(v1.X() + v2.X(), v1.Y() + v2.Y());
    }

    inline Vec2d operator*(double scal, const Vec2d& v)
    {
        return Vec2d(scal * v.X(), scal * v.Y());
    }

}

#endif
