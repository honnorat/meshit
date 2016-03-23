#ifndef FILE_GEOM2D_HPP
#define FILE_GEOM2D_HPP

/* *************************************************************************/
/* File:   geom2d.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   5. Aug. 95                                                      */
/* *************************************************************************/

#include <iostream>
#include <algorithm>

#include "../general/index.hpp"
#include "geomfuncs.hpp"

namespace meshit
{
    /* Geometric Algorithms */

    class Point2d;

    class Vec2d;

    Vec2d operator-(const Point2d& p1, const Point2d& p2);
    Point2d operator-(const Point2d& p1, const Vec2d& v);
    Point2d operator+(const Point2d& p1, const Vec2d& v);

    Vec2d operator-(const Point2d& p1, const Point2d& p2);
    Point2d operator-(const Point2d& p1, const Vec2d& v);
    Point2d operator+(const Point2d& p1, const Vec2d& v);
    Vec2d operator-(const Vec2d& p1, const Vec2d& v);
    Vec2d operator+(const Vec2d& p1, const Vec2d& v);
    Vec2d operator*(double scal, const Vec2d& v);

    std::ostream& operator<<(std::ostream& s, const Point2d& p);
    std::ostream& operator<<(std::ostream& s, const Vec2d& v);

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

        explicit Point2d(const Point<2>& p2)
            : px(p2[0]), py(p2[1]) { }

        Point2d& operator=(const Point2d& p2)
        {
            px = p2.px;
            py = p2.py;
            return *this;
        }

        bool operator==(const Point2d& p2) const
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

        friend Vec2d operator-(const Point2d& p1, const Point2d& p2);
        friend Point2d operator-(const Point2d& p1, const Vec2d& v);
        friend Point2d operator+(const Point2d& p1, const Vec2d& v);

        friend double Dist2(const Point2d& p1, const Point2d& p2);

        friend bool CW(const Point2d& p1, const Point2d& p2, const Point2d& p3);
        friend bool CCW(const Point2d& p1, const Point2d& p2, const Point2d& p3, double eps);

        friend std::ostream& operator<<(std::ostream& s, const Point2d& p);
    };


    inline double Dist2(const Point2d& p1, const Point2d& p2)
    {
        return ((p1.px - p2.px) * (p1.px - p2.px) +
                (p1.py - p2.py) * (p1.py - p2.py));
    }

    /** Are the points (p1, p2, p3) clock-wise ?
     */
    inline bool CW(const Point2d& p1, const Point2d& p2, const Point2d& p3)
    {
        // return Cross (p2 - p1, p3 - p2) < 0;
        return (p2.px - p1.px) * (p3.py - p2.py) -
               (p2.py - p1.py) * (p3.px - p2.px) < 0;
    }

    /**  Are the points (p1, p2, p3) counter-clock-wise ?
     */
    inline bool CCW(const Point2d& p1, const Point2d& p2, const Point2d& p3, double eps)
    {
        // return Cross (p2 - p1, p3 - p2) > 0;
        double ax = p2.px - p1.px;
        double ay = p2.py - p1.py;
        double bx = p3.px - p2.px;
        double by = p3.py - p2.py;

        return ax * by - ay * bx > eps * eps * std::max(ax * ax + ay * ay, bx * bx + by * by);
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

        Vec2d& operator+=(const Vec2d& v2);
        Vec2d& operator-=(const Vec2d& v2);
        Vec2d& operator*=(double s);
        Vec2d& operator/=(double s);

        friend Vec2d operator-(const Point2d& p1, const Point2d& p2);
        friend Point2d operator-(const Point2d& p1, const Vec2d& v);
        friend Point2d operator+(const Point2d& p1, const Vec2d& v);
        friend Vec2d operator-(const Vec2d& p1, const Vec2d& v);
        friend Vec2d operator+(const Vec2d& p1, const Vec2d& v);
        friend Vec2d operator*(double scal, const Vec2d& v);
        friend double operator*(const Vec2d& v1, const Vec2d& v2);

        friend double Cross(const Vec2d& v1, const Vec2d& v2);

        friend std::ostream& operator<<(std::ostream& s, const Vec2d& v);
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
        } else {
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

    inline double operator*(const Vec2d& v1, const Vec2d& v2)
    {
        return v1.X() * v2.X() + v1.Y() * v2.Y();
    }

    inline double Cross(const Vec2d& v1, const Vec2d& v2)
    {
        return v1.X() * v2.Y() - v1.Y() * v2.X();
    }

}  // namespace meshit

#endif
