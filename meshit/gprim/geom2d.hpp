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
#include "geomobjects.hpp"

namespace meshit
{
    class Vec2d;

    class Vec3d;

    class Point2d
    {
        friend class Vec2d;

        friend class Vec3d;

     protected:
        double px, py;

     public:
        Point2d() { /* px = py = 0; */ }

        Point2d(double ax, double ay)
            : px(ax), py(ay) { }

        Point2d(const Point2d& p2)
            : px(p2.px), py(p2.py) { }

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

        double operator[](size_t i) const
        {
            return (i == 0) ? px : (i == 1) ? py : 0.0;
        }

        const Point2d& SetToMin(const Point2d& p2)
        {
            px = std::min(px, p2.px);
            py = std::min(py, p2.py);
            return *this;
        }

        const Point2d& SetToMax(const Point2d& p2)
        {
            px = std::max(px, p2.px);
            py = std::max(py, p2.py);
            return *this;
        }

        friend Vec2d operator-(const Point2d& p1, const Point2d& p2);
        friend Point2d operator-(const Point2d& p1, const Vec2d& v);
        friend Point2d operator+(const Point2d& p1, const Vec2d& v);

        friend double Dist(const Point2d& p1, const Point2d& p2);
        friend double Dist2(const Point2d& p1, const Point2d& p2);

        friend Point2d Center(const Point2d& p1, const Point2d& p2);
        friend Point2d Center(const Point2d& p1, const Point2d& p2, const Point2d& p3);

        friend bool CW(const Point2d& p1, const Point2d& p2, const Point2d& p3);
        friend bool CCW(const Point2d& p1, const Point2d& p2, const Point2d& p3, double eps);

        friend std::ostream& operator<<(std::ostream& s, const Point2d& p);
    };


    inline double Dist(const Point2d& p1, const Point2d& p2)
    {
        return sqrt((p1.px - p2.px) * (p1.px - p2.px) +
                    (p1.py - p2.py) * (p1.py - p2.py));
    }

    inline double Dist2(const Point2d& p1, const Point2d& p2)
    {
        return (p1.px - p2.px) * (p1.px - p2.px) +
               (p1.py - p2.py) * (p1.py - p2.py);
    }

    inline Point2d Center(const Point2d& p1, const Point2d& p2)
    {
        return Point2d(0.5 * (p1.px + p2.px),
                       0.5 * (p1.py + p2.py));
    }

    inline Point2d Center(const Point2d& p1, const Point2d& p2, const Point2d& p3)
    {
        return Point2d(1.0 / 3.0 * (p1.px + p2.px + p3.px),
                       1.0 / 3.0 * (p1.py + p2.py + p3.py));
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
        Vec2d() { vx = vy = 0; }

        Vec2d(const Vec2d& v2)
            : vx{v2.vx}, vy{v2.vy} { }

        Vec2d(double ax, double ay)
            : vx{ax}, vy{ay} { }

        explicit Vec2d(double d)
            : vx{d}, vy{d} { }

        explicit Vec2d(const Point2d& p2)
            : vx{p2.X()}, vy{p2.Y()} { }

        Vec2d(const Point2d& p1, const Point2d& p2)
        {
            vx = p2.px - p1.px;
            vy = p2.py - p1.py;
        }

        Vec2d& operator=(double val)
        {
            vx = vy = val;
            return *this;
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

        friend double Angle(const Vec2d& v1, const Vec2d& v2);
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

    inline double Angle(const Vec2d& v1, const Vec2d& v2)
    {
        double ps = v1.vx * v2.vx + v1.vy * v2.vy;
        double n1 = v1.vx * v1.vx + v1.vy * v1.vy;
        double n2 = v2.vx * v2.vx + v2.vy * v2.vy;
        double norm = sqrt(n1 * n2);
        return acos(std::max(-1.0, std::min(1.0, ps / norm)));
    }

    class Box2d
    {
     protected:
        Point2d pmin, pmax;

     public:
        Box2d() { }

        Box2d(const Point2d& p1, const Point2d& p2)
        {
            pmin.X() = std::min(p1.X(), p2.X());
            pmax.X() = std::max(p1.X(), p2.X());
            pmin.Y() = std::min(p1.Y(), p2.Y());
            pmax.Y() = std::max(p1.Y(), p2.Y());
        }

        const Point2d& PMin() const
        {
            return pmin;
        }

        const Point2d& PMax() const
        {
            return pmax;
        }

        const Box2d& operator+=(const Box2d& b);

        void SetPoint(const Point2d& p)
        {
            pmin = pmax = p;
        }

        void AddPoint(const Point2d& p)
        {
            pmin.X() = std::min(pmin.X(), p.X());
            pmin.Y() = std::min(pmin.Y(), p.Y());
            pmax.X() = std::max(pmax.X(), p.X());
            pmax.Y() = std::max(pmax.Y(), p.Y());
        }

        inline double Diam() const
        {
            double dx = pmax.X() - pmin.X();
            double dy = pmax.Y() - pmin.Y();
            return sqrt(dx * dx + dy * dy);
        }

        inline double LargestSide() const
        {
            return std::max(pmax.X() - pmin.X(),
                            pmax.Y() - pmin.Y());
        }
    };

    inline std::ostream& operator<<(std::ostream& ost, const Box2d& b)
    {
        ost << b.PMin() << " - " << b.PMax();
        return ost;
    }


}  // namespace meshit

#endif
