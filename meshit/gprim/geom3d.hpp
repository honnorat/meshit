#ifndef FILE_GEOM3D_HPP
#define FILE_GEOM3D_HPP

/* *************************************************************************/
/* File:   geom3d.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   5. Aug. 95                                                      */
/* *************************************************************************/

#include <iostream>
#include <fstream>

#include "geom2d.hpp"

namespace meshit
{
    class Vec3d
    {
     protected:
        double x[3];

     public:
        Vec3d()
        {
            x[0] = x[1] = x[2] = 0.0;
        }

        Vec3d(double ax, double ay)
        {
            x[0] = ax;
            x[1] = ay;
            x[2] = 0.0;
        }

        Vec3d(double ax, double ay, double az)
        {
            x[0] = ax;
            x[1] = ay;
            x[2] = az;
        }

        explicit Vec3d(double ax[3])
        {
            x[0] = ax[0];
            x[1] = ax[1];
            x[2] = ax[2];
        }

        Vec3d(const Point2d& p1, const Point2d& p2)
        {
            x[0] = p2.px - p1.px;
            x[1] = p2.py - p1.py;
            x[2] = 0.0;
        }

        Vec3d(const Vec3d& v2)
        {
            x[0] = v2.x[0];
            x[1] = v2.x[1];
            x[2] = v2.x[2];
        }

        Vec3d(const Vec2d& v2)
        {
            x[0] = v2.X();
            x[1] = v2.Y();
            x[2] = 0.0;
        }

        Vec3d& operator=(const Vec3d& v2)
        {
            x[0] = v2.x[0];
            x[1] = v2.x[1];
            x[2] = v2.x[2];
            return *this;
        }

        Vec3d& operator=(double val)
        {
            x[0] = x[1] = x[2] = val;
            return *this;
        }

        double& X()
        {
            return x[0];
        }

        double& Y()
        {
            return x[1];
        }

        double& Z()
        {
            return x[2];
        }

        double X() const
        {
            return x[0];
        }

        double Y() const
        {
            return x[1];
        }

        double Z() const
        {
            return x[2];
        }

        double Length() const
        {
            return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        }

        Vec3d& operator+=(const Vec3d& v2);
        Vec3d& operator-=(const Vec3d& v2);
        Vec3d& operator*=(double s);
        Vec3d& operator/=(double s);

        friend Vec3d operator-(const Vec3d& p1, const Vec3d& v);
        friend Vec3d operator+(const Vec3d& p1, const Vec3d& v);
        friend Vec3d operator*(double scal, const Vec3d& v);
        friend Vec3d operator*(const Mat3x3& m, const Vec3d& v);
        friend double operator*(const Vec3d& v1, const Vec3d& v2);

        friend std::ostream& operator<<(std::ostream& s, const Vec3d& v);
    };

//    inline Vec3d& Vec3d::operator+=(const Vec3d& v2)
//    {
//        x[0] += v2.X();
//        x[1] += v2.Y();
//        x[2] += v2.Z();
//        return *this;
//    }
//
//    inline Vec3d& Vec3d::operator-=(const Vec3d& v2)
//    {
//        x[0] -= v2.X();
//        x[1] -= v2.Y();
//        x[2] -= v2.Z();
//        return *this;
//    }
//
//    inline Vec3d& Vec3d::operator*=(double s)
//    {
//        x[0] *= s;
//        x[1] *= s;
//        x[2] *= s;
//        return *this;
//    }
//
//    inline Vec3d& Vec3d::operator/=(double s)
//    {
//        if (s != 0) {
//            x[0] /= s;
//            x[1] /= s;
//            x[2] /= s;
//        }
//        return *this;
//    }
//
//    inline Vec3d operator-(const Vec3d& v1, const Vec3d& v2)
//    {
//        return Vec3d(v1.x[0] - v2.x[0],
//                     v1.x[1] - v2.x[1],
//                     v1.x[2] - v2.x[2]);
//    }
//
//    inline Vec3d operator+(const Vec3d& v1, const Vec3d& v2)
//    {
//        return Vec3d(v1.x[0] + v2.x[0],
//                     v1.x[1] + v2.x[1],
//                     v1.x[2] + v2.x[2]);
//    }
//
//    inline Vec3d operator*(double scal, const Vec3d& v)
//    {
//        return Vec3d(scal * v.x[0],
//                     scal * v.x[1],
//                     scal * v.x[2]);
//    }
//
//    inline double operator*(const Vec3d& v1, const Vec3d& v2)
//    {
//        return v1.x[0] * v2.x[0] +
//               v1.x[1] * v2.x[1] +
//               v1.x[2] * v2.x[2];
//    }
//
//
//    inline Vec3d operator*(const Mat3x3& m, const Vec3d& v)
//    {
//        Vec3d res;
//        res.X() = m(0, 0) * v.X() + m(0, 1) * v.Y() + m(0, 2) * v.Z();
//        res.Y() = m(1, 0) * v.X() + m(1, 1) * v.Y() + m(1, 2) * v.Z();
//        res.Z() = m(2, 0) * v.X() + m(2, 1) * v.Y() + m(2, 2) * v.Z();
//        return res;
//    }

}  // namespace meshit


#endif
