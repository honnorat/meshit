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
    class Point3d;

    class Vec3d;

    inline Vec3d operator-(const Point3d& p1, const Point3d& p2);
    inline Point3d operator-(const Point3d& p1, const Vec3d& v);
    inline Point3d operator+(const Point3d& p1, const Vec3d& v);
    inline Point3d Center(const Point3d& p1, const Point3d& p2);
    inline Point3d Center(const Point3d& p1, const Point3d& p2, const Point3d& p3);
    std::ostream& operator<<(std::ostream& s, const Point3d& p);
    inline Vec3d operator-(const Vec3d& p1, const Vec3d& v);
    inline Vec3d operator+(const Vec3d& p1, const Vec3d& v);
    inline Vec3d operator*(double scal, const Vec3d& v);
    inline double operator*(const Vec3d& v1, const Vec3d& v2);
    inline Vec3d Cross(const Vec3d& v1, const Vec3d& v2);
    inline Vec3d Cross(const Vec3d& v1, const Vec3d& v2);
    double Angle(const Vec3d& v1, const Vec3d& v2);
    std::ostream& operator<<(std::ostream& s, const Vec3d& v);
    void Transpose(Vec3d& v1, Vec3d& v2, Vec3d& v3);
    int SolveLinearSystem(const Vec3d& col1, const Vec3d& col2, const Vec3d& col3,
                          const Vec3d& rhs, Vec3d& sol);
    double Determinant(const Vec3d& col1, const Vec3d& col2, const Vec3d& col3);

    inline double Dist2(const Point3d& p1, const Point3d& p2);

    class Point3d
    {
     protected:
        double x[3];

     public:
        Point3d()
        {
            x[0] = x[1] = x[2] = 0;
        }

        Point3d(double ax, double ay, double az = 0.0)
        {
            x[0] = ax;
            x[1] = ay;
            x[2] = az;
        }

        Point3d(const Point3d& p2)
        {
            x[0] = p2.x[0];
            x[1] = p2.x[1];
            x[2] = p2.x[2];
        }

        explicit Point3d(const Point<3>& p2)
        {
            x[0] = p2[0];
            x[1] = p2[1];
            x[2] = p2[2];
        }

        Point3d& operator=(const Point3d& p2)
        {
            x[0] = p2.x[0];
            x[1] = p2.x[1];
            x[2] = p2.x[2];
            return *this;
        }

        int operator==(const Point3d& p) const
        {
            return (x[0] == p.x[0] && x[1] == p.x[1] && x[2] == p.x[2]);
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

        double operator[](size_t i) const
        {
            return x[i];
        }

        double& operator[](size_t i)
        {
            return x[i];
        }

        double& X(int i)
        {
            return x[i - 1];
        }

        double X(int i) const
        {
            return x[i - 1];
        }

        const Point3d& SetToMin(const Point3d& p2)
        {
            if (p2.x[0] < x[0]) x[0] = p2.x[0];
            if (p2.x[1] < x[1]) x[1] = p2.x[1];
            if (p2.x[2] < x[2]) x[2] = p2.x[2];
            return *this;
        }

        const Point3d& SetToMax(const Point3d& p2)
        {
            if (p2.x[0] > x[0]) x[0] = p2.x[0];
            if (p2.x[1] > x[1]) x[1] = p2.x[1];
            if (p2.x[2] > x[2]) x[2] = p2.x[2];
            return *this;
        }

        friend inline Vec3d operator-(const Point3d& p1, const Point3d& p2);
        friend inline Point3d operator-(const Point3d& p1, const Vec3d& v);
        friend inline Point3d operator+(const Point3d& p1, const Vec3d& v);
        inline Point3d& operator+=(const Vec3d& v);
        inline Point3d& operator-=(const Vec3d& v);

        friend inline double Dist(const Point3d& p1, const Point3d& p2)
        {
            return sqrt(
                (p1.x[0] - p2.x[0]) * (p1.x[0] - p2.x[0]) +
                (p1.x[1] - p2.x[1]) * (p1.x[1] - p2.x[1]) +
                (p1.x[2] - p2.x[2]) * (p1.x[2] - p2.x[2]));
        }

        inline friend double Dist2(const Point3d& p1, const Point3d& p2)
        {
            return ((p1.x[0] - p2.x[0]) * (p1.x[0] - p2.x[0]) +
                    (p1.x[1] - p2.x[1]) * (p1.x[1] - p2.x[1]) +
                    (p1.x[2] - p2.x[2]) * (p1.x[2] - p2.x[2]));
        }

        friend inline Point3d Center(const Point3d& p1, const Point3d& p2);
        friend inline Point3d Center(const Point3d& p1, const Point3d& p2, const Point3d& p3);
        friend std::ostream& operator<<(std::ostream& s, const Point3d& p);

        friend class Vec3d;

        friend class Box3d;

        operator Point<3>() const
        {
            return Point<3>(x[0], x[1], x[2]);
        }
    };

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

        explicit Vec3d(const Point3d& p1)
        {
            x[0] = p1.x[0];
            x[1] = p1.x[1];
            x[2] = p1.x[2];
        }

        explicit Vec3d(const Vec<3>& v2)
        {
            x[0] = v2[0];
            x[1] = v2[1];
            x[2] = v2[2];
        }

        Vec3d(const Point3d& p1, const Point3d& p2)
        {
            x[0] = p2.x[0] - p1.x[0];
            x[1] = p2.x[1] - p1.x[1];
            x[2] = p2.x[2] - p1.x[2];
        }

        Vec3d(const Vec3d& v2)
        {
            x[0] = v2.x[0];
            x[1] = v2.x[1];
            x[2] = v2.x[2];
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

        double& X(int i)
        {
            return x[i - 1];
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

        double X(int i) const
        {
            return x[i - 1];
        }

        double Length() const
        {
            return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        }

        double Length2() const
        {
            return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        }

        Vec3d& operator+=(const Vec3d& v2);
        Vec3d& operator-=(const Vec3d& v2);
        Vec3d& operator*=(double s);
        Vec3d& operator/=(double s);
        inline Vec3d& Add(double d, const Vec3d& v);

        friend inline Vec3d operator-(const Point3d& p1, const Point3d& p2);
        friend inline Point3d operator-(const Point3d& p1, const Vec3d& v);
        friend inline Point3d operator+(const Point3d& p1, const Vec3d& v);
        friend inline Vec3d operator-(const Vec3d& p1, const Vec3d& v);
        friend inline Vec3d operator+(const Vec3d& p1, const Vec3d& v);
        friend inline Vec3d operator*(double scal, const Vec3d& v);

        friend inline double operator*(const Vec3d& v1, const Vec3d& v2);
        friend inline Vec3d Cross(const Vec3d& v1, const Vec3d& v2);

        /// Returns one normal-vector to n
        void GetNormal(Vec3d& n) const;
        friend double Angle(const Vec3d& v1, const Vec3d& v2);

        void Normalize()
        {
            double len = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            if (len == 0) return;
            len = sqrt(len);
            x[0] /= len;
            x[1] /= len;
            x[2] /= len;
        }

        friend std::ostream& operator<<(std::ostream& s, const Vec3d& v);

        friend class Point3d;

        friend void Transpose(Vec3d& v1, Vec3d& v2, Vec3d& v3);
        friend int SolveLinearSystem(const Vec3d& col1,
                                     const Vec3d& col2,
                                     const Vec3d& col3,
                                     const Vec3d& rhs,
                                     Vec3d& sol);
        friend double Determinant(const Vec3d& col1,
                                  const Vec3d& col2,
                                  const Vec3d& col3);
    };

    inline Point3d Center(const Point3d& p1, const Point3d& p2)
    {
        return Point3d(
            0.5 * (p1.x[0] + p2.x[0]),
            0.5 * (p1.x[1] + p2.x[1]),
            0.5 * (p1.x[2] + p2.x[2]));
    }

    inline Point3d Center(const Point3d& p1, const Point3d& p2,
                          const Point3d& p3)
    {
        return Point3d(
            1.0 / 3.0 * (p1.x[0] + p2.x[0] + p3.x[0]),
            1.0 / 3.0 * (p1.x[1] + p2.x[1] + p3.x[1]),
            1.0 / 3.0 * (p1.x[2] + p2.x[2] + p3.x[2]));
    }

    inline Vec3d& Vec3d::operator+=(const Vec3d& v2)
    {
        x[0] += v2.X();
        x[1] += v2.Y();
        x[2] += v2.Z();
        return *this;
    }

    inline Vec3d& Vec3d::operator-=(const Vec3d& v2)
    {
        x[0] -= v2.X();
        x[1] -= v2.Y();
        x[2] -= v2.Z();
        return *this;
    }

    inline Vec3d& Vec3d::operator*=(double s)
    {
        x[0] *= s;
        x[1] *= s;
        x[2] *= s;
        return *this;
    }

    inline Vec3d& Vec3d::operator/=(double s)
    {
        if (s != 0) {
            x[0] /= s;
            x[1] /= s;
            x[2] /= s;
        }
        return *this;
    }

    inline Vec3d& Vec3d::Add(double d, const Vec3d& v)
    {
        x[0] += d * v.x[0];
        x[1] += d * v.x[1];
        x[2] += d * v.x[2];
        return *this;
    }

    inline Vec3d operator-(const Point3d& p1, const Point3d& p2)
    {
        return Vec3d(p1.x[0] - p2.x[0], p1.x[1] - p2.x[1], p1.x[2] - p2.x[2]);
    }

    inline Point3d operator-(const Point3d& p1, const Vec3d& v)
    {
        return Point3d(p1.x[0] - v.x[0], p1.x[1] - v.x[1], p1.x[2] - v.x[2]);
    }

    inline Point3d operator+(const Point3d& p1, const Vec3d& v)
    {
        return Point3d(p1.x[0] + v.x[0], p1.x[1] + v.x[1], p1.x[2] + v.x[2]);
    }

    inline Point3d& Point3d::operator+=(const Vec3d& v)
    {
        x[0] += v.x[0];
        x[1] += v.x[1];
        x[2] += v.x[2];
        return *this;
    }

    inline Point3d& Point3d::operator-=(const Vec3d& v)
    {
        x[0] -= v.x[0];
        x[1] -= v.x[1];
        x[2] -= v.x[2];
        return *this;
    }

    inline Vec3d operator-(const Vec3d& v1, const Vec3d& v2)
    {
        return Vec3d(v1.x[0] - v2.x[0], v1.x[1] - v2.x[1], v1.x[2] - v2.x[2]);
    }

    inline Vec3d operator+(const Vec3d& v1, const Vec3d& v2)
    {
        return Vec3d(v1.x[0] + v2.x[0], v1.x[1] + v2.x[1], v1.x[2] + v2.x[2]);
    }

    inline Vec3d operator*(double scal, const Vec3d& v)
    {
        return Vec3d(scal * v.x[0], scal * v.x[1], scal * v.x[2]);
    }

    inline double operator*(const Vec3d& v1, const Vec3d& v2)
    {
        return v1.x[0] * v2.x[0] + v1.x[1] * v2.x[1] + v1.x[2] * v2.x[2];
    }


    inline Vec3d operator*(const Mat<3, 3>& m, const Vec3d& v)
    {
        Vec3d res;
        for (int i = 0; i < 3; i++) {
            res.X(i + 1) = 0;
            for (int j = 0; j < 3; j++)
                res.X(i + 1) += m(i, j) * v.X(j + 1);
        }
        return res;
    }

    inline Vec3d Cross(const Vec3d& v1, const Vec3d& v2)
    {
        return Vec3d(
            v1.Y() * v2.Z() - v1.Z() * v2.Y(),
            v1.Z() * v2.X() - v1.X() * v2.Z(),
            v1.X() * v2.Y() - v1.Y() * v2.X());
    }

    inline double Determinant(
        const Vec3d& col1,
        const Vec3d& col2,
        const Vec3d& col3)
    {
        return
            col1.x[0] * (col2.x[1] * col3.x[2] - col2.x[2] * col3.x[1]) +
            col1.x[1] * (col2.x[2] * col3.x[0] - col2.x[0] * col3.x[2]) +
            col1.x[2] * (col2.x[0] * col3.x[1] - col2.x[1] * col3.x[0]);
    }

    class Box3d
    {
     protected:
        double minx[3], maxx[3];

     public:
        Box3d() { }
        Box3d(double aminx, double amaxx,
              double aminy, double amaxy,
              double aminz, double amaxz);
        Box3d(const Box3d& b2);
        Box3d(const Point3d& p1, const Point3d& p2);
        Box3d(const Box<3>& b2);

        Point3d PMin() const
        {
            return Point3d(minx[0], minx[1], minx[2]);
        }

        Point3d PMax() const
        {
            return Point3d(maxx[0], maxx[1], maxx[2]);
        }

        inline void SetPoint(const Point3d& p)
        {
            minx[0] = maxx[0] = p.X();
            minx[1] = maxx[1] = p.Y();
            minx[2] = maxx[2] = p.Z();
        }

        inline void AddPoint(const Point3d& p)
        {
            if (p.x[0] < minx[0]) minx[0] = p.x[0];
            if (p.x[0] > maxx[0]) maxx[0] = p.x[0];
            if (p.x[1] < minx[1]) minx[1] = p.x[1];
            if (p.x[1] > maxx[1]) maxx[1] = p.x[1];
            if (p.x[2] < minx[2]) minx[2] = p.x[2];
            if (p.x[2] > maxx[2]) maxx[2] = p.x[2];
        }

        const Box3d& operator+=(const Box3d& b);
    };

}  // namespace meshit


#endif
