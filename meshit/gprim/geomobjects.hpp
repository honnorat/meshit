#ifndef FILE_OBJECTS_HPP
#define FILE_OBJECTS_HPP

/* *************************************************************************/
/* File:   geomobjects.hpp                                                 */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace meshit {

class Mat2x2
{
 public:
    Mat2x2() { }

    Mat2x2(const Mat2x2& b)
    {
        x[0] = b.x[0];
        x[1] = b.x[1];
        x[2] = b.x[2];
        x[3] = b.x[3];
    }

    Mat2x2(double x0, double x1, double x2, double x3)
    {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
        x[3] = x3;
    }

    double& operator()(size_t i, size_t j) { return x[i * 2 + j]; }
    double operator()(size_t i, size_t j) const { return x[i * 2 + j]; }

 protected:
    double x[4];
};

class Mat3x3
{
 public:
    Mat3x3() { }

    Mat3x3(const Mat3x3& b)
    {
        for (size_t i = 0; i < 9; i++) x[i] = b.x[i];
    }

    Mat3x3& operator=(double s)
    {
        for (size_t i = 0; i < 9; i++) x[i] = s;
        return *this;
    }

    Mat3x3& operator=(const Mat3x3& b)
    {
        for (size_t i = 0; i < 9; i++) x[i] = b.x[i];
        return *this;
    }

    double& operator()(size_t i, size_t j) { return x[i * 3 + j]; }

    const double& operator()(size_t i, size_t j) const { return x[i * 3 + j]; }

    void CalcInverse(Mat3x3& inv) const;
    double Det() const;

 protected:
    double x[9];
};

}  // namespace meshit

#endif
