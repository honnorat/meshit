#ifndef MESHIT_OBJECTS_HPP
#define MESHIT_OBJECTS_HPP
/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2015-2016 Marc Honnorat <marc.honnorat@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library in the file LICENSE.LGPL; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

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
