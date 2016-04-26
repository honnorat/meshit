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

#include "densemat.hpp"

namespace meshit {

DenseMatrix::DenseMatrix(size_t h, size_t w)
{
    if (w == 0) w = h;

    width = w;
    height = h;

    if (h * w > 0) {
        data = new double[h * w];
    } else {
        data = nullptr;
    }
    for (size_t i = 0; i < (h * w); i++) {
        data[i] = 0.0;
    }
}

DenseMatrix::DenseMatrix(const DenseMatrix& m2)
{
    data = nullptr;
    height = width = 0;
    SetSize(m2.Height(), m2.Width());
    memcpy(data, m2.data, sizeof(double) * Height() * Width());
}

DenseMatrix::~DenseMatrix()
{
    delete[] data;
}

void DenseMatrix::SetSize(size_t h, size_t w)
{
    if (w == 0) w = h;

    if (height == h && width == w) return;

    height = h;
    width = w;
    delete[] data;

    if (h * w > 0) {
        data = new double[h * w];
    } else {
        data = nullptr;
    }
}

DenseMatrix& DenseMatrix::operator=(const DenseMatrix& m2)
{
    SetSize(m2.Height(), m2.Width());

    if (data) memcpy(data, m2.data, sizeof(double) * m2.Height() * m2.Width());
    return *this;
}

DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& m2)
{
    if (height != m2.Height() || width != m2.Width()) {
        std::cerr << "DenseMatrix::Operator+=: Sizes don't fit" << std::endl;
        return *this;
    }

    if (data) {
        size_t wh = width * height;
        double* p = data;
        double* q = m2.data;

        for (size_t i = 0; i < wh; i++) {
            *p += *q;
            p++;
            q++;
        }
    } else {
        std::cerr << "DenseMatrix::Operator+=: Matrix not allocated" << std::endl;
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator-=(const DenseMatrix& m2)
{
    if (height != m2.Height() || width != m2.Width()) {
        std::cerr << "DenseMatrix::Operator-=: Sizes don't fit" << std::endl;
        return *this;
    }

    if (data) {
        size_t wh = width * height;
        double* p = data;
        double* q = m2.data;

        for (size_t i = 0; i < wh; i++) {
            *p -= *q;
            p++;
            q++;
        }
    } else {
        std::cerr << "DenseMatrix::Operator-=: Matrix not allocated" << std::endl;
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator=(double v)
{
    if (data) {
        size_t wh = width * height;
        double* p = data;
        for (size_t i = 0; i < wh; i++) {
            *p = v;
            p++;
        }
    }
    return *this;
}

DenseMatrix& DenseMatrix::operator*=(double v)
{
    if (data) {
        size_t wh = width * height;
        double* p = data;
        for (size_t i = 0; i < wh; i++) {
            *p *= v;
            p++;
        }
    }
    return *this;
}

void DenseMatrix::Mult(const FlatVector& v, FlatVector& prod) const
{
    const double* mp = data;
    double* dp = &prod[0];
    for (size_t i = 0; i < height; i++) {
        double sum = 0;
        const double* sp = &v[0];

        for (size_t j = 0; j < width; j++) {
            sum += *mp * *sp;
            mp++;
            sp++;
        }

        *dp = sum;
        dp++;
    }
}

}  // namespace meshit
