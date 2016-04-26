#ifndef MESHIT_VECTOR_HPP
#define MESHIT_VECTOR_HPP
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

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>

namespace meshit {

class FlatVector
{
 public:
    FlatVector() { }

    explicit FlatVector(size_t as, double* adata = nullptr)
        : size_{as}, data{adata} { }

    FlatVector& operator=(const FlatVector& v)
    {
        memcpy(data, v.data, size_ * sizeof(double));
        return *this;
    }

    FlatVector& operator=(double scal)
    {
        for (size_t i = 0; i < size_; i++) data[i] = scal;
        return *this;
    }

    size_t Size() const { return size_; }
    double& operator[](size_t i) { return data[i]; }
    const double& operator[](size_t i) const { return data[i]; }

    FlatVector& operator*=(double scal)
    {
        for (size_t i = 0; i < size_; i++) data[i] *= scal;
        return *this;
    }

    FlatVector& Add(double scal, const FlatVector& v2)
    {
        for (size_t i = 0; i < size_; i++) {
            data[i] += scal * v2[i];
        }
        return *this;
    }

    FlatVector& Set2(double scal1, const FlatVector& v1, double scal2, const FlatVector& v2)
    {
        for (size_t i = 0; i < size_; i++) {
            data[i] = scal1 * v1[i] + scal2 * v2[i];
        }
        return *this;
    }

    friend double operator*(const FlatVector& v1, const FlatVector& v2);

 protected:
    size_t size_;
    double* data;
};


class Vector : public FlatVector
{
    bool ownmem;

 public:
    Vector()
        : FlatVector{0, nullptr}, ownmem{false} { }

    explicit Vector(size_t as)
        : FlatVector{as}
    {
        data = new double[size_];
        ownmem = true;
    }

    Vector(size_t as, double* mem)
        : FlatVector{as, mem}, ownmem{false} { }

    ~Vector()
    {
        if (ownmem) delete[] data;
    }

    Vector& operator=(const FlatVector& v)
    {
        memcpy(data, &v[0], size_ * sizeof(double));
        return *this;
    }

    Vector& operator=(double scal)
    {
        for (size_t i = 0; i < size_; i++) data[i] = scal;
        return *this;
    }

    void SetSize(size_t as)
    {
        if (size_ != as) {
            size_ = as;
            if (ownmem) delete[] data;
            data = new double[size_];
            ownmem = true;
        }
    }
};

inline double operator*(const FlatVector& v1, const FlatVector& v2)
{
    double sum = 0;
    for (size_t i = 0; i < v1.Size(); i++) {
        sum += v1.data[i] * v2.data[i];
    }
    return sum;
}

inline std::ostream& operator<<(std::ostream& ost, const FlatVector& v)
{
    for (size_t i = 0; i < v.Size(); i++) {
        ost << " " << std::setw(7) << v[i];
    }
    return ost;
}

}  // namespace meshit

#endif
