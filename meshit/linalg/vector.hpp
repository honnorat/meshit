#ifndef FILE_VECTOR_HPP
#define FILE_VECTOR_HPP

/* *************************************************************************/
/* File:   vector.hpp                                                      */
/* Author: Joachim Schoeberl                                               */
/* Date:   01. Oct. 94                                                     */
/* *************************************************************************/
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>

namespace meshit
{
    class FlatVector
    {
     protected:
        size_t size_;
        double* data;

     public:
        FlatVector() { }

        explicit FlatVector(size_t as, double* adata = nullptr)
            : size_{as}, data{adata} { }

        size_t Size() const { return size_; }

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

        FlatVector& Set2(double scal1, const FlatVector& v1,
                         double scal2, const FlatVector& v2)
        {
            for (size_t i = 0; i < size_; i++) {
                data[i] = scal1 * v1[i] + scal2 * v2[i];
            }
            return *this;
        }

        friend double operator*(const FlatVector& v1, const FlatVector& v2);
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


