#ifndef FILE_VECTOR
#define FILE_VECTOR

/* *************************************************************************/
/* File:   vector.hpp                                                      */
/* Author: Joachim Schoeberl                                               */
/* Date:   01. Oct. 94                                                     */
/* *************************************************************************/
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>

namespace meshit {

    class FlatVector
    {
     protected:
        size_t s;
        double* data;

     public:
        FlatVector() { }

        explicit FlatVector(size_t as, double* adata = nullptr)
                : s{as}, data{adata} { }

        size_t Size() const
        {
            return s;
        }

        FlatVector& operator=(const FlatVector& v)
        {
            memcpy(data, v.data, s * sizeof(double));
            return *this;
        }

        FlatVector& operator=(double scal)
        {
            for (size_t i = 0; i < s; i++) data[i] = scal;
            return *this;
        }

        double& operator[](size_t i)
        {
            return data[i];
        }

        const double& operator[](size_t i) const
        {
            return data[i];
        }

        double& operator()(size_t i)
        {
            return data[i];
        }

        const double& operator()(size_t i) const
        {
            return data[i];
        }

        FlatVector& operator*=(double scal)
        {
            for (size_t i = 0; i < s; i++) data[i] *= scal;
            return *this;
        }

        FlatVector& Add(double scal, const FlatVector& v2)
        {
            for (size_t i = 0; i < s; i++) {
                data[i] += scal * v2[i];
            }
            return *this;
        }

        FlatVector& Set(double scal, const FlatVector& v2)
        {
            for (size_t i = 0; i < s; i++) {
                data[i] = scal * v2[i];
            }
            return *this;
        }

        FlatVector& Set2(double scal1, const FlatVector& v1,
                         double scal2, const FlatVector& v2)
        {
            for (size_t i = 0; i < s; i++) {
                data[i] = scal1 * v1[i] + scal2 * v2[i];
            }
            return *this;
        }

        double L2Norm() const
        {
            double sum = 0;
            for (size_t i = 0; i < s; i++) {
                sum += data[i] * data[i];
            }
            return sqrt(sum);
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
            data = new double[s];
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
            memcpy(data, &v(0), s * sizeof(double));
            return *this;
        }

        Vector& operator=(double scal)
        {
            for (size_t i = 0; i < s; i++) data[i] = scal;
            return *this;
        }

        void SetSize(size_t as)
        {
            if (s != as) {
                s = as;
                if (ownmem) delete[] data;
                data = new double[s];
                ownmem = true;
            }
        }
    };

    template<size_t S>
    class VectorMem : public Vector
    {
        double mem[S];

     public:
        VectorMem() : Vector(S, &mem[0]) { }

        VectorMem& operator=(const FlatVector& v)
        {
            memcpy(data, &v(0), S * sizeof(double));
            return *this;
        }

        VectorMem& operator=(double scal)
        {
            for (size_t i = 0; i < S; i++) data[i] = scal;
            return *this;
        }
    };

    inline double operator*(const FlatVector& v1, const FlatVector& v2)
    {
        double sum = 0;
        for (size_t i = 0; i < v1.s; i++) {
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


