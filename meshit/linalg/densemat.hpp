#ifndef FILE_DENSEMAT_HPP
#define FILE_DENSEMAT_HPP

/**************************************************************************/
/* File:   densemat.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/**************************************************************************/

/** 
    Data type dense matrix
 */

#include "vector.hpp"

namespace meshit {

    class DenseMatrix
    {
     protected:
        size_t height;
        size_t width;
        double* data;

     public:
        DenseMatrix()
                : height{0}, width{0}, data{nullptr} { }

        explicit DenseMatrix(size_t h, size_t w = 0);
        DenseMatrix(const DenseMatrix& m2);
        ~DenseMatrix();

        void SetSize(size_t h, size_t w = 0);

        size_t Height() const
        {
            return height;
        }

        size_t Width() const
        {
            return width;
        }

        double& operator()(size_t i, size_t j)
        {
            return data[i * width + j];
        }

        double operator()(size_t i, size_t j) const
        {
            return data[i * width + j];
        }

        double& operator()(size_t i)
        {
            return data[i];
        }

        double operator()(size_t i) const
        {
            return data[i];
        }

        DenseMatrix& operator=(const DenseMatrix& m2);
        DenseMatrix& operator+=(const DenseMatrix& m2);
        DenseMatrix& operator-=(const DenseMatrix& m2);

        DenseMatrix& operator=(double v);
        DenseMatrix& operator*=(double v);

        void Mult(const FlatVector& v, FlatVector& prod) const
        {
            const double* mp = data;
            double* dp = &prod(0);
            for (size_t i = 0; i < height; i++) {
                double sum = 0;
                const double* sp = &v(0);

                for (size_t j = 0; j < width; j++) {
                    sum += *mp * *sp;
                    mp++;
                    sp++;
                }

                *dp = sum;
                dp++;
            }
        }

        void MultTrans(const Vector& v, Vector& prod) const;
        void Residuum(const Vector& x, const Vector& b, Vector& res) const;
        double Det() const;

        friend DenseMatrix operator*(const DenseMatrix& m1, const DenseMatrix& m2);
        friend DenseMatrix operator+(const DenseMatrix& m1, const DenseMatrix& m2);

        friend void Transpose(const DenseMatrix& m1, DenseMatrix& m2);
        friend void Mult(const DenseMatrix& m1, const DenseMatrix& m2, DenseMatrix& m3);
        friend void CalcAAt(const DenseMatrix& a, DenseMatrix& m2);
        friend void CalcABt(const DenseMatrix& a, const DenseMatrix& b, DenseMatrix& m2);
        friend void CalcAtB(const DenseMatrix& a, const DenseMatrix& b, DenseMatrix& m2);
        void Solve(const Vector& b, Vector& x) const;
        void SolveDestroy(const Vector& b, Vector& x);

        const double& Get(size_t i, size_t j) const
        {
            return data[(i - 1) * width + j - 1];
        }

        const double& Get(size_t i) const
        {
            return data[i - 1];
        }

        void Set(size_t i, size_t j, double v)
        {
            data[(i - 1) * width + j - 1] = v;
        }

        double& Elem(size_t i, size_t j)
        {
            return data[(i - 1) * width + j - 1];
        }

        const double& ConstElem(size_t i, size_t j) const
        {
            return data[(i - 1) * width + j - 1];
        }
    };

    extern std::ostream& operator<<(std::ostream& ost, const DenseMatrix& m);

    template<size_t WIDTH>
    class MatrixFixWidth
    {
     protected:
        size_t height;
        double* data;
        bool ownmem;

     public:
        MatrixFixWidth()
                : height{0}, ownmem{false}
        {
            data = nullptr;
        }

        explicit MatrixFixWidth(size_t h)
                : height{h}, ownmem{true}
        {
            data = new double[WIDTH * height];
        }

        explicit MatrixFixWidth(size_t h, double* adata)
                : height{h}, ownmem{false}, data{adata} { }

        ~MatrixFixWidth()
        {
            if (ownmem) delete[] data;
        }

        void SetSize(size_t h)
        {
            if (h != height) {
                if (ownmem) delete data;
                height = h;
                data = new double[WIDTH * height];
                ownmem = true;
            }
        }

        size_t Height() const
        {
            return height;
        }

        size_t Width() const
        {
            return WIDTH;
        }

        MatrixFixWidth& operator=(double v)
        {
            for (size_t i = 0; i < height * WIDTH; i++) {
                data[i] = v;
            }
            return *this;
        }

        void Mult(const FlatVector& v, FlatVector& prod) const
        {
            const double* mp = data;
            double* dp = &prod[0];

            for (size_t i = 0; i < height; i++) {
                double sum = 0;
                const double* sp = &v[0];

                for (size_t j = 0; j < WIDTH; j++) {
                    sum += *mp * *sp;
                    mp++;
                    sp++;
                }

                *dp = sum;
                dp++;
            }
        }

        double& operator()(size_t i, size_t j)
        {
            return data[i * WIDTH + j];
        }

        const double& operator()(size_t i, size_t j) const
        {
            return data[i * WIDTH + j];
        }

        MatrixFixWidth& operator*=(double v)
        {
            if (data) {
                for (size_t i = 0; i < height * WIDTH; i++) {
                    data[i] *= v;
                }
            }
            return *this;
        }

        const double& Get(size_t i, size_t j) const
        {
            return data[(i - 1) * WIDTH + j - 1];
        }

        const double& Get(size_t i) const
        {
            return data[i - 1];
        }

        void Set(size_t i, size_t j, double v)
        {
            data[(i - 1) * WIDTH + j - 1] = v;
        }

        double& Elem(size_t i, size_t j)
        {
            return data[(i - 1) * WIDTH + j - 1];
        }
    };

    template<size_t WIDTH>
    extern std::ostream& operator<<(std::ostream& ost, const MatrixFixWidth<WIDTH>& m)
    {
        for (size_t i = 0; i < m.Height(); i++) {
            for (size_t j = 0; j < m.Width(); j++) {
                ost << m.Get(i + 1, j + 1) << " ";
            }
            ost << std::endl;
        }
        return ost;
    };

    extern void CalcAtA(const DenseMatrix& a, DenseMatrix& m2);

}  // namespace meshit

#endif
