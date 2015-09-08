#ifndef FILE_DENSEMAT
#define FILE_DENSEMAT

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
        int height;
        int width;
        double * data;

      public:
        DLL_HEADER DenseMatrix();
        DLL_HEADER DenseMatrix(int h, int w = 0);
        DLL_HEADER DenseMatrix(const DenseMatrix & m2);
        DLL_HEADER ~DenseMatrix();

        DLL_HEADER void SetSize(int h, int w = 0);

        int Height() const
        {
            return height;
        }

        int Width() const
        {
            return width;
        }

        double & operator()(int i, int j)
        {
            return data[i * width + j];
        }

        double operator()(int i, int j) const
        {
            return data[i * width + j];
        }

        double & operator()(int i)
        {
            return data[i];
        }

        double operator()(int i) const
        {
            return data[i];
        }

        DLL_HEADER DenseMatrix & operator=(const DenseMatrix & m2);
        DLL_HEADER DenseMatrix & operator+=(const DenseMatrix & m2);
        DLL_HEADER DenseMatrix & operator-=(const DenseMatrix & m2);

        DLL_HEADER DenseMatrix & operator=(double v);
        DLL_HEADER DenseMatrix & operator*=(double v);

        DLL_HEADER void Mult(const FlatVector & v, FlatVector & prod) const
        {
            const double * mp, * sp;
            double * dp;

            mp = data;
            dp = &prod(0);
            for (int i = 0; i < height; i++) {
                double sum = 0;
                sp = &v(0);

                for (int j = 0; j < width; j++) {
                    sum += *mp * *sp;
                    mp++;
                    sp++;
                }

                *dp = sum;
                dp++;
            }
        }

        DLL_HEADER void MultTrans(const Vector & v, Vector & prod) const;
        DLL_HEADER void Residuum(const Vector & x, const Vector & b, Vector & res) const;
        DLL_HEADER double Det() const;

        friend DenseMatrix operator*(const DenseMatrix & m1, const DenseMatrix & m2);
        friend DenseMatrix operator+(const DenseMatrix & m1, const DenseMatrix & m2);

        friend void Transpose(const DenseMatrix & m1, DenseMatrix & m2);
        friend void Mult(const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3);
        friend void CalcAAt(const DenseMatrix & a, DenseMatrix & m2);
        friend void CalcABt(const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
        friend void CalcAtB(const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
        DLL_HEADER void Solve(const Vector & b, Vector & x) const;
        void SolveDestroy(const Vector & b, Vector & x);

        const double & Get(int i, int j) const
        {
            return data[(i - 1) * width + j - 1];
        }

        const double & Get(int i) const
        {
            return data[i - 1];
        }

        void Set(int i, int j, double v)
        {
            data[(i - 1) * width + j - 1] = v;
        }

        double & Elem(int i, int j)
        {
            return data[(i - 1) * width + j - 1];
        }

        const double & ConstElem(int i, int j) const
        {
            return data[(i - 1) * width + j - 1];
        }
    };

    extern std::ostream & operator<<(std::ostream & ost, const DenseMatrix & m);

    template <int WIDTH>
    class MatrixFixWidth
    {
      protected:

        int height;
        double * data;
        bool ownmem;

      public:

        MatrixFixWidth()
        {
            height = 0;
            data = 0;
            ownmem = false;
        }

        MatrixFixWidth(int h)
        {
            height = h;
            data = new double[WIDTH * height];
            ownmem = true;
        }

        MatrixFixWidth(int h, double * adata)
        {
            height = h;
            data = adata;
            ownmem = false;
        }

        ~MatrixFixWidth()
        {
            if (ownmem) delete [] data;
        }

        void SetSize(int h)
        {
            if (h != height) {
                if (ownmem) delete data;
                height = h;
                data = new double[WIDTH * height];
                ownmem = true;
            }
        }

        int Height() const
        {
            return height;
        }

        int Width() const
        {
            return WIDTH;
        }

        MatrixFixWidth & operator=(double v)
        {
            for (int i = 0; i < height * WIDTH; i++)
                data[i] = v;
            return *this;
        }

        void Mult(const FlatVector & v, FlatVector & prod) const
        {
            double sum;
            const double * mp, * sp;
            double * dp;

            mp = data;
            dp = &prod[0];
            for (int i = 0; i < height; i++) {
                sum = 0;
                sp = &v[0];

                for (int j = 0; j < WIDTH; j++) {
                    sum += *mp * *sp;
                    mp++;
                    sp++;
                }

                *dp = sum;
                dp++;
            }
        }

        double & operator()(int i, int j)
        {
            return data[i * WIDTH + j];
        }

        const double & operator()(int i, int j) const
        {
            return data[i * WIDTH + j];
        }

        MatrixFixWidth & operator*=(double v)
        {
            if (data)
                for (int i = 0; i < height * WIDTH; i++)
                    data[i] *= v;
            return *this;
        }

        const double & Get(int i, int j) const
        {
            return data[(i - 1) * WIDTH + j - 1];
        }

        const double & Get(int i) const
        {
            return data[i - 1];
        }

        void Set(int i, int j, double v)
        {
            data[(i - 1) * WIDTH + j - 1] = v;
        }

        double & Elem(int i, int j)
        {
            return data[(i - 1) * WIDTH + j - 1];
        }

        const double & ConstElem(int i, int j) const
        {
            return data[(i - 1) * WIDTH + j - 1];
        }
    };

    template <int WIDTH>
    extern std::ostream & operator<<(std::ostream & ost, const MatrixFixWidth<WIDTH> & m)
    {
        for (int i = 0; i < m.Height(); i++) {
            for (int j = 0; j < m.Width(); j++)
                ost << m.Get(i + 1, j + 1) << " ";
            ost << std::endl;
        }
        return ost;
    };

    extern DLL_HEADER void CalcAtA(const DenseMatrix & a, DenseMatrix & m2);
    extern DLL_HEADER void CalcInverse(const DenseMatrix & m1, DenseMatrix & m2);
}

#endif