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
 public:
    DenseMatrix()
        : height{0}, width{0}, data{nullptr} { }

    explicit DenseMatrix(size_t h, size_t w = 0);
    DenseMatrix(const DenseMatrix& m2);
    ~DenseMatrix();

    void SetSize(size_t h, size_t w = 0);

    size_t Height() const { return height; }
    size_t Width() const { return width; }

    double& operator()(size_t i, size_t j) { return data[i * width + j]; }
    double operator()(size_t i, size_t j) const { return data[i * width + j]; }
    double& Elem(size_t i, size_t j) { return data[i * width + j]; }

    DenseMatrix& operator=(const DenseMatrix& m2);
    DenseMatrix& operator+=(const DenseMatrix& m2);
    DenseMatrix& operator-=(const DenseMatrix& m2);

    DenseMatrix& operator=(double v);
    DenseMatrix& operator*=(double v);

    void Mult(const FlatVector& v, FlatVector& prod) const;

 protected:
    size_t height;
    size_t width;
    double* data;
};

template<size_t WIDTH>
class MatrixFixWidth
{
 protected:
    size_t height;
    bool ownmem;
    double* data;

 public:
    MatrixFixWidth()
        : height{0}, ownmem{false}, data{nullptr} { }

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
            if (ownmem) delete[] data;
            height = h;
            data = new double[WIDTH * height];
            ownmem = true;
        }
    }

    size_t Height() const { return height; }
    size_t Width() const { return WIDTH; }

    MatrixFixWidth& operator=(double v)
    {
        for (size_t i = 0; i < height * WIDTH; i++) {
            data[i] = v;
        }
        return *this;
    }

    double& operator()(size_t i, size_t j) { return data[i * WIDTH + j]; }
    const double& operator()(size_t i, size_t j) const { return data[i * WIDTH + j]; }

    MatrixFixWidth& operator*=(double v)
    {
        if (data) {
            for (size_t i = 0; i < height * WIDTH; i++) {
                data[i] *= v;
            }
        }
        return *this;
    }
};

template<size_t WIDTH>
extern std::ostream& operator<<(std::ostream& ost, const MatrixFixWidth<WIDTH>& m)
{
    for (size_t i = 0; i < m.Height(); i++) {
        for (size_t j = 0; j < m.Width(); j++) {
            ost << m(i, j) << " ";
        }
        ost << std::endl;
    }
    return ost;
};

}  // namespace meshit

#endif
