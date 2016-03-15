#include "densemat.hpp"

#include "../general/array.hpp"

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

        if (height == h && width == w)
            return;

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

    double DenseMatrix::Det() const
    {
        if (width != height) {
            std::cerr << "DenseMatrix :: Det: width != height" << std::endl;
            return 0.0;
        }

        switch (width) {
            case 1:
                return data[0];
            case 2:
                return data[0] * data[3] - data[1] * data[2];
            case 3:
                return data[0] * data[4] * data[8]
                       + data[1] * data[5] * data[6]
                       + data[2] * data[3] * data[7]
                       - data[0] * data[5] * data[7]
                       - data[1] * data[3] * data[8]
                       - data[2] * data[4] * data[6];
            default: {
                std::cerr << "Matrix :: Det:  general size not implemented (size=" << width << ")" << std::endl;
                return 0.0;
            }
        }
    }

    void CalcABt(const DenseMatrix& a, const DenseMatrix& b, DenseMatrix& m2)
    {
        size_t n1 = a.Height();
        size_t n2 = a.Width();
        size_t n3 = b.Height();

        if (m2.Height() != n1 || m2.Width() != n3 || b.Width() != n2) {
            std::cerr << "CalcABt: sizes don't fit" << std::endl;
            return;
        }

        double* pm2 = m2.DataP();
        const double* pa1 = a.DataP();

        for (size_t i = 0; i < n1; i++) {
            const double* pb = b.DataP();
            for (size_t j = 0; j < n3; j++) {
                double sum = 0.0;
                const double* pa = pa1;
                for (size_t k = 0; k < n2; k++) {
                    sum += (*pa) * (*pb);
                    pa++;
                    pb++;
                }

                *pm2 = sum;
                pm2++;
            }
            pa1 += n2;
        }
    }

}  // namespace meshit
