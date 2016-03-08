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

    void CalcAAt(const DenseMatrix& a, DenseMatrix& m2)
    {
        size_t n1 = a.Height();
        size_t n2 = a.Width();

        if (m2.Height() != n1 || m2.Width() != n1) {
            std::cerr << "CalcAAt: sizes don't fit" << std::endl;
            return;
        }

        for (size_t i = 1; i <= n1; i++) {
            double sum = 0;
            const double* p = &a.ConstElem(i, 1);
            for (size_t k = 1; k <= n2; k++) {
                sum += *p * *p;
                p++;
            }
            m2.Set(i, i, sum);

            const double* p0 = &a.ConstElem(i, 1);
            const double* q = a.data;
            for (size_t j = 1; j < i; j++) {
                sum = 0;
                p = p0;

                for (size_t k = 1; k <= n2; k++) {
                    sum += *p * *q;
                    p++;
                    q++;
                }
                m2.Set(i, j, sum);
                m2.Set(j, i, sum);
            }
        }
    }

    void CalcAtA(const DenseMatrix& a, DenseMatrix& m2)
    {
        size_t n1 = a.Height();
        size_t n2 = a.Width();

        if (m2.Height() != n2 || m2.Width() != n2) {
            std::cerr << "CalcAtA: sizes don't fit" << std::endl;
            return;
        }

        for (size_t i = 1; i <= n2; i++) {
            for (size_t j = 1; j <= n2; j++) {
                double sum = 0;
                for (size_t k = 1; k <= n1; k++) {
                    sum += a.Get(k, i) * a.Get(k, j);
                }
                m2.Elem(i, j) = sum;
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

        double* pm2 = &m2.Elem(1, 1);
        const double* pa1 = &a.Get(1, 1);

        for (size_t i = 1; i <= n1; i++) {
            const double* pb = &b.Get(1, 1);
            for (size_t j = 1; j <= n3; j++) {
                double sum = 0.0;
                const double* pa = pa1;

                for (size_t k = 1; k <= n2; k++) {
                    sum += *pa * *pb;
                    pa++;
                    pb++;
                }

                *pm2 = sum;
                pm2++;
            }
            pa1 += n2;
        }
    }

    void CalcAtB(const DenseMatrix& a, const DenseMatrix& b, DenseMatrix& m2)
    {
        size_t n1 = a.Height();
        size_t n2 = a.Width();
        size_t n3 = b.Width();

        if (m2.Height() != n2 || m2.Width() != n3 || b.Height() != n1) {
            std::cerr << "CalcAtB: sizes don't fit" << std::endl;
            return;
        }

        for (size_t i = 1; i <= n2 * n3; i++)
            m2.data[i - 1] = 0;

        for (size_t i = 1; i <= n1; i++)
            for (size_t j = 1; j <= n2; j++) {
                const double va = a.Get(i, j);
                double* pm2 = &m2.Elem(j, 1);
                const double* pb = &b.Get(i, 1);

                for (size_t k = 1; k <= n3; ++k, ++pm2, ++pb)
                    *pm2 += va * *pb;
            }
    }

    DenseMatrix operator*(const DenseMatrix& m1, const DenseMatrix& m2)
    {
        DenseMatrix temp(m1.Height(), m2.Width());

        if (m1.Width() != m2.Height()) {
            std::cerr << "DenseMatrix :: operator*: Matrix Size does not fit" << std::endl;
        } else if (temp.Height() != m1.Height()) {
            std::cerr << "DenseMatrix :: operator*: temp not allocated" << std::endl;
        } else {
            Mult(m1, m2, temp);
        }
        return temp;
    }

    void Mult(const DenseMatrix& m1, const DenseMatrix& m2, DenseMatrix& m3)
    {
        if (m1.Width() != m2.Height() || m1.Height() != m3.Height() ||
            m2.Width() != m3.Width()) {
            std::cerr << "DenseMatrix :: Mult: Matrix Size does not fit" << std::endl;
            std::cerr << "m1: " << m1.Height() << " x " << m1.Width() << std::endl;
            std::cerr << "m2: " << m2.Height() << " x " << m2.Width() << std::endl;
            std::cerr << "m3: " << m3.Height() << " x " << m3.Width() << std::endl;
            return;
        } else {
            size_t n1 = m1.Height();
            size_t n2 = m2.Width();
            size_t n3 = m1.Width();

            double* p3 = m3.data;
            double* p1s = m1.data;
            double* p2sn = m2.data + n2;
            double* p1snn = p1s + n1 * n3;

            while (p1s != p1snn) {
                double* p1sn = p1s + n3;
                double* p2s = m2.data;

                while (p2s != p2sn) {
                    double sum = 0.0;
                    double* p1 = p1s;
                    double* p2 = p2s;
                    p2s++;

                    while (p1 != p1sn) {
                        sum += *p1 * *p2;
                        p1++;
                        p2 += n2;
                    }
                    *p3++ = sum;
                }
                p1s = p1sn;
            }
        }
    }

    DenseMatrix operator+(const DenseMatrix& m1, const DenseMatrix& m2)
    {
        DenseMatrix temp(m1.Height(), m1.Width());

        if (m1.Width() != m2.Width() || m1.Height() != m2.Height()) {
            std::cerr << "BaseMatrix :: operator+: Matrix Size does not fit" << std::endl;
        } else if (temp.Height() != m1.Height()) {
            std::cerr << "BaseMatrix :: operator+: temp not allocated" << std::endl;
        } else {
            for (size_t i = 1; i <= m1.Height(); i++) {
                for (size_t j = 1; j <= m1.Width(); j++) {
                    temp.Set(i, j, m1.Get(i, j) + m2.Get(i, j));
                }
            }
        }
        return temp;
    }

    void Transpose(const DenseMatrix& m1, DenseMatrix& m2)
    {
        size_t w = m1.Width();
        size_t h = m1.Height();

        m2.SetSize(w, h);

        double* pm2 = &m2.Elem(1, 1);
        for (size_t j = 1; j <= w; j++) {
            const double* pm1 = &m1.Get(1, j);
            for (size_t i = 1; i <= h; i++) {
                *pm2 = *pm1;
                pm2++;
                pm1 += w;
            }
        }
    }

    void DenseMatrix::MultTrans(const Vector& v, Vector& prod) const
    {
        if (prod.Size() != width)
            prod.SetSize(width);

        const double* pmat = &Get(1, 1);
        const double* pv = &v(0);

        prod = 0;

        for (size_t i = 1; i <= height; i++) {
            double val = *pv;
            ++pv;

            double* pprod = &prod(0);

            for (size_t j = width - 1; j >= 0; --j, ++pmat, ++pprod) {
                *pprod += val * *pmat;
            }
        }
    }

    void DenseMatrix::Residuum(const Vector& x, const Vector& b, Vector& res) const
    {
        res.SetSize(height);

        if (width != x.Size() || height != b.Size()) {
            std::cerr << "\nMatrix and Vector don't fit" << std::endl;
        } else if (height != res.Size()) {
            std::cerr << "Base_Matrix::operator*(Vector): prod vector not ok" << std::endl;
        } else {
            const double* mp = &Get(1, 1);

            for (size_t i = 1; i <= height; i++) {
                double sum = b(i - 1);
                const double* xp = &x(0);

                for (size_t j = 1; j <= width; ++j, ++mp, ++xp) {
                    sum -= *mp * *xp;
                }
                res(i - 1) = sum;
            }
        }
    }

    void DenseMatrix::Solve(const Vector& v, Vector& sol) const
    {
        DenseMatrix temp(*this);
        temp.SolveDestroy(v, sol);
    }

    void DenseMatrix::SolveDestroy(const Vector& v, Vector& sol)
    {
        double q;

        if (width != height) {
            std::cerr << "SolveDestroy: Matrix not square";
            return;
        }
        if (width != v.Size()) {
            std::cerr << "SolveDestroy: Matrix and Vector don't fit";
            return;
        }

        sol = v;
        if (height != sol.Size()) {
            std::cerr << "SolveDestroy: Solution Vector not ok";
            return;
        }

        for (size_t i = 1; i <= height; i++) {
            for (size_t j = i + 1; j <= height; j++) {
                q = Get(j, i) / Get(i, i);
                if (q) {
                    const double* pik = &Get(i, i + 1);
                    double* pjk = &Elem(j, i + 1);

                    for (size_t k = i + 1; k <= height; ++k, ++pik, ++pjk) {
                        *pjk -= q * *pik;
                    }
                    sol(j - 1) -= q * sol(i - 1);
                }
            }
        }

        for (size_t i = height; i >= 1; i--) {
            q = sol(i - 1);
            for (size_t j = i + 1; j <= height; j++) {
                q -= Get(i, j) * sol(j - 1);
            }
            sol(i - 1) = q / Get(i, i);
        }
    }

    std::ostream& operator<<(std::ostream& ost, const DenseMatrix& m)
    {
        for (size_t i = 0; i < m.Height(); i++) {
            for (size_t j = 0; j < m.Width(); j++) {
                ost << m.Get(i + 1, j + 1) << " ";
            }
            ost << std::endl;
        }
        return ost;
    }

}  // namespace meshit
