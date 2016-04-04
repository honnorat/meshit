#include "geomfuncs.hpp"

namespace meshit
{
    inline void CalcInverse(const Mat<2, 2>& m, Mat<2, 2>& inv)
    {
        double det = m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
        if (det == 0) {
            inv = 0;
            return;
        }

        double idet = 1.0 / det;
        inv(0, 0) = idet * m(1, 1);
        inv(0, 1) = -idet * m(0, 1);
        inv(1, 0) = -idet * m(1, 0);
        inv(1, 1) = idet * m(0, 0);
    }

    void CalcInverse(const Mat<3, 3>& m, Mat<3, 3>& inv)
    {
        double det = m(0, 0) * m(1, 1) * m(2, 2)
                     + m(1, 0) * m(2, 1) * m(0, 2)
                     + m(2, 0) * m(0, 1) * m(1, 2)
                     - m(0, 0) * m(2, 1) * m(1, 2)
                     - m(1, 0) * m(0, 1) * m(2, 2)
                     - m(2, 0) * m(1, 1) * m(0, 2);;

        if (det == 0) {
            inv = 0;
            return;
        }

        double idet = 1.0 / det;
        inv(0, 0) = idet * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1));
        inv(1, 0) = -idet * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0));
        inv(2, 0) = idet * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));

        inv(0, 1) = -idet * (m(0, 1) * m(2, 2) - m(0, 2) * m(2, 1));
        inv(1, 1) = idet * (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0));
        inv(2, 1) = -idet * (m(0, 0) * m(2, 1) - m(0, 1) * m(2, 0));

        inv(0, 2) = idet * (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1));
        inv(1, 2) = -idet * (m(0, 0) * m(1, 2) - m(0, 2) * m(1, 0));
        inv(2, 2) = idet * (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0));
    }

    inline void CalcInverse(const Mat<2, 3>& m, Mat<3, 2>& inv)
    {
        Mat<2, 2> a = m * Trans(m);
        Mat<2, 2> ainv;
        CalcInverse(a, ainv);
        inv = Trans(m) * ainv;
    }

    inline void CalcInverse(const Mat<3, 2>& m, Mat<2, 3>& inv)
    {
        Mat<2, 2> a = Trans(m) * m;
        Mat<2, 2> ainv;
        CalcInverse(a, ainv);
        inv = ainv * Trans(m);
    }

    inline double Det(const Mat<2, 2>& m)
    {
        return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
    }

    double Det(const Mat<3, 3>& m)
    {
        return m(0, 0) * m(1, 1) * m(2, 2)
               + m(1, 0) * m(2, 1) * m(0, 2)
               + m(2, 0) * m(0, 1) * m(1, 2)
               - m(0, 0) * m(2, 1) * m(1, 2)
               - m(1, 0) * m(0, 1) * m(2, 2)
               - m(2, 0) * m(1, 1) * m(0, 2);
    }

}  // namespace meshit
