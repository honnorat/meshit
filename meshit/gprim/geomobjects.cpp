#include "geomobjects.hpp"

namespace meshit {

double Mat3x3::Det() const
{
    return x[0] * (x[4] * x[8] - x[7] * x[5]) +  // m(0,0) * ( m(1,1)*m(2,2) - m(2,1)*m(1,2) )
           x[3] * (x[7] * x[2] - x[1] * x[8]) +  // m(1,0) * ( m(2,1)*m(0,2) - m(0,1)*m(2,2) )
           x[6] * (x[1] * x[5] - x[4] * x[2]);   // m(2,0) * ( m(0,1)*m(1,2) - m(1,1)*m(0,2) )
}

void Mat3x3::CalcInverse(Mat3x3& inv) const
{
    double det = Det();

    if (det == 0) {
        inv = 0;
        return;
    }

    double idet = 1.0 / det;

    inv(0, 0) = +idet * (x[4] * x[8] - x[5] * x[7]);  // (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1));
    inv(1, 0) = -idet * (x[3] * x[8] - x[5] * x[6]);  // (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0));
    inv(2, 0) = +idet * (x[3] * x[7] - x[4] * x[6]);  // (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));

    inv(0, 1) = -idet * (x[1] * x[8] - x[2] * x[7]);  // (m(0, 1) * m(2, 2) - m(0, 2) * m(2, 1));
    inv(1, 1) = +idet * (x[0] * x[8] - x[2] * x[6]);  // (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0));
    inv(2, 1) = -idet * (x[0] * x[7] - x[1] * x[6]);  // (m(0, 0) * m(2, 1) - m(0, 1) * m(2, 0));

    inv(0, 2) = +idet * (x[1] * x[5] - x[2] * x[4]);  // (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1));
    inv(1, 2) = -idet * (x[0] * x[5] - x[2] * x[3]);  // (m(0, 0) * m(1, 2) - m(0, 2) * m(1, 0));
    inv(2, 2) = +idet * (x[0] * x[4] - x[1] * x[3]);  // (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0));
}

}  // namespace meshit
