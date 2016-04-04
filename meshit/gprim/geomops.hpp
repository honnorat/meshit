#ifndef FILE_GEOMOPS_HPP
#define FILE_GEOMOPS_HPP

/* *************************************************************************/
/* File:   geomops.hpp                                                     */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

#include <iostream>
#include "geomobjects.hpp"

namespace meshit
{
    // Matrix - Matrix operations

    inline Mat<2, 2> operator*(const Mat<2, 2>& a, const Mat<2, 2>& b)
    {
        Mat<2, 2> m;
        for (size_t i = 0; i < 2; i++) {
            for (size_t j = 0; j < 2; j++) {
                double sum = 0;
                for (size_t k = 0; k < 2; k++) {
                    sum += a(i, k) * b(k, j);
                }
                m(i, j) = sum;
            }
        }
        return m;
    }

    inline Mat<2, 2> operator*(const Mat<2, 3>& a, const Mat<3, 2>& b)
    {
        Mat<2, 2> m;
        for (size_t i = 0; i < 2; i++) {
            for (size_t j = 0; j < 2; j++) {
                double sum = 0;
                for (size_t k = 0; k < 3; k++) {
                    sum += a(i, k) * b(k, j);
                }
                m(i, j) = sum;
            }
        }
        return m;
    }

    inline Mat<3, 2> operator*(const Mat<3, 2>& a, const Mat<2, 2>& b)
    {
        Mat<3, 2> m;
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 2; j++) {
                double sum = 0;
                for (size_t k = 0; k < 2; k++) {
                    sum += a(i, k) * b(k, j);
                }
                m(i, j) = sum;
            }
        }
        return m;
    }

    inline Mat<2, 3> operator*(const Mat<2, 2>& a, const Mat<2, 3>& b)
    {
        Mat<2, 3> m;
        for (size_t i = 0; i < 2; i++) {
            m(i, 0) = a(i, 0) * b(0, 0) + a(i, 1) * b(1, 0);
            m(i, 1) = a(i, 0) * b(0, 1) + a(i, 1) * b(1, 1);
            m(i, 2) = a(i, 0) * b(0, 2) + a(i, 2) * b(1, 2);
        }
        return m;
    }

    template<int H, int W>
    inline Mat <W, H> Trans(const Mat <H, W>& m)
    {
        Mat<W, H> res;
        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++)
                res(j, i) = m(i, j);
        return res;
    }

    template<int H, int W>
    inline std::ostream& operator<<(std::ostream& ost, const Mat <H, W>& m)
    {
        ost << "(";
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++)
                ost << m(i, j) << "   ";
            ost << std::endl;
        }
        return ost;
    }

}  // namespace meshit

#endif
