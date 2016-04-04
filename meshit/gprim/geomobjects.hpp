#ifndef FILE_OBJECTS_HPP
#define FILE_OBJECTS_HPP

/* *************************************************************************/
/* File:   geomobjects.hpp                                                 */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

#include <algorithm>
#include <cstdio>
#include <cmath>

namespace meshit
{
    template<int D>
    class Vec;

    template<int H, int W = H>
    class Mat
    {
     protected:
        double x[H * W];

     public:
        Mat() { }

        Mat(const Mat& b)
        {
            for (size_t i = 0; i < H * W; i++) x[i] = b.x[i];
        }

        Mat& operator=(double s)
        {
            for (size_t i = 0; i < H * W; i++) x[i] = s;
            return *this;
        }

        Mat& operator=(const Mat& b)
        {
            for (size_t i = 0; i < H * W; i++) x[i] = b.x[i];
            return *this;
        }

        double& operator()(size_t i, size_t j)
        {
            return x[i * W + j];
        }

        const double& operator()(size_t i, size_t j) const
        {
            return x[i * W + j];
        }

        double& operator()(size_t i)
        {
            return x[i];
        }

        const double& operator()(size_t i) const
        {
            return x[i];
        }

        Vec<H> Col(size_t i) const
        {
            Vec<H> hv;
            for (size_t j = 0; j < H; j++) {
                hv(j) = x[j * W + i];
            }
            return hv;
        }

        Vec<W> Row(size_t i) const
        {
            Vec<W> hv;
            for (size_t j = 0; j < W; j++) {
                hv(j) = x[i * W + j];
            }
            return hv;
        }

        void Solve(const Vec<H>& rhs, Vec<W>& sol) const
        {
            Mat<W, H> inv;
            CalcInverse(*this, inv);
            sol = inv * rhs;
        }
    };



}  // namespace meshit

#endif
