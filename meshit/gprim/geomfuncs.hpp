#ifndef FILE_GEOMFUNCS_HPP
#define FILE_GEOMFUNCS_HPP

/* *************************************************************************/
/* File:   geomfuncs.hpp                                                   */
/* Author: Joachim Schoeberl                                               */
/* Date:   20. Jul. 02                                                     */
/* *************************************************************************/

#include "geomobjects.hpp"
#include "geomops.hpp"

namespace meshit
{
    void CalcInverse(const Mat<2, 2>& m, Mat<2, 2>& inv);
    void CalcInverse(const Mat<3, 3>& m, Mat<3, 3>& inv);
    void CalcInverse(const Mat<2, 3>& m, Mat<3, 2>& inv);
    void CalcInverse(const Mat<3, 2>& m, Mat<2, 3>& inv);

    double Det(const Mat<2, 2>& m);
    double Det(const Mat<3, 3>& m);

}  // namespace meshit

#endif
