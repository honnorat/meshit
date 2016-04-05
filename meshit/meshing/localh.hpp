#ifndef LOCALH_HPP
#define LOCALH_HPP

/**************************************************************************/
/* File:   localh.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   29. Jan. 97                                                    */
/**************************************************************************/

#include <iostream>
#include "../general/optmem.hpp"
#include "../gprim/geom3d.hpp"

namespace meshit
{
    class GradingBox
    {
     public:
        GradingBox(const double* ax1, const double* ax2);
        void DeleteChilds();

        void SetBox(const double* ax1, const double* ax2);

        friend class LocalH;

        static BlockAllocator ball;
        void* operator new(size_t);
        void operator delete(void*);

     protected:
        double xmid[2];
        double h2;  // half edgelength

        GradingBox* childs[4];

        double hopt;
    };

    /**
       Control of 3D mesh grading
     */
    class LocalH
    {
        GradingBox* root;
        double grading;

     public:
        LocalH()
            : root{nullptr} { }

        ~LocalH();

        void Init(const Point2d& pmin, const Point2d& pmax, double grading);
        void CleanRoot();

        void SetH(const Point2d& x, double h);
        double GetH(const Point2d& x) const;
        double GetMinH(const Point2d& pmin, const Point2d& pmax) const;

     private:
        double GetMinHRec(const Point2d& pmin, const Point2d& pmax, const GradingBox* box) const;
    };
}  // namespace meshit

#endif
