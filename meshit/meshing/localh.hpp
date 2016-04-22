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

namespace meshit {

class GradingBox
{
 public:
    GradingBox(const Point2d& px1, const Point2d& px2);

    void SetBox(const Point2d& px1, const Point2d& px2);

    void DeleteChilds();

    // Allocation :
    static BlockAllocator ball;
    void* operator new(size_t);
    void operator delete(void* p);

    friend class LocalH;

 protected:
    Point2d xmid;
    double h2;  // half edgelength
    GradingBox* childs[4];
    double hopt;
};

/**
   Control of 3D mesh grading
 */
class LocalH
{
 public:
    LocalH()
        : root{nullptr} { }

    ~LocalH();

    void Init(const Box2d& bbox, double grading);
    void CleanRoot();

    void SetH(const Point2d& x, double h);
    double GetH(const Point2d& x) const;
    double GetMinH(const Point2d& pmin, const Point2d& pmax) const;

 private:
    double GetMinHRec(const Point2d& pmin, const Point2d& pmax, const GradingBox* box) const;

 protected:
    GradingBox* root;
    double grading;
};

}  // namespace meshit

#endif
