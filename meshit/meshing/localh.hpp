#ifndef MESHIT_LOCALH_HPP
#define MESHIT_LOCALH_HPP
/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2015-2016 Marc Honnorat <marc.honnorat@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library in the file LICENSE.LGPL; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

#include <iostream>

#include "../general/block_allocator.hpp"
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
