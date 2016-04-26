#ifndef MESHIT_ADTREE_HPP
#define MESHIT_ADTREE_HPP
/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 1998 Wolfram Muehlhuber
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
#include "geom3d.hpp"
#include "geomobjects.hpp"
#include "../meshing/mesh_types.hpp"

namespace meshit {
/**
  Alternating Digital Tree
 */
class ADTreeNode3
{
 public:
    ADTreeNode3()
        : left_{nullptr}, right_{nullptr}, pi_{CONST<GenericIndex>::undefined} { }

    void DeleteChilds();
    void SetData(const Point2d& point, GenericIndex pi);

    void* operator new(size_t);
    void operator delete(void* p);

    friend class ADTree3;

 protected:
    ADTreeNode3* left_;
    ADTreeNode3* right_;
    double sep_;
    double data_[2];
    GenericIndex pi_;

    static BlockAllocator ball_;
};

class ADTree3
{
 public:
    ADTree3(const Point2d& cmin, const Point2d& cmax);
    ~ADTree3();

    void Insert(const Point2d& p, GenericIndex pi);
    void GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<GenericIndex>& pis) const;
    void DeleteElement(GenericIndex pi);

 protected:
    ADTreeNode3* root_;
    Point2d cmin_;
    Point2d cmax_;
    std::vector<ADTreeNode3*> nodes_;
};

class ADTreeNode6
{
 public:
    ADTreeNode6()
        : left_{nullptr}, right_{nullptr}, pi_{CONST<GenericIndex>::undefined} { }

    void DeleteChilds();
    void SetData(const Point2d& pmin, const Point2d& pmax, GenericIndex pi);

    friend class ADTree6;

    void* operator new(size_t);
    void operator delete(void* p);

 protected:
    ADTreeNode6* left_;
    ADTreeNode6* right_;
    double sep_;
    Point2d pmin_, pmax_;
    GenericIndex pi_;

    static BlockAllocator ball_;
};

class ADTree6
{
 public:
    ADTree6(const Point2d& pmin, const Point2d& pmax);
    ~ADTree6();

    void Insert(const Point2d& bmin, const Point2d& bmax, GenericIndex pi);
    void GetIntersecting(const Point2d& bmin, const Point2d& bmax,
                         const Point2d& pmin, const Point2d& pmax,
                         std::vector<GenericIndex>& pis) const;
    void DeleteElement(GenericIndex pi);

 protected:
    ADTreeNode6* root_;
    Point2d pmin_;
    Point2d pmax_;
    std::vector<ADTreeNode6*> nodes_;
};

class Point3dTree
{
 public:
    Point3dTree(const Point2d& pmin, const Point2d& pmax);
    ~Point3dTree();

    void Insert(const Point2d& p, PointIndex pi);
    void DeleteElement(PointIndex pi) { tree_->DeleteElement(pi); }
    void GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<PointIndex>& pis) const;

 protected:
    ADTree3* tree_;
};

class Box3dTree
{
 public:
    Box3dTree(const Point2d& apmin, const Point2d& apmax);
    ~Box3dTree();

    void Insert(const Point2d& bmin, const Point2d& bmax, GenericIndex pi);
    void DeleteElement(GenericIndex pi) { tree_->DeleteElement(pi); }
    void GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<GenericIndex>& pis) const;

 protected:
    ADTree6* tree_;
    Point2d box_pmin_, box_pmax_;
};

}  // namespace meshit

#endif
