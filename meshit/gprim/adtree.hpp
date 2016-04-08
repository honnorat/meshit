#ifndef FILE_ADTREE_HPP
#define FILE_ADTREE_HPP

/* *************************************************************************/
/* File:   adtree.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   16. Feb. 98                                                     */
/* Redesigned by Wolfram Muehlhuber, May 1998                              */
/* *************************************************************************/

#include <iostream>

#include "../meshit.hpp"
#include "../general/optmem.hpp"
#include "geom3d.hpp"
#include "geomobjects.hpp"

namespace meshit
{
    /**
      Alternating Digital Tree
     */
    class ADTreeNode3
    {
     public:
        ADTreeNode3* left, * right, * father;
        double sep;
        double data[2];
        int pi;
        int nchilds;

        ADTreeNode3();
        void DeleteChilds();

        friend class ADTree3;

        static BlockAllocator ball;
        void* operator new(size_t);
        void operator delete(void* p);
    };

    class ADTree3
    {
     protected:
        ADTreeNode3* root;
        double cmin[2], cmax[2];
        std::vector<ADTreeNode3*> ela;

     public:
        ADTree3(const Point2d& acmin, const Point2d& acmax);
        ~ADTree3();

        void Insert(const Point2d& p, int pi);

        void GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const;

        void DeleteElement(int pi);

        void PrintRec(std::ostream& ost, const ADTreeNode3* node) const;
    };

    class ADTreeNode6
    {
     public:
        ADTreeNode6* left, * right, * father;
        double sep_;
        Point2d pmin_, pmax_;
        int pi_;
        int nchilds;

        ADTreeNode6();
        void DeleteChilds();
        void SetData(const Point2d& pmin, const Point2d& pmax, int pi);

        friend class ADTree6;

        static BlockAllocator ball;
        void* operator new(size_t);
        void operator delete(void* p);
    };

    class ADTree6
    {
        ADTreeNode6* root;
        Point2d pmin_;
        Point2d pmax_;
        std::vector<ADTreeNode6*> ela;

     public:
        ADTree6(const Point2d& pmin, const Point2d& pmax);
        ~ADTree6();

        void Insert(const Point2d& bmin, const Point2d& bmax, int pi);
        void GetIntersecting(const Point2d& bmin, const Point2d& bmax,
                             const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const;
        void DeleteElement(int pi);
    };

    class Point3dTree
    {
        ADTree3* tree;

     public:
        Point3dTree(const Point2d& pmin, const Point2d& pmax);
        ~Point3dTree();

        void Insert(const Point2d& p, int pi);

        void DeleteElement(int pi)
        {
            tree->DeleteElement(pi);
        }

        void GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const;
    };

    class Box3dTree
    {
     protected:
        ADTree6* tree;
        Point2d boxpmin, boxpmax;

     public:
        Box3dTree(const Point2d& apmin, const Point2d& apmax);
        ~Box3dTree();

        void Insert(const Point2d& bmin, const Point2d& bmax, int pi);

        void DeleteElement(int pi)
        {
            tree->DeleteElement(pi);
        }

        void GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const;
    };

}  // namespace meshit

#endif
