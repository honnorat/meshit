#ifndef FILE_ADTREE
#define FILE_ADTREE

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

namespace meshit {

    /**
      Alternating Digital Tree
     */

    class ADTreeNode3
    {
     public:
        ADTreeNode3* left, * right, * father;
        double sep;
        double data[3];
        int pi;
        int nchilds;

        ADTreeNode3();
        void DeleteChilds();

        friend class ADTree3;

        static BlockAllocator ball;
        void* operator new(size_t);
        void operator delete(void*);
    };

    class ADTree3
    {
     protected:
        ADTreeNode3* root;
        double cmin[3], cmax[3];
        std::vector<ADTreeNode3*> ela;

     public:
        ADTree3(const double* acmin, const double* acmax);
        ~ADTree3();

        void Insert(const Point<3>& p, int pi);

        void GetIntersecting(const double* bmin, const double* bmax, std::vector<size_t>& pis) const;

        void DeleteElement(int pi);

        void PrintRec(std::ostream& ost, const ADTreeNode3* node) const;
    };

    class ADTreeNode6
    {
     public:
        ADTreeNode6* left, * right, * father;
        double sep;
        double data[6];
        int pi;
        int nchilds;

        ADTreeNode6();
        void DeleteChilds();

        friend class ADTree6;

        static BlockAllocator ball;
        void* operator new(size_t);
        void operator delete(void*);
    };

    class ADTree6
    {
        ADTreeNode6* root;
        double cmin[6], cmax[6];
        std::vector<ADTreeNode6*> ela;

     public:
        ADTree6(const double* acmin,
                const double* acmax);
        ~ADTree6();

        void Insert(const Point3d& bmin, const Point3d& bmax, int pi);
        void GetIntersecting(const double* bmin, const double* bmax, std::vector<size_t>& pis) const;

        void DeleteElement(int pi);

        int DepthRec(const ADTreeNode6* node) const;
        int ElementsRec(const ADTreeNode6* node) const;

    };

    class Point3dTree
    {
        ADTree3* tree;

     public:
        Point3dTree(const Point<3>& pmin, const Point<3>& pmax);
        ~Point3dTree();
        void Insert(const Point<3>& p, int pi);

        void DeleteElement(int pi)
        {
            tree->DeleteElement(pi);
        }

        void GetIntersecting(const Point<3>& pmin, const Point<3>& pmax, std::vector<size_t>& pis) const;
    };

    class Box3dTree
    {
     protected:
        ADTree6* tree;
        Point3d boxpmin, boxpmax;

     public:
        explicit Box3dTree(const Box<3>& abox);
        Box3dTree(const Point3d& apmin, const Point3d& apmax);
        ~Box3dTree();

        void Insert(const Point3d& bmin, const Point3d& bmax, int pi);

        void Insert(const Box<3>& box, int pi)
        {
            Insert(box.PMin(), box.PMax(), pi);
        }

        void DeleteElement(int pi)
        {
            tree->DeleteElement(pi);
        }

        void GetIntersecting(const Point3d& pmin, const Point3d& pmax, std::vector<size_t>& pis) const;
    };

}  // namespace meshit

#endif
