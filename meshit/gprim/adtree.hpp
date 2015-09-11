#ifndef FILE_ADTREE
#define FILE_ADTREE

/* *************************************************************************/
/* File:   adtree.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   16. Feb. 98                                                     */
/* Redesigned by Wolfram Muehlhuber, May 1998                              */
/* *************************************************************************/

#include <iostream>

#include <meshit.hpp>

#include "../general/optmem.hpp"
#include "geomobjects.hpp"

namespace meshit {

    /**
      Alternating Digital Tree
     */

    class ADTreeNode
    {
      public:
        ADTreeNode *left, *right, *father;
        int dim;
        double sep;
        double *data;
        double *boxmin;
        double *boxmax;
        int pi;
        int nchilds;

        ADTreeNode(int adim);
        ~ADTreeNode();

        friend class ADTree;
    };

    class ADTreeCriterion
    {
      public:

        ADTreeCriterion() { }
        virtual int Eval(const ADTreeNode * node) const = 0;
    };

    class ADTree
    {
        int dim;
        ADTreeNode * root;
        double *cmin, *cmax;
        Array<ADTreeNode*> ela;
        const ADTreeCriterion * criterion;

        Array<ADTreeNode*> stack;
        Array<int> stackdir;
        int stackindex;

      public:
        ADTree(int adim, const double * acmin, const double * acmax);
        ~ADTree() {}

        void Insert(const double * p, int pi);

        void SetCriterion(ADTreeCriterion & acriterion);
        void Reset();
        int Next();
        void GetMatch(Array<int> & matches);

        void DeleteElement(int pi);

        void Print(std::ostream & ost) const
        {
            PrintRec(ost, root);
        }

        void PrintRec(std::ostream & ost, const ADTreeNode * node) const;
    };

    class ADTreeNode3
    {
      public:
        ADTreeNode3 *left, *right, *father;
        double sep;
        double data[3];
        int pi;
        int nchilds;

        ADTreeNode3();
        void DeleteChilds();
        friend class ADTree3;

        static BlockAllocator ball;
        void * operator new(size_t);
        void operator delete (void *);
    };

    class ADTree3
    {
        ADTreeNode3 * root;
        double cmin[3], cmax[3];
        Array<ADTreeNode3*> ela;

      public:
        ADTree3(const double * acmin,
                const double * acmax);
        ~ADTree3();

        void Insert(const double * p, int pi);
        void GetIntersecting(const double * bmin, const double * bmax,
                Array<int> & pis) const;

        void DeleteElement(int pi);

        void Print(std::ostream & ost) const
        {
            PrintRec(ost, root);
        }

        void PrintRec(std::ostream & ost, const ADTreeNode3 * node) const;
    };

    class ADTreeNode6
    {
      public:
        ADTreeNode6 *left, *right, *father;
        double sep;
        double data[6];
        int pi;
        int nchilds;

        ADTreeNode6();
        void DeleteChilds();
        friend class ADTree6;

        static BlockAllocator ball;
        void * operator new(size_t);
        void operator delete (void *);
    };

    class ADTree6
    {
        ADTreeNode6 * root;
        double cmin[6], cmax[6];
        Array<ADTreeNode6*> ela;

      public:
        ADTree6(const double * acmin,
                const double * acmax);
        ~ADTree6();

        void Insert(const double * p, int pi);
        void GetIntersecting(const double * bmin, const double * bmax,
                Array<int> & pis) const;

        void DeleteElement(int pi);

        void Print(std::ostream & ost) const
        {
            PrintRec(ost, root);
        }

        int Depth() const
        {
            return DepthRec(root);
        }

        int Elements() const
        {
            return ElementsRec(root);
        }

        void PrintRec(std::ostream & ost, const ADTreeNode6 * node) const;
        int DepthRec(const ADTreeNode6 * node) const;
        int ElementsRec(const ADTreeNode6 * node) const;

        void PrintMemInfo(std::ostream & ost) const;
    };

    class Point3dTree
    {
        ADTree3 * tree;

      public:
        Point3dTree(const Point<3> & pmin, const Point<3> & pmax);
        ~Point3dTree();
        void Insert(const Point<3> & p, int pi);

        void DeleteElement(int pi)
        {
            tree->DeleteElement(pi);
        }
        void GetIntersecting(const Point<3> & pmin, const Point<3> & pmax,
                Array<int> & pis) const;

        const ADTree3 & Tree() const
        {
            return *tree;
        };
    };

    class Box3dTree
    {
        ADTree6 * tree;
        Point<3> boxpmin, boxpmax;
      public:
        Box3dTree(const Box<3> & abox);
        Box3dTree(const Point<3> & apmin, const Point<3> & apmax);
        ~Box3dTree();
        void Insert(const Point<3> & bmin, const Point<3> & bmax, int pi);

        void Insert(const Box<3> & box, int pi)
        {
            Insert(box.PMin(), box.PMax(), pi);
        }

        void DeleteElement(int pi)
        {
            tree->DeleteElement(pi);
        }
        void GetIntersecting(const Point<3> & pmin, const Point<3> & pmax,
                Array<int> & pis) const;

        const ADTree6 & Tree() const
        {
            return *tree;
        };
    };

}

#endif
