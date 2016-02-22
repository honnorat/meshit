#ifndef LOCALH
#define LOCALH

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
        double xmid[3];
        double h2; // half edgelength

        GradingBox * childs[8];
        GradingBox * father;

        double hopt;

      public:

        struct
        {
            unsigned int cutboundary : 1;
            unsigned int isinner : 1;
            unsigned int oldcell : 1;
            unsigned int pinner : 1;
        } flags;

        GradingBox(const double * ax1, const double * ax2);
        void DeleteChilds();

        friend class LocalH;

        static BlockAllocator ball;
        void * operator new(size_t);
        void operator delete (void *);
    };

    /**
       Control of 3D mesh grading
     */
    class LocalH
    {
        GradingBox * root;
        double grading;
        Array<GradingBox*> boxes;

     public:

        LocalH(const Point3d & pmin, const Point3d & pmax, double grading);
        LocalH(const Box<3> & box, double grading);
        ~LocalH();

        void SetH(const Point3d & x, double h);
        double GetH(const Point3d & x) const;
        double GetMinH(const Point3d & pmin, const Point3d & pmax) const;


      private:

        double GetMinHRec(const Point3d & pmin, const Point3d & pmax, const GradingBox * box) const;

        friend std::ostream & operator<<(std::ostream & ost, const LocalH & loch);
    };
}

#endif
