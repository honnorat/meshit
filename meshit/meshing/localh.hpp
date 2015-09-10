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

        Point<3> PMid() const
        {
            return Point<3> (xmid[0], xmid[1], xmid[2]);
        }

        double H2() const
        {
            return h2;
        }

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
        Box3d boundingbox;

      public:

        LocalH(const Point3d & pmin, const Point3d & pmax, double grading);
        LocalH(const Box<3> & box, double grading);
        ~LocalH();
        void Delete();

        void SetGrading(double agrading)
        {
            grading = agrading;
        }

        void SetH(const Point3d & x, double h);
        double GetH(const Point3d & x) const;
        double GetMinH(const Point3d & pmin, const Point3d & pmax) const;

        void CutBoundary(const Box<3> & box)
        {
            CutBoundaryRec(box.PMin(), box.PMax(), root);
        }

        void FindInnerBoxes(class AdFront3 * adfront, int (*testinner)(const Point3d & p1));

        void FindInnerBoxes(class AdFront2 * adfront, int (*testinner)(const Point<2> & p1));

        void ClearFlags()
        {
            ClearFlagsRec(root);
        }

        void GetInnerPoints(Array<Point<3> > & points);
        void GetOuterPoints(Array<Point<3> > & points);

        void Convexify();

        int GetNBoxes()
        {
            return boxes.size();
        }

        const Box3d & GetBoundingBox() const
        {
            return boundingbox;
        }

      private:

        double GetMinHRec(const Point3d & pmin, const Point3d & pmax, const GradingBox * box) const;

        void CutBoundaryRec(const Point3d & pmin, const Point3d & pmax, GradingBox * box);

        void FindInnerBoxesRec(int (*inner)(const Point3d & p), GradingBox * box);

        void FindInnerBoxesRec2(
                GradingBox * box,
                class AdFront3 * adfront,
                Array<Box3d> & faceboxes,
                Array<int> & finds, int nfinbox);

        void SetInnerBoxesRec(GradingBox * box);

        void ClearFlagsRec(GradingBox * box);

        friend std::ostream & operator<<(std::ostream & ost, const LocalH & loch);
    };
}

#endif
