#ifndef FILE_ADFRONT2_HPP
#define FILE_ADFRONT2_HPP

/**************************************************************************/
/* File:   adfront2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include <climits>
#include <iostream>
#include <vector>

#include "../gprim/geomobjects.hpp"
#include "../gprim/adtree.hpp"
#include "meshtype.hpp"

/**
    Advancing front class for surfaces
 */
namespace meshit
{
    class FrontPoint2
    {
        /// coordinates
        Point3d p;
        /// global node index
        PointIndex globalindex;
        /// number of front lines connected to point
        int nlinetopoint;
        /// distance to original boundary
        int frontnr;

     public:
        FrontPoint2()
        {
            nlinetopoint = 0;
            frontnr = INT_MAX - 10;  // attention: overflow on calculating  INT_MAX + 1
        }

        FrontPoint2(const Point3d& ap, PointIndex agi);

        ~FrontPoint2() { }

        const Point3d& P() const
        {
            return p;
        }

        operator const Point3d&() const
        {
            return p;
        }

        PointIndex GlobalIndex() const
        {
            return globalindex;
        }

        void AddLine()
        {
            nlinetopoint++;
        }

        void RemoveLine()
        {
            nlinetopoint--;
            if (nlinetopoint == 0)
                nlinetopoint = -1;
        }

        bool Valid() const
        {
            return nlinetopoint >= 0;
        }

        void DecFrontNr(int afrontnr)
        {
            if (frontnr > afrontnr) frontnr = afrontnr;
        }

        int FrontNr() const
        {
            return frontnr;
        }
    };

    class FrontLine
    {
     private:
        /// Point Indizes
        INDEX_2 l;
        /// quality class
        int lineclass;

     public:
        FrontLine()
        {
            lineclass = 1;
        }

        explicit FrontLine(const INDEX_2& al)
        {
            l = al;
            lineclass = 1;
        }

        const INDEX_2& L() const
        {
            return l;
        }

        int LineClass() const
        {
            return lineclass;
        }

        void IncrementClass()
        {
            lineclass++;
        }

        bool Valid() const
        {
            return l.I1() != -1;
        }

        void Invalidate()
        {
            l.I1() = -1;
            l.I2() = -1;
            lineclass = 1000;
        }

        friend class AdFront2;
    };

    class AdFront2
    {
        std::vector<FrontPoint2> points;  // front points
        std::vector<FrontLine> lines;     // front lines

        Box3d boundingbox;
        Box3dTree linesearchtree;      // search tree for lines
        Point3dTree pointsearchtree;   // search tree for points
        Point3dTree cpointsearchtree;  // search tree for cone points (not used ???)

        std::vector<int> delpointl;  // list of deleted front points
        std::vector<int> dellinel;   // list of deleted front lines

        std::vector<size_t> nearlines;
        std::vector<size_t> nearpoints;


        int nfl;  // number of front lines;
        INDEX_2_map<int> allflines;  // all front lines ever have been

        std::vector<int> invpindex;

        int minval;
        size_t starti;

     public:
        explicit AdFront2(const Box3d& aboundingbox);

        ~AdFront2() { }

        bool Empty() const
        {
            return nfl == 0;
        }

        int GetNFL() const
        {
            return nfl;
        }

        int SelectBaseLine(Point3d& p1, Point3d& p2, int& qualclass);

        int GetLocals(int baseline,
                      std::vector<Point3d>& locpoints,
                      std::vector<INDEX_2>& loclines,  // local index
                      std::vector<int>& pindex,
                      std::vector<int>& lindex,
                      double xh);

        void DeleteLine(int li);

        int AddPoint(const Point3d& p, PointIndex globind);

        int AddLine(int pi1, int pi2);

        int ExistsLine(int gpi1, int gpi2);

        void IncrementClass(int li)
        {
            lines[li].IncrementClass();
        }

        PointIndex GetGlobalIndex(int pi) const
        {
            return points[pi].GlobalIndex();
        }

        void SetStartFront();
        void PrintOpenSegments(std::ostream& ost) const;
    };
}  // namespace meshit

#endif

