#ifndef FILE_ADFRONT2
#define FILE_ADFRONT2

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
namespace meshit {

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
            frontnr = INT_MAX - 10; // attention: overflow on calculating  INT_MAX + 1
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
        /// geometry specific data
        PointGeomInfo geominfo[2];

     public:
        FrontLine()
        {
            lineclass = 1;
        }

        FrontLine(const INDEX_2& al)
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

        void SetGeomInfo(const PointGeomInfo& gi1, const PointGeomInfo& gi2)
        {
            geominfo[0] = gi1;
            geominfo[1] = gi2;
        }

        const PointGeomInfo& GetGeomInfo(int endp) const
        {
            return geominfo[endp - 1];
        }

        friend class AdFront2;
    };

    class AdFront2
    {
        std::vector<FrontPoint2> points; // front points
        std::vector<FrontLine> lines;    // front lines

        Box3d boundingbox;
        Box3dTree linesearchtree;     // search tree for lines
        Point3dTree pointsearchtree;  // search tree for points
        Point3dTree cpointsearchtree; // search tree for cone points (not used ???)

        Array<int> delpointl;  // list of deleted front points
        Array<int> dellinel;   // list of deleted front lines

        int nfl; /// number of front lines;
        INDEX_2_HASHTABLE<int>* allflines; /// all front lines ever have been

        std::vector<int> invpindex;

        int minval;
        int starti;

     public:

        AdFront2(const Box3d& aboundingbox);
        ~AdFront2();

        bool Empty() const
        {
            return nfl == 0;
        }

        int GetNFL() const
        {
            return nfl;
        }

        int SelectBaseLine(Point<3>& p1, Point<3>& p2,
                           const PointGeomInfo*& geominfo1,
                           const PointGeomInfo*& geominfo2,
                           int& qualclass);

        int GetLocals(int baseline,
                      Array<Point3d>& locpoints,
                      Array<INDEX_2>& loclines, // local index
                      Array<int>& pindex,
                      Array<int>& lindex,
                      double xh);

        void DeleteLine(int li);

        int AddPoint(const Point3d& p, PointIndex globind);

        int AddLine(int pi1, int pi2, const PointGeomInfo& gi1, const PointGeomInfo& gi2);

        int ExistsLine(int gpi1, int gpi2);

        void IncrementClass(int li)
        {
            lines[li].IncrementClass();
        }

        const PointGeomInfo& GetLineGeomInfo(int li, int lend) const
        {
            return lines[li].GetGeomInfo(lend);
        }

        PointIndex GetGlobalIndex(int pi) const
        {
            return points[pi].GlobalIndex();
        }

        void SetStartFront();
        void PrintOpenSegments(std::ostream& ost) const;
    };

}

#endif

