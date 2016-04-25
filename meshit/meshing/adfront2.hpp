#ifndef FILE_ADFRONT2_HPP
#define FILE_ADFRONT2_HPP

/**************************************************************************/
/* File:   adfront2.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

#include <algorithm>
#include <climits>
#include <iostream>
#include <vector>

#include "../gprim/adtree.hpp"
#include "../gprim/geomobjects.hpp"
#include "mesh_types.hpp"

/**
    Advancing front class for surfaces
 */
namespace meshit {

class FrontPoint2
{
 public:
    FrontPoint2(const Point2d& point, PointIndex gi);

    ~FrontPoint2() { }

    const Point2d& P() const { return point_; }
    operator const Point2d&() const { return point_; }
    int FrontNr() const { return frontnr; }
    PointIndex GlobalIndex() const { return global_index_; }

    void AddLine() { nlinetopoint++; }

    void RemoveLine()
    {
        nlinetopoint--;
        if (nlinetopoint == 0) nlinetopoint = -1;
    }

    bool Valid() const { return nlinetopoint >= 0; }
    void DecFrontNr(int afrontnr) { frontnr = std::min(frontnr, afrontnr); }

 protected:
    Point2d point_;            // coordinates
    PointIndex global_index_;  // global node index
    int nlinetopoint;          // number of front lines connected to point
    int frontnr;               // distance to original boundary
};

class FrontLine
{
 public:
    FrontLine()
        : lineclass{1} { }

    explicit FrontLine(const INDEX_2& al)
        : l{al}, lineclass{1} { }

    const INDEX_2& L() const { return l; }
    int LineClass() const { return lineclass; }

    void IncrementClass() { lineclass++; }

    bool Valid() const { return l.I1() != -1; }

    void Invalidate()
    {
        l.I1() = -1;
        l.I2() = -1;
        lineclass = 1000;
    }

    friend class AdFront2;

 private:
    INDEX_2 l;      // Point Indizes
    int lineclass;  // quality class
};

class AdFront2
{
 public:
    explicit AdFront2(const Box2d& aboundingbox);

    ~AdFront2() { }

    bool Empty() const { return nfl == 0; }
    size_t GetNFL() const { return nfl; }

    int SelectBaseLine(Point2d& p1, Point2d& p2, int& qualclass);

    void GetLocals(int baseline, std::vector<Point2d>& locpoints,
                   std::vector<INDEX_2>& loclines,  // local index
                   std::vector<int>& pindex, std::vector<int>& lindex, double xh);

    void DeleteLine(int li);

    int AddPoint(const Point2d& p, PointIndex globind);
    int AddLine(PointIndex pi1, PointIndex pi2);

    bool LineExists(PointIndex gpi1, PointIndex gpi2);

    void IncrementClass(int li) { lines[li].IncrementClass(); }

    PointIndex GetGlobalIndex(PointIndex pi) const { return points[pi].GlobalIndex(); }

    void SetStartFront();
    void PrintOpenSegments(std::ostream& ost) const;

 protected:
    std::vector<FrontPoint2> points;  // front points
    std::vector<FrontLine> lines;     // front lines

    Box2d boundingbox;
    Box3dTree linesearchtree;      // search tree for lines
    Point3dTree pointsearchtree;   // search tree for points
    Point3dTree cpointsearchtree;  // search tree for cone points (not used ???)

    std::vector<PointIndex> delpointl;  // list of deleted front points
    std::vector<int> dellinel;          // list of deleted front lines

    std::vector<size_t> nearlines;
    std::vector<size_t> nearpoints;

    size_t nfl;                        // number of front lines;
    INDEX_2_map<uint8_t> all_flines_;  // all front lines ever have been

    std::vector<int> invpindex;

    int minval;
    size_t starti;
};

}  // namespace meshit

#endif
