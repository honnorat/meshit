#ifndef MESHIT_ADFRONT2_HPP
#define MESHIT_ADFRONT2_HPP
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

typedef GenericIndex FrontLineIndex;

class FrontPoint2
{
 public:
    FrontPoint2(const Point2d& point, PointIndex gi);

    ~FrontPoint2() { }

    const Point2d& P() const { return point_; }
    operator const Point2d&() const { return point_; }
    uint32_t FrontNr() const { return frontnr; }
    PointIndex GlobalIndex() const { return global_index_; }

    void AddLine() { nlinetopoint++; }

    void RemoveLine()
    {
        nlinetopoint--;
        if (nlinetopoint == 0) nlinetopoint = -1;
    }

    bool Valid() const { return nlinetopoint >= 0; }
    void DecFrontNr(uint32_t afrontnr) { frontnr = std::min(frontnr, afrontnr); }

 protected:
    Point2d point_;            // coordinates
    PointIndex global_index_;  // global node index
    int nlinetopoint;          // number of front lines connected to point
    uint32_t frontnr;               // distance to original boundary
};

class FrontLine
{
 public:
    FrontLine()
        : line_class_{1} { }

    FrontLine(PointIndex pi1, PointIndex pi2)
        : indices_{static_cast<GenericIndex>(pi1), static_cast<GenericIndex>(pi2)},
          line_class_{1} { }

    const IndexPair& L() const { return indices_; }
    uint32_t LineClass() const { return line_class_; }

    void IncrementClass() { line_class_++; }

    bool Valid() const { return indices_.I1() != CONST<PointIndex>::undefined; }

    void Invalidate()
    {
        indices_.I1() = CONST<PointIndex>::undefined;
        indices_.I2() = CONST<PointIndex>::undefined;
        line_class_ = 1000;
    }

    friend class AdFront2;

 private:
    IndexPair indices_;    // Point Indizes
    uint32_t line_class_;  // quality class
};

class AdFront2
{
 public:
    explicit AdFront2(const Box2d& aboundingbox);

    ~AdFront2() { }

    bool Empty() const { return nfl == 0; }
    size_t GetNFL() const { return nfl; }

    FrontLineIndex SelectBaseLine(Point2d& p1, Point2d& p2, uint32_t& qualclass);

    void GetLocals(FrontLineIndex baseline,
                   std::vector<Point2d>& locpoints,
                   std::vector<IndexPair>& loclines,  // local index
                   std::vector<PointIndex>& pindex,
                   std::vector<FrontLineIndex>& lindex, double xh);

    void DeleteLine(FrontLineIndex li);

    PointIndex AddPoint(const Point2d& p, PointIndex globind);
    int AddLine(PointIndex pi1, PointIndex pi2);

    bool LineExists(PointIndex gpi1, PointIndex gpi2);

    void IncrementClass(FrontLineIndex li) { lines[li].IncrementClass(); }

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

    std::vector<PointIndex> delpointl;     // list of deleted front points
    std::vector<FrontLineIndex> dellinel;  // list of deleted front lines

    std::vector<PointIndex> nearpoints;
    std::vector<FrontLineIndex> nearlines;

    size_t nfl;                          // number of front lines;
    IndexPair_map<uint8_t> all_flines_;  // all front lines ever have been

    std::vector<PointIndex> inv_pindex;

    uint32_t minval;
    FrontLineIndex starti;
};

}  // namespace meshit

#endif
