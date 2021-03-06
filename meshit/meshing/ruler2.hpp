#ifndef MESHIT_RULER_HPP
#define MESHIT_RULER_HPP
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

#include <iostream>
#include <vector>

#include "../gprim/geom2d.hpp"
#include "../linalg/densemat.hpp"
#include "mesh_types.hpp"

namespace meshit {

class netrule
{
 public:
    netrule();
    ~netrule();

    struct threefloat
    {
        double f1, f2, f3;
    };

    struct threeint
    {
        int i1, i2, i3;
    };

    size_t GetNP() const { return points.size(); }
    size_t GetNL() const { return lines.size(); }
    size_t GetNE() const { return elements.size(); }
    size_t GetNOldP() const { return noldp; }
    size_t GetNOldL() const { return noldl; }
    size_t GetNDelL() const { return dellines.size(); }
    size_t GetNOrientations() const { return orientations.size(); }

    int GetQuality() const { return quality; }
    const char* Name() const { return name; }

    uint32_t GetLNearness(size_t li) const { return lnearness[li]; }
    const Point2d& GetPoint(size_t i) const { return points[i]; }
    const IndexPair& GetLine(size_t i) const { return lines[i]; }
    const Element2d& GetElement(size_t i) const { return elements[i]; }
    const threeint& GetOrientation(size_t i) const { return orientations[i]; }

    int GetDelLine(size_t i) const { return dellines[i]; }
    const DenseMatrix& GetOldUToNewU() const { return oldutonewu; }

    double CalcPointDist(size_t pi, const Point2d& p) const
    {
        double dx = p.X() - points[pi].X();
        double dy = p.Y() - points[pi].Y();
        const threefloat* tfp = &tolerances[pi];
        return tfp->f1 * dx * dx + tfp->f2 * dx * dy + tfp->f3 * dy * dy;
    }

    double CalcLineError(size_t li, const Vec2d& v) const;

    void SetFreeZoneTransformation(const Vector& u, int tolclass);

    bool IsInFreeZone(const Point2d& p) const;
    bool IsLineInFreeZone(const Point2d& p1, const Point2d& p2) const;

    int ConvexFreeZone() const;

    GenericIndex GetPointNr1(size_t ln) const { return lines[ln].I1(); }
    GenericIndex GetPointNr2(size_t ln) const { return lines[ln].I2(); }

    void LoadRule(std::istream& ist);

 private:
    int quality;
    char* name;
    std::vector<Point2d> points;
    std::vector<IndexPair> lines;
    std::vector<Point2d> freezone, freezonelimit;
    std::vector<std::vector<Point2d>*> freezone_i;
    std::vector<Point2d> transfreezone;

    std::vector<int> dellines;
    std::vector<Element2d> elements;
    std::vector<threefloat> tolerances, linetolerances;
    std::vector<threeint> orientations;
    DenseMatrix oldutonewu, oldutofreearea, oldutofreearealimit;
    std::vector<DenseMatrix*> oldutofreearea_i;
    MatrixFixWidth<3> freesetinequ;

    std::vector<Vec2d> linevecs;

    size_t noldp, noldl;
    double fzminx, fzmaxx, fzminy, fzmaxy;

    /// topological distance of line to base element
    std::vector<uint32_t> lnearness;
};

}  // namespace meshit

#endif
