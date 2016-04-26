#ifndef MESHIT_MESH_GENERATOR_HPP
#define MESHIT_MESH_GENERATOR_HPP
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

#include "adfront2.hpp"
#include "mesh_class.hpp"
#include "ruler2.hpp"

namespace meshit {

/*
    The basis class for 2D mesh generation.
    Has the method GenerateMesh
 */

class MeshGenerator
{
 public:
    explicit MeshGenerator(Mesh& amesh, const Box2d& aboundingbox);
    ~MeshGenerator();

    // Load rules, either from file, or compiled rules
    void LoadRules(const char* filename);
    void Reset();
    bool GenerateMesh(const MeshingParameters& mp, double gh, DomainIndex facenr);

    void AddPoint(const Point2d& p, PointIndex globind);
    void AddBoundaryElement(PointIndex i1, PointIndex i2);
    void SetMaxArea(double amaxarea);

 protected:
    void DefineTransformation(const Point2d& p1, const Point2d& p2);
    void TransformToPlain(const Point2d& locpoint, Point2d& plainpoint, double h);
    void TransformFromPlain(Point2d& plainpoint, Point2d& locpoint, double h);

    /** Applies 2D rules.
     Tests all 2D rules */
    int ApplyRules(std::vector<Point2d>& lpoints,
                   std::vector<bool>& legalpoints, size_t maxlegalpoint,
                   std::vector<IndexPair>& llines, size_t maxlegelline,
                   std::vector<Element2d>& elements,
                   std::vector<uint32_t>& dellines,
                   int tolerance, const MeshingParameters& mp);

 protected:
    Mesh& mesh_;
    const Box2d& boundingbox_;
    AdFront2* adfront;  // the current advancing front
    std::vector<netrule*> rules;  // rules for mesh generation

    // statistics
    std::vector<int> ruleused, canuse, foundmap;
    double max_area;

    Vec2d ex, ey;
    Point2d glob_p1;
};

}  // namespace meshit

#endif
