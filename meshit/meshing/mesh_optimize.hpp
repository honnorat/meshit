#ifndef MESHIT_MESH_OPTIMIZE_HPP
#define MESHIT_MESH_OPTIMIZE_HPP
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

#include "mesh_class.hpp"

namespace meshit {

class MeshOptimize
{
 public:
    explicit MeshOptimize(Mesh& mesh)
        : mesh_{mesh} { }

    ~MeshOptimize() { }

    void EdgeSwapping(bool use_metric, double metric_weight = 0.0);
    void ImproveMesh(double metric_weight);
    void CombineImprove();

    friend class Opti2SurfaceMinFunction;

 protected:
    void EdgeSwapping(DomainIndex face_index, bool use_metric, double metric_weight);
    void ImproveMesh(DomainIndex face_index, double metric_weight);
    void CombineImprove(DomainIndex face_index);

 protected:
    Mesh& mesh_;
};

double CalcTriangleBadness(const Point2d& p1, const Point2d& p2, const Point2d& p3,
                           double metric_weight, double h);

double CalcTriangleBadness_2(const Point2d& p1, const Point2d& p2, const Point2d& p3, double n_z);

}  // namespace meshit

#endif
