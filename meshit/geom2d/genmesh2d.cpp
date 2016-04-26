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
#include <sstream>

#include "../gprim/geom3d.hpp"
#include "../meshing/mesh_optimize.hpp"
#include "../meshing/mesh_generator.hpp"
#include "geometry2d.hpp"

namespace meshit {

void Mesh::Optimize2d(MeshingParameters& mp)
{
    IndexBoundaryEdges();

    const char* optstr = mp.optimize2d;
    int optsteps = mp.optsteps2d;
    MeshOptimize meshopt(*this);

    for (int i = 1; i <= optsteps; i++) {
        for (size_t j = 0; j < strlen(optstr); j++) {
            switch (optstr[j]) {
            case 's': {  // topological swap
                meshopt.EdgeSwapping(false);
                mp.n_steps++;
                break;
            }
            case 'S': {  // metric swap
                meshopt.EdgeSwapping(true, 0.0);
                mp.n_steps++;
                break;
            }
            case 'm': {
                meshopt.ImproveMesh(1.0);
                mp.n_steps++;
                break;
            }
            case 'c': {
                meshopt.CombineImprove();
                mp.n_steps++;
                break;
            }
            case 'p': {
                // print mesh
                std::stringstream mesh_name;
                mesh_name << "mesh_debug_" << std::setfill('0') << std::setw(2) << mp.n_steps << ".msh";
                Export(mesh_name.str());
                break;
            }
            default:
                MESHIT_LOG_ERROR("Optimization code " << optstr[j] << " not defined");
            }
        }
    }
}

// partitionizes spline curve

void SplineGeometry::PartitionBoundary(Mesh& mesh2d, MeshingParameters& mp)
{
    // mesh size restrictions ...
    for (size_t i = 0; i < splines.size(); i++) {
        const SplineSeg& spline = *splines[i];
        const GeomPoint& p1 = spline.StartPI();
        const GeomPoint& p2 = spline.EndPI();

        double h1 = std::min(p1.hmax, mp.maxh / p1.refatpoint);
        double h2 = std::min(p2.hmax, mp.maxh / p2.refatpoint);
        mesh2d.RestrictLocalH(p1, h1);
        mesh2d.RestrictLocalH(p2, h2);

        if (mp.segments_per_edge > 1) {
            double len = spline.Length();
            mesh2d.RestrictLocalHLine(p1, p2, len / mp.segments_per_edge);
        }

        double hcurve = std::min(spline.hmax_, mp.maxh / spline.ref_fac_);
        double hl = GetDomainMaxh(spline.dom_left);
        double hr = GetDomainMaxh(spline.dom_right);
        if (hl > 0.0) hcurve = std::min(hcurve, hl);
        if (hr > 0.0) hcurve = std::min(hcurve, hr);

        uint32_t np = 1000;
        for (double t = 0.5 / np; t < 1; t += 1.0 / np) {
            Point2d x = spline.GetPoint(t);
            double hc = 1.0 / mp.curvature_safety / (1e-99 + spline.CalcCurvature(t));
            mesh2d.RestrictLocalH(x, std::min(hc, hcurve));
        }
    }

    Box2d bbox;
    GetBoundingBox(bbox);
    SplineSegmenter sseg(mesh2d, mp, elto0, bbox);
    sseg.Partition(splines);
}

}  // namespace meshit
