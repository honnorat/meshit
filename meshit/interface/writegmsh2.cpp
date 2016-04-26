/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2008 Philippose Rajan
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

#include "writeuser.hpp"

#include "../meshing/mesh_class.hpp"

namespace meshit {

/*! GMSH v2.xx mesh format export function
 */
void WriteGmsh2Format(const Mesh& mesh, const std::string& filename)
{
    std::ofstream ofs(filename.c_str());
    WriteGmsh2Format(mesh, ofs);
}

void WriteGmsh2Format(const Mesh& mesh, std::ostream& os)
{
    os.precision(6);
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.setf(std::ios::showpoint);

    size_t np = mesh.GetNbPoints();    // number of points in mesh
    size_t ns = mesh.GetNbSegments();  // number of boundary (edge) segments in mesh
    size_t ne = mesh.GetNbElements();  // number of 2d elements

    /// Prepare GMSH 2.2 file (See GMSH 2.2 Documentation)
    os << "$MeshFormat\n";
    os << "2.2 0 " << sizeof(double) << "\n";
    os << "$EndMeshFormat\n";

    /// Write nodes
    os << "$Nodes\n";
    os << np << "\n";
    int cnt = 0;
    for (size_t i = 0; i < np; i++) {
        const MeshPoint& p = mesh.Point(i);
        os << ++cnt << " ";  // number of nodes
        os << p.X() << " ";
        os << p.Y() << " ";
        os << 0.0 << "\n";
    }
    os << "$EndNodes\n";

    // Write triangles
    os << "$Elements\n";
    os << ne + ns << "\n";  // Edge segments are elements too.

    cnt = 0;
    for (size_t i = 0; i < ns; i++) {
        const Segment& seg = mesh.LineSegment(i);
        os << ++cnt << " 1";  // GMSH type for a segment
        os << " 2 " << seg.edge_id << " " << seg.edge_id << " ";
        os << seg[0] + 1 << " " << seg[1] + 1 << std::endl;
    }
    for (size_t k = 0; k < ne; k++) {
        const Element2d& el = mesh.Element(k);
        int domain = mesh.GetDomainNumber(el);
        os << ++cnt << " 2";  // GMSH type for a 3 node triangle
        os << " 2 " << domain << " " << domain << " ";
        os << el.PointID(0) + 1 << " ";
        os << el.PointID(1) + 1 << " ";
        os << el.PointID(2) + 1 << "\n";
    }
    os << "$EndElements\n";
}

}  // namespace meshit
