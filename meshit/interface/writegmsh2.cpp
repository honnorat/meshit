/*! \file writegmsh2.cpp
 *  \brief Export Netgen Mesh in the GMSH v2.xx File format
 *  \author Philippose Rajan
 *  \date 02 November 2008
 *
 *  This function extends the export capabilities of
 *  Netgen to include the GMSH v2.xx File Format.
 *
 *  Current features of this function include:
 *
 *  1. Exports Triangles, Quadrangles and Tetrahedra \n
 *  2. Supports upto second order elements of each type
 *
 */
#include "writeuser.hpp"

#include "../meshing/meshclass.hpp"

namespace meshit
{
    /*! GMSH v2.xx mesh format export function
     *
     *  This function extends the export capabilities of
     *  Netgen to include the GMSH v2.xx File Format.
     *
     *  Current features of this function include:
     *
     *  1. Exports Triangles, Quadrangles and Tetrahedra \n
     *  2. Supports upto second order elements of each type
     *
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

        size_t np = mesh.GetNP();    // number of points in mesh
        size_t ns = mesh.GetNSeg();  // number of segments in mesh
        size_t nse = mesh.GetNSE();  // number of surface elements (BC)

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

        /*
         * 2D section : available for triangles and quadrangles
         *              upto 2nd Order
         */
        // Write triangles & quadrangles
        os << "$Elements\n";
        os << nse + ns << "\n";

        cnt = 0;
        for (size_t i = 0; i < ns; i++) {
            const Segment& seg = mesh.LineSegment(i);
            os << ++cnt << " 1";   // GMSH type for a segment
            os << " 2 " << seg.edge_id << " " << seg.edge_id << " ";
            os << seg[0] + 1 << " " << seg[1] + 1 << std::endl;
        }
        for (size_t k = 0; k < nse; k++) {
            const Element2d& el = mesh.Element(k);
            int domain = mesh.GetDomainNumber(el);
            os << ++cnt << " 2";   // GMSH type for a 3 node triangle
            os << " 2 " << domain << " " << domain << " ";
            os << el.PointID(0) + 1 << " ";
            os << el.PointID(1) + 1 << " ";
            os << el.PointID(2) + 1 << "\n";
        }
        os << "$EndElements\n";
        /*
         * End of 2D section
         */
    }
}  // namespace meshit
