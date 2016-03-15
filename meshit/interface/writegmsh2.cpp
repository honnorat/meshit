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

namespace meshit {

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
        os << "2.2 0 "
        << (int) sizeof(double) << "\n";
        os << "$EndMeshFormat\n";

        /// Write nodes
        os << "$Nodes\n";
        os << np << "\n";
        int cnt = 0;
        for (PointIndex i = 0; i < np; i++) {
            const Point3d& p = mesh.Point(i);
            os << cnt++ << " ";  // number of nodes
            os << p.X() << " ";
            os << p.Y() << " ";
            os << p.Z() << "\n";
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
            os << cnt++ << " 1 2 " << seg.si << " " << seg.si << " " << seg[0] << " " << seg[1] << std::endl;
        }
        //            cnt += ns;
        for (size_t k = 0; k < nse; k++) {
            const Element2d& el = mesh.SurfaceElement(k);
            os << cnt++ << " 2 2 ";   // GMSH Type for a 3 node triangle
            os << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
            os << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
            os << el.PNum(1) << " ";
            os << el.PNum(2) << " ";
            os << el.PNum(3) << "\n";
        }
        os << "$EndElements\n";
        /*
         * End of 2D section
         */
    }
}  // namespace meshit
