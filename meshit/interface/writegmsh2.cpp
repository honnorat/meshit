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

    enum GMSH_ELEMENTS
    {
        GMSH_TRIG = 2, GMSH_TRIG6 = 9,
        GMSH_QUAD = 3, GMSH_QUAD8 = 16
    };
    const int triGmsh[7] = {0, 1, 2, 3, 6, 4, 5};
    const int quadGmsh[9] = {0, 1, 2, 3, 4, 5, 8, 6, 7};

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
    void WriteGmsh2Format(
            const Mesh& mesh,
            const std::string& filename)
    {
        std::ofstream ofs(filename.c_str());
        WriteGmsh2Format(mesh, ofs);
    }

    void WriteGmsh2Format(
            const Mesh& mesh,
            std::ostream& os)
    {
        os.precision(6);
        os.setf(std::ios::fixed, std::ios::floatfield);
        os.setf(std::ios::showpoint);

        int np = mesh.GetNP();    // number of points in mesh
        int ns = mesh.GetNSeg();  // number of segments in mesh
        int nse = mesh.GetNSE();  // number of surface elements (BC)

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
        for (SegmentIndex i = 0; i < ns; i++) {
            const Segment& seg = mesh.LineSegment(i);
            os << cnt++ << " 1 2 " << seg.si << " " << seg.si << " " << seg[0] << " " << seg[1] << std::endl;
        }
        //            cnt += ns;
        for (SurfaceElementIndex k = 0; k < nse; k++) {
            const Element2d& el = mesh.SurfaceElement(k);

            int elType = 0;
            if (el.GetNP() == 3) elType = GMSH_TRIG;   // GMSH Type for a 3 node triangle
            if (el.GetNP() == 6) elType = GMSH_TRIG6;  // GMSH Type for a 6 node triangle
            if (el.GetNP() == 4) elType = GMSH_QUAD;   // GMSH Type for a 4 node quadrangle
            if (el.GetNP() == 8) elType = GMSH_QUAD8;  // GMSH Type for an 8 node quadrangle
            if (elType == 0) {
                std::cerr << " Invalid surface element type for Gmsh 2.0 2D-Mesh Export Format !\n";
                return;
            }

            os << cnt++ << " " << elType << " 2 ";
            os << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
            os << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty();
            for (size_t l = 1; l <= el.GetNP(); l++) {
                os << " ";
                if ((elType == GMSH_TRIG) || (elType == GMSH_TRIG6)) {
                    os << el.PNum(triGmsh[l]);
                } else if ((elType == GMSH_QUAD) || (elType == GMSH_QUAD8)) {
                    os << el.PNum(quadGmsh[l]);
                }
            }
            os << "\n";
        }
        os << "$EndElements\n";
        /*
         * End of 2D section
         */
    }
}  // namespace meshit
