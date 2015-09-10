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
#include <meshit.hpp>
#include "writeuser.hpp"

#include "../meshing/meshclass.hpp"

namespace meshit {

    enum GMSH_ELEMENTS
    {
        GMSH_TRIG = 2, GMSH_TRIG6 = 9,
        GMSH_QUAD = 3, GMSH_QUAD8 = 16,
        GMSH_TET = 4, GMSH_TET10 = 11
    };
    const int triGmsh[7] = {0, 1, 2, 3, 6, 4, 5};
    const int quadGmsh[9] = {0, 1, 2, 3, 4, 5, 8, 6, 7};
    const int tetGmsh[11] = {0, 1, 2, 3, 4, 5, 8, 6, 7, 10, 9};

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
            const Mesh & mesh,
            const std::string & filename)
    {
        std::ofstream outfile(filename.c_str());
        outfile.precision(6);
        outfile.setf(std::ios::fixed, std::ios::floatfield);
        outfile.setf(std::ios::showpoint);

        int np = mesh.GetNP(); /// number of points in mesh
        int ns = mesh.GetNSeg(); /// number of segments in mesh
        int nse = mesh.GetNSE(); /// number of surface elements (BC)

        /// Prepare GMSH 2.2 file (See GMSH 2.2 Documentation)
        outfile << "$MeshFormat\n";
        outfile << "2.2 0 "
                << (int) sizeof (double) << "\n";
        outfile << "$EndMeshFormat\n";

        /// Write nodes
        outfile << "$Nodes\n";
        outfile << np << "\n";
        int cnt = PointIndex::BASE;
        for (int i = 1; i <= np; i++) {
            const Point3d & p = mesh.Point(i);
            outfile << cnt++ << " "; /// node number
            outfile << p.X() << " ";
            outfile << p.Y() << " ";
            outfile << p.Z() << "\n";
        }
        outfile << "$EndNodes\n";

        /*
         * 2D section : available for triangles and quadrangles
         *              upto 2nd Order
         */
        /// write triangles & quadrangles
        outfile << "$Elements\n";
        outfile << nse + ns << "\n";

        cnt = PointIndex::BASE;
        for (int i = 1; i <= ns; i++) {
            const Segment & seg = mesh.LineSegment(i);
            outfile << cnt++ << " 1 2 " << seg.si << " " << seg.si << " " << seg[0] << " " << seg[1] << std::endl;
        }
        //            cnt += ns;
        for (int k = 1; k <= nse; k++) {
            const Element2d & el = mesh.SurfaceElement(k);

            int elType = 0;
            if (el.GetNP() == 3) elType = GMSH_TRIG; //// GMSH Type for a 3 node triangle
            if (el.GetNP() == 6) elType = GMSH_TRIG6; //// GMSH Type for a 6 node triangle
            if (el.GetNP() == 4) elType = GMSH_QUAD; //// GMSH Type for a 4 node quadrangle
            if (el.GetNP() == 8) elType = GMSH_QUAD8; //// GMSH Type for an 8 node quadrangle
            if (elType == 0) {
                std::cerr << " Invalid surface element type for Gmsh 2.0 2D-Mesh Export Format !\n";
                return;
            }

            outfile << cnt++ << " " << elType << " 2 ";
            outfile << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty() << " ";
            outfile << mesh.GetFaceDescriptor(el.GetIndex()).BCProperty();
            for (int l = 1; l <= el.GetNP(); l++) {
                outfile << " ";
                if ((elType == GMSH_TRIG) || (elType == GMSH_TRIG6)) {
                    outfile << el.PNum(triGmsh[l]);
                }
                else if ((elType == GMSH_QUAD) || (elType == GMSH_QUAD8)) {
                    outfile << el.PNum(quadGmsh[l]);
                }
            }
            outfile << "\n";
        }
        outfile << "$EndElements\n";
        /*
         * End of 2D section
         */
    }
}
