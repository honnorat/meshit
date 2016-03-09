#ifndef MESHTOOL_HPP
#define MESHTOOL_HPP

#include "meshclass.hpp"

namespace meshit {

    void MeshQuality2d(const Mesh & mesh);

    class Surface;

    int CheckSurfaceMesh(const Mesh & mesh);
    int CheckSurfaceMesh2(const Mesh & mesh);

}
#endif
