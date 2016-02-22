#ifndef FILE_MESHTOOL
#define FILE_MESHTOOL

#include "meshclass.hpp"

namespace meshit {

    void MeshQuality2d(const Mesh & mesh);

    void SaveEdges(const Mesh & mesh, const char * geomfile, double h, char * filename);

    void SaveSurfaceMesh(const Mesh & mesh, double h, char * filename);

    class Surface;

    int CheckSurfaceMesh(const Mesh & mesh);
    int CheckSurfaceMesh2(const Mesh & mesh);

}
#endif
