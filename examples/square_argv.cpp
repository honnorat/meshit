#include <iostream>
#include <sstream>

#include <meshit/meshing.hpp>
#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/geom2d/genmesh2d.hpp>
#include <meshit/meshing/meshtool.hpp>

int main(int argc, char ** argv) {

    std::cout << "MeshIt Square_argv" << std::endl;

    meshit::MeshingParameters mp;
    meshit::Mesh * mesh = nullptr;

    // creates geometry structure
    meshit::SplineGeometry2d geom;
    geom.Load(argv[1]);

    std::cout << "start meshing" << std::endl;

    mp.optsteps2d = 5;
    geom.SetGrading(2.0);
    meshit::MeshFromSpline2D(geom, mesh, mp);
    std::cout << "meshing done" << std::endl;
    meshit::MeshQuality2d(*mesh);
    mesh->Export(geom, "square_argv.msh", "Gmsh2 Format");
    mesh->Save("square_argv.meshit");
    delete mesh;

    return 0;
}
