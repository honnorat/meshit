#include <iostream>
#include <sstream>

#include <meshit/meshing/meshclass.hpp>
#include <meshit/meshing/meshtool.hpp>
#include <meshit/geom2d/geometry2d.hpp>

int main(int argc, char ** argv) {

    std::cout << "MeshIt Square_argv" << std::endl;

    meshit::MeshingParameters mp;
    meshit::Mesh mesh;

    // creates geometry structure
    meshit::SplineGeometry2d geom;
    geom.Load(argv[1]);

    std::cout << "start meshing" << std::endl;

    mp.optsteps2d = 1;
    mp.quad = 0;
    geom.SetGrading(2.0);
    mesh.BuildFromSpline2D(geom, mp);
    std::cout << "meshing done" << std::endl;
    mesh.Export("square_argv.msh", "Gmsh2 Format");
    mesh.Save("square_argv.meshit");

    return 0;
}
