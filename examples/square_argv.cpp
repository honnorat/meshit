#include <iostream>
#include <sstream>

#include <meshgen/meshing/meshclass.hpp>
#include <meshgen/meshing/meshtool.hpp>
#include <meshgen/geom2d/geometry2d.hpp>

int main(int argc, char ** argv) {

    std::cout << "MeshIt Square_argv" << std::endl;

    netgen::MeshingParameters mp;
    netgen::Mesh mesh;

    // creates geometry structure
    netgen::SplineGeometry2d geom;
    geom.Load(argv[1]);

    std::cout << "start meshing" << std::endl;

    mp.optsteps2d = 5;
    geom.SetGrading(2.0);
    mesh.BuildFromSpline2D(geom, mp);
    std::cout << "meshing done" << std::endl;
    mesh.Export("square_argv.msh", "Gmsh2 Format");
    mesh.Save("square_argv.meshit");

    return 0;
}
