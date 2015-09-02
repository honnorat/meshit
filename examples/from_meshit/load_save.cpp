#include <iostream>
#include <sstream>

#include <meshit/meshing.hpp>
#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/geom2d/genmesh2d.hpp>
#include <meshit/meshing/meshtool.hpp>

int main(int argc, char ** argv) {

    std::cout << "MeshIt Load_Save" << std::endl;

    meshit::Mesh * mesh = new meshit::Mesh;
    meshit::SplineGeometry2d geom;
    mesh->Load(argv[1]);
    meshit::RemoveProblem(*mesh, 0);
    meshit::RemoveProblem(*mesh, 1);
    mesh->Save(std::string(argv[1]) + "-copy");
    mesh->Export(geom, "load_save.msh", "Gmsh2 Format");
    delete mesh;
    
    return 0;
}
