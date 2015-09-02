#include <iostream>
#include <sstream>

#include <meshit/meshing.hpp>
#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/geom2d/genmesh2d.hpp>
#include <meshit/meshing/meshtool.hpp>

int main(int argc, char ** argv) {

    std::cout << "MeshIt Square_Inline" << std::endl;

    meshit::MeshingParameters mp;
    meshit::Mesh * mesh = nullptr;

    meshit::SplineGeometry2d geom;

    std::vector<meshit::Point2d> list_outer;
    list_outer.push_back(meshit::Point2d(-1, -1));
    list_outer.push_back(meshit::Point2d( 1, -1));
    list_outer.push_back(meshit::Point2d( 1,  1));
    list_outer.push_back(meshit::Point2d(-1,  1));
    geom.AddLine(list_outer, 0.25);
    geom.FakeData();

    std::cout << "start meshing" << std::endl;

    mp.optsteps2d = 5;
    geom.SetGrading(2.0);
    meshit::MeshFromSpline2D(geom, mesh, mp);
    std::cout << "meshing done" << std::endl;
    meshit::MeshQuality2d(*mesh);
    mesh->Export(geom, "square_inline.msh", "Gmsh2 Format");
    mesh->Save("square_inline.meshit");
    delete mesh;

    return 0;
}
