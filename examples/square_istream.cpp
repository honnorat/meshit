#include <iostream>
#include <sstream>

#include <meshit/meshing/meshclass.hpp>
#include <meshit/meshing/meshtool.hpp>
#include <meshit/geom2d/geometry2d.hpp>

int main(int argc, char ** argv) {

    std::cout << "MeshIt Square_istream" << std::endl;

    meshit::MeshingParameters mp;
    meshit::Mesh mesh;

    // creates geometry structure
    meshit::SplineGeometry2d geom;

    std::stringstream ss(std::ios::in | std::ios::out);
    ss << "2" << std::endl;
    ss << "points" << std::endl;
    ss << "1   -1   -1" << std::endl;
    ss << "2    1   -1" << std::endl;
    ss << "3    1    1" << std::endl;
    ss << "4   -1    1" << std::endl;
    ss << "segments" << std::endl;
    ss << "1    0    2    1    2  -bc=1 -maxh=0.05" << std::endl;
    ss << "1    0    2    2    3  -bc=2 " << std::endl;
    ss << "1    0    2    3    4  -bc=3 " << std::endl;
    ss << "1    0    2    4    1  -bc=4 " << std::endl;
    ss << "materials" << std::endl;
    ss << "1    domain1   -maxh=0.5" << std::endl;
    geom.LoadData(ss);

    std::cout << "start meshing" << std::endl;

    mp.optsteps2d = 5;
    geom.SetGrading(2.0);
    mesh.BuildFromSpline2D(geom, mp);
    std::cout << "meshing done" << std::endl;
    meshit::MeshQuality2d(mesh);
    mesh.Export("square_istream.msh", "Gmsh2 Format");
    mesh.Save("square_istream.meshit");

    return 0;
}
