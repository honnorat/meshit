#include <iostream>
#include <sstream>
#include <vector>

#include <meshit/meshing/meshclass.hpp>
#include <meshit/meshing/meshtool.hpp>
#include <meshit/geom2d/geometry2d.hpp>

int main(int argc, char ** argv)
{
    std::cout << "MeshGen Square_hole" << std::endl;

    meshit::MeshingParameters mp;
    meshit::Mesh mesh;

    // creates geometry structure
    meshit::SplineGeometry2d geom;

    std::vector<meshit::Point2d> list_outer;
    list_outer.push_back(meshit::Point2d(-1, -1));
    list_outer.push_back(meshit::Point2d(1, -1));
    list_outer.push_back(meshit::Point2d(1, 1));
    list_outer.push_back(meshit::Point2d(-1, 1));
    geom.AddLine(list_outer, 0.25);

//    std::vector<meshit::Point2d> line_inner;
//    line_inner.push_back(meshit::Point2d(-0.5, -0.5));
//    line_inner.push_back(meshit::Point2d(-0.5, 0.5));
//    line_inner.push_back(meshit::Point2d(0.5, 0.5));
//    line_inner.push_back(meshit::Point2d(0.5, -0.5));
//    geom.AddLine(line_inner, 0.05, true, 2);
//
    geom.FakeData();

    std::cout << "start meshing" << std::endl;
//
    mp.optsteps2d = 5;
    geom.SetGrading(0.02);
    mesh.BuildFromSpline2D(geom, mp);
    std::cout << "meshing done" << std::endl;
    meshit::MeshQuality2d(mesh);
    meshit::CheckSurfaceMesh(mesh);
    meshit::CheckSurfaceMesh2(mesh);
    meshit::RemoveProblem(mesh, 0);
    meshit::RemoveProblem(mesh, 1);
    mesh.CheckConsistentBoundary();
    mesh.Export("square_hole.msh", "Gmsh2 Format");
    mesh.Save("square_hole.meshit");

    return 0;
}
