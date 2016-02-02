#include <iostream>
#include <sstream>
#include <vector>

#include <meshit/meshit.hpp>
#include <meshit/meshing/meshclass.hpp>
#include <meshit/meshing/meshtool.hpp>
#include <meshit/geom2d/geometry2d.hpp>

int main(int argc, char ** argv)
{
    MESHIT_LOG_INFO("== MeshGen Square_hole");

    meshit::MeshingParameters mp;
    meshit::Mesh mesh;

    // creates geometry structure
    meshit::SplineGeometry2d geom;

    std::vector<meshit::Point2d> list_outer;
    list_outer.push_back(meshit::Point2d(-1, -1));
    list_outer.push_back(meshit::Point2d(1, -1));
    list_outer.push_back(meshit::Point2d(1, 1));
    list_outer.push_back(meshit::Point2d(-1, 1));
    geom.AddLine(list_outer, 0.25, false, 2);

    std::vector<meshit::Point2d> line_inner;
    line_inner.push_back(meshit::Point2d(-0.5, -0.5));
    line_inner.push_back(meshit::Point2d(-0.5, 0.5));
    line_inner.push_back(meshit::Point2d(0.5, 0.5));
    line_inner.push_back(meshit::Point2d(0.5, -0.5));
    geom.AddLine(line_inner, 0.05, true, 3);

    geom.FakeData();

    MESHIT_LOG_INFO("== start meshing");
    mp.optsteps2d = 5;
    geom.SetGrading(0.02);
    mesh.BuildFromSpline2D(geom, mp);
    MESHIT_LOG_INFO("== meshing done");
    meshit::MeshQuality2d(mesh);
    meshit::CheckSurfaceMesh(mesh);
    meshit::CheckSurfaceMesh2(mesh);
    mesh.CheckConsistentBoundary();
    MESHIT_LOG_INFO("== export mesh");
    mesh.Export("square_hole.msh", "Gmsh2 Format");
    mesh.Save("square_hole.meshit");

    return 0;
}
