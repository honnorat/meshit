#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/meshing/meshclass.hpp>
#include <meshit/meshing/meshtool.hpp>

#include <iostream>
#include <sstream>
#include <vector>

inline double uclock(void)
{
    return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
}

int main(int argc, char** argv)
{
    MESHIT_LOG_INFO("== MeshGen Square_hole");

    meshit::MeshingParameters mp;
    meshit::Mesh mesh;

    // creates geometry structure
    meshit::SplineGeometry geom;

    int face1 = geom.AddFace("face 1");

    std::vector<meshit::Point2d> list_outer;
    list_outer.push_back(meshit::Point2d(-1, -1));
    list_outer.push_back(meshit::Point2d(1, -1));
    list_outer.push_back(meshit::Point2d(1, 1));
    list_outer.push_back(meshit::Point2d(-1, 1));
    geom.AddLine(list_outer, 0.25, 2, face1);

    std::vector<meshit::Point2d> line_inner;
    line_inner.push_back(meshit::Point2d(-0.5, -0.5));
    line_inner.push_back(meshit::Point2d(0.5, -0.5));
    line_inner.push_back(meshit::Point2d(0.5, 0.5));
    line_inner.push_back(meshit::Point2d(-0.5, 0.5));
    geom.AddHole(line_inner, 0.05, 3, face1);

    MESHIT_LOG_INFO("== start meshing");
    mp.optsteps2d = 5;
    geom.SetGrading(0.1);
    mesh.BuildFromSplineGeometry(geom, mp);
    MESHIT_LOG_INFO("== meshing done");
    meshit::MeshQuality2d(mesh);
    meshit::CheckSurfaceMesh(mesh);
    mesh.CheckOverlappingBoundary();
    mesh.PrintMemInfo(std::cout);
    MESHIT_LOG_INFO(" - averageH(0) = " << mesh.AverageH(0));
    MESHIT_LOG_INFO(" - averageH(1) = " << mesh.AverageH(1));

    MESHIT_LOG_INFO("== export mesh to GMSH format");
    mesh.Export("square_hole.msh");

    MESHIT_LOG_INFO("== save mesh to MESHIT format");
    mesh.Save("square_hole.meshit");

    MESHIT_LOG_INFO("== refine mesh");
    mesh.Refine();
    meshit::MeshQuality2d(mesh);
    meshit::CheckSurfaceMesh(mesh);
    mesh.CheckOverlappingBoundary();
    mesh.PrintMemInfo(std::cout);
    MESHIT_LOG_INFO(" - averageH(0) = " << mesh.AverageH(0));
    MESHIT_LOG_INFO(" - averageH(1) = " << mesh.AverageH(1));

    MESHIT_LOG_INFO("== export mesh to GMSH format");
    mesh.Export("square_hole_fine.msh");
    MESHIT_LOG_INFO("== save mesh to MESHIT format");
    mesh.Save("square_hole_fine.meshit");

    return 0;
}
