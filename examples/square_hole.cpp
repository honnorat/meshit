#include <iostream>
#include <sstream>
#include <vector>

#include <meshit/meshing/meshclass.hpp>
#include <meshit/meshing/meshtool.hpp>
#include <meshit/geom2d/geometry2d.hpp>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/consoleappender.h>

int main(int argc, char ** argv)
{
    log4cplus::initialize();
    log4cplus::SharedAppenderPtr console(new log4cplus::ConsoleAppender(false, true));
    log4cplus::Logger::getRoot().addAppender(console);
    log4cplus::Logger logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("main"));

    LOG4CPLUS_INFO(logger, "== MeshGen Square_hole");

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

    LOG4CPLUS_INFO(logger, "== start meshing");
    mp.optsteps2d = 5;
    geom.SetGrading(0.02);
    mesh.BuildFromSpline2D(geom, mp);
    LOG4CPLUS_INFO(logger, "== meshing done");
    meshit::MeshQuality2d(mesh);
    meshit::CheckSurfaceMesh(mesh);
    meshit::CheckSurfaceMesh2(mesh);
    meshit::RemoveProblem(mesh, 0);
    meshit::RemoveProblem(mesh, 1);
    mesh.CheckConsistentBoundary();
    LOG4CPLUS_INFO(logger, "== export mesh");
    mesh.Export("square_hole.msh", "Gmsh2 Format");
    mesh.Save("square_hole.meshit");

    return 0;
}
