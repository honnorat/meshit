#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/meshing/meshtool.hpp>

#include <ctime>

inline double uclock(void)
{
    return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
}

int main(int argc, char** argv)
{
    meshit::SetLogLevel(MESHIT_DEBUG_LOG_LEVEL);

    // creates geometry structure
    meshit::SplineGeometry geom;
    int bc_num = 1;

    // outer boundary
    int face1 = geom.AddFace("face1");
    geom.AddCircle({0.0, 1.5}, 3.0, 1e99, ++bc_num, face1);

    // add hole
    geom.AddCircle({0.0, 0.0}, 1.0, 0.1, ++bc_num, 0, face1);

    // Add inclusion
    int face2 = geom.AddFace("face2");
    std::vector<meshit::Point2d> ellipse = {{0.0,  2.25},
                                            {-1.5, 2.25},
                                            {-1.5, 1.75},
                                            {-1.5, 1.25},
                                            {+0.0, 1.25},
                                            {+1.5, 1.25},
                                            {+1.5, 1.75},
                                            {+1.5, 2.25}};
    geom.AddSpline(ellipse, 0.05, ++bc_num, face2, face1);

    meshit::Mesh mesh;
    meshit::MeshingParameters mp;
    mp.optsteps2d = 5;
    mp.maxh = 0.2;

    double cc;
    cc = uclock();
    mesh.BuildFromSplineGeometry(geom, mp);
    MESHIT_LOG_INFO("~ meshing  : " << uclock() - cc << " s.");

    cc = uclock();
    meshit::MeshQuality2d(mesh);
    MESHIT_LOG_INFO("~ quality  : " << uclock() - cc << " s.");

    cc = uclock();
    mesh.CheckOverlappingBoundary();
    MESHIT_LOG_INFO("~ checkover: " << uclock() - cc << " s.");

    mesh.PrintMemInfo(std::cout);
    MESHIT_LOG_INFO("AverageH(0) = " << mesh.AverageH(0));
    MESHIT_LOG_INFO("AverageH(1) = " << mesh.AverageH(face1));
    MESHIT_LOG_INFO("AverageH(2) = " << mesh.AverageH(face2));

    cc = uclock();
    mesh.Export("circles.msh");
    MESHIT_LOG_INFO("~ export   : " << uclock() - cc << " s.");

    cc = uclock();
    mesh.Save("circles.meshit");
    MESHIT_LOG_INFO("~ save     : " << uclock() - cc << " s.");

    mesh.Refine();
    meshit::MeshQuality2d(mesh);
    meshit::CheckSurfaceMesh(mesh);
    mesh.CheckOverlappingBoundary();
    mesh.PrintMemInfo(std::cout);
    MESHIT_LOG_INFO("AverageH(0) = " << mesh.AverageH(0));
    MESHIT_LOG_INFO("AverageH(1) = " << mesh.AverageH(face1));
    MESHIT_LOG_INFO("AverageH(2) = " << mesh.AverageH(face2));
    mesh.Export("circles_fine.msh");

    return 0;
}
