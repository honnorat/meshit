#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/meshing/meshtool.hpp>

#include <ctime>

inline double uclock(void)
{
    return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
}

int main(int argc, char** argv)
{
    meshit::SetLogLevel(MESHIT_INFO_LOG_LEVEL);

    // creates geometry structure
    meshit::SplineGeometry2d geom;
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
    geom.AddSpline(ellipse, 0.2, ++bc_num, face2, face1);

    meshit::Mesh mesh;
    meshit::MeshingParameters mp;
    mp.optsteps2d = 5;
    mp.maxh = 0.2;

    double cc;

    cc = uclock();
    mesh.BuildFromSpline2D(geom, mp);
    MESHIT_LOG_INFO("~ meshing  : " << uclock() - cc << " s.");

    cc = uclock();
    mesh.FindOpenElements(0);
    mesh.FindOpenElements(face1);
    mesh.FindOpenElements(face2);
    MESHIT_LOG_INFO("~ findopen : " << uclock() - cc << " s.");

    cc = uclock();
    meshit::MeshQuality2d(mesh);
    MESHIT_LOG_INFO("~ quality  : " << uclock() - cc << " s.");

    cc = uclock();
    meshit::CheckSurfaceMesh2(mesh);
    MESHIT_LOG_INFO("~ checksurf: " << uclock() - cc << " s.");

    cc = uclock();
    mesh.CheckOverlappingBoundary();
    MESHIT_LOG_INFO("~ checkover: " << uclock() - cc << " s.");

    cc = uclock();
    mesh.CheckConsistentBoundary();
    MESHIT_LOG_INFO("~ checkcons: " << uclock() - cc << " s.");

    mesh.PrintMemInfo(std::cout);
    MESHIT_LOG_INFO("AverageH(0) = " << mesh.AverageH(0));
    MESHIT_LOG_INFO("AverageH(1) = " << mesh.AverageH(face1));
    MESHIT_LOG_INFO("AverageH(2) = " << mesh.AverageH(face2));

    cc = uclock();
    mesh.UpdateTopology();
    MESHIT_LOG_INFO("~ updatetop: " << uclock() - cc << " s.");

    cc = uclock();
    mesh.Export("circles.msh");
    MESHIT_LOG_INFO("~ export   : " << uclock() - cc << " s.");

    meshit::Refinement ref;
    ref.Refine(mesh);
    meshit::MeshQuality2d(mesh);
    meshit::CheckSurfaceMesh2(mesh);
    mesh.CheckOverlappingBoundary();
    mesh.CheckConsistentBoundary();
    mesh.PrintMemInfo(std::cout);
    MESHIT_LOG_INFO("AverageH(0) = " << mesh.AverageH(0));
    MESHIT_LOG_INFO("AverageH(1) = " << mesh.AverageH(face1));
    MESHIT_LOG_INFO("AverageH(2) = " << mesh.AverageH(face2));
    mesh.Export("circles_fine.msh");

    return 0;
}
