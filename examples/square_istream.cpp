#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/meshing/mesh_class.hpp>

#include <iostream>
#include <sstream>

inline double uclock(void)
{
    return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
}

int main(int argc, char** argv)
{
    MESHIT_LOG_INFO("MeshIt square_istream");

    meshit::MeshingParameters mp;
    meshit::Mesh mesh;

    // creates geometry structure
    meshit::SplineGeometry geom;

    double cc;

    cc = uclock();
    std::stringstream ss(std::ios::in | std::ios::out);
    ss << "2" << std::endl;
    ss << "points" << std::endl;
    ss << "1   -1   -1" << std::endl;
    ss << "2    1   -1" << std::endl;
    ss << "3    1    1" << std::endl;
    ss << "4   -1    1" << std::endl;
    ss << "5   -0.5   -0.5" << std::endl;
    ss << "6    0.5   -0.5" << std::endl;
    ss << "7    0.5    0.5" << std::endl;
    ss << "8   -0.5    0.5" << std::endl;
    ss << "segments" << std::endl;
    ss << "1    0    2    1    2  -id=1 -maxh=0.0125" << std::endl;
    ss << "1    0    2    2    3  -id=2 " << std::endl;
    ss << "1    0    2    3    4  -id=3 " << std::endl;
    ss << "1    0    2    4    1  -id=4 " << std::endl;
    ss << "2    1    2    5    6  -id=5 " << std::endl;
    ss << "2    1    2    6    7  -id=6 -maxh=0.007 " << std::endl;
    ss << "2    1    2    7    8  -id=7 " << std::endl;
    ss << "2    1    2    8    5  -id=8 " << std::endl;
    ss << "materials" << std::endl;
    ss << "1    domain1   -maxh=0.25" << std::endl;
    ss << "2    domain2   -maxh=0.0125" << std::endl;
    geom.LoadData(ss);
    MESHIT_LOG_INFO(" . geometry loaded  in " << uclock() - cc << " s.");

    cc = uclock();
    mp.optsteps2d = 5;
    geom.SetGrading(0.3);
    mesh.BuildFromSplineGeometry(geom, mp);
    MESHIT_LOG_INFO(" . mesh created     in " << uclock() - cc << " s.");

    cc = uclock();
    mesh.PrintQuality();
    MESHIT_LOG_INFO(" . mesh diagnostics in " << uclock() - cc << " s.");

    std::string filename_export = "square_istream.msh";
    std::string filename_save = "square_istream.meshit";

    cc = uclock();
    mesh.Export(filename_export);
    MESHIT_LOG_INFO(" . export to GMSH   in " << uclock() - cc << " s."
                                              << " -> '" << filename_export << "'");
    mesh.Save("square_istream.meshit");

    cc = uclock();
    mesh.Save(filename_save);
    MESHIT_LOG_INFO(" . save to MESHIT   in " << uclock() - cc << " s."
                                              << " -> '" << filename_save << "'");

    return 0;
}
