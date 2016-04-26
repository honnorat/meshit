#include <meshit/geom2d/geometry2d.hpp>
#include <meshit/meshing/mesh_class.hpp>

#include <sstream>
#include <string>

inline double uclock(void)
{
    return static_cast<double>(std::clock()) / static_cast<double>(CLOCKS_PER_SEC);
}

inline bool ends_with(const std::string& str, const std::string& end)
{
    if (end.size() > str.size()) return false;

    return std::equal(end.rbegin(), end.rend(), str.rbegin());
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        MESHIT_LOG_FATAL("No input file\n");
        MESHIT_LOG_INFO("Usage : " << argv[0] << " FILENAME");
        MESHIT_LOG_INFO("  where FILENAME is either :");
        MESHIT_LOG_INFO("   - a 2d geometry (.in2d) file, or");
        MESHIT_LOG_INFO("   - a MESHIT (.meshit) save file.");
        exit(1);
    }

    std::string filename = argv[1];
    std::string filename_export;
    std::string filename_save;

    MESHIT_LOG_INFO("MeshIt read_argv with '" << filename << "'");

    meshit::MeshingParameters mp;
    meshit::Mesh mesh;
    double cc;

    if (ends_with(filename, ".in2d")) {
        MESHIT_LOG_INFO("File '" << filename << "' is a geometry file. Let's try to mesh it !");

        // creates geometry structure
        cc = uclock();
        meshit::SplineGeometry geom;
        geom.Load(filename);
        MESHIT_LOG_INFO(" . geometry read in " << uclock() - cc << " s.");

        cc = uclock();
        mp.optsteps2d = 5;
        mesh.BuildFromSplineGeometry(geom, mp);
        MESHIT_LOG_INFO(" . meshing done  in " << uclock() - cc << " s.");

        filename_export = "square_argv.msh";
        filename_save = "square_argv.meshit";
    } else if (ends_with(filename, ".meshit")) {
        MESHIT_LOG_INFO("File '" << filename << "' is a save file in MESHIT format. Let's try to load it !");

        cc = uclock();
        mesh.Load(filename);
        MESHIT_LOG_INFO(" . mesh loaded in " << uclock() - cc << " s.");

        filename_export = "square_argv_copy.msh";
        filename_save = "square_argv_copy.meshit";
    } else {
        MESHIT_LOG_INFO("File format not recognized for'"
                        << filename
                        << "'. I will handle only MESHIT (.meshit) and 2d geometry (.in2d) formats.");
        return 1;
    }

    cc = uclock();
    mesh.Save(filename_save);
    MESHIT_LOG_INFO(" . save to MESHIT done in " << uclock() - cc << " s."
                                                 << " -> '" << filename_save << "'");

    cc = uclock();
    mesh.Export(filename_export);
    MESHIT_LOG_INFO(" . export to GMSH done in " << uclock() - cc << " s."
                                                 << " -> '" << filename_export << "'");

    return 0;
}
