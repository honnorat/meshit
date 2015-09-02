//
//  Write user dependent output file
//

#include <iostream>
#include <meshgen.hpp>
#include "writeuser.hpp"

namespace netgen
{
    bool WriteUserFormat(
            const std::string & format,
            const Mesh & mesh,
            const std::string & filename) {

        std::cout << "Export mesh to file '" << filename << "', format is " << format << std::endl;

        if (format == "Gmsh2 Format")
            WriteGmsh2Format(mesh, filename);
        else {
            return 1;
        }
        return 0;
    }
}

