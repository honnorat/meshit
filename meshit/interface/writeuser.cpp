//
//  Write user dependent output file
//

#include <iostream>
#include "writeuser.hpp"

namespace meshit {
    bool WriteUserFormat(const std::string& format,
                         const Mesh& mesh,
                         const std::string& filename)
    {
        if (format == "Gmsh2 Format") {
            WriteGmsh2Format(mesh, filename);
        } else {
            return 1;
        }
        return 0;
    }

    bool WriteUserFormat(const std::string& format,
                         const Mesh& mesh,
                         std::ostream& os)
    {
        if (format == "Gmsh2 Format") {
            WriteGmsh2Format(mesh, os);
        } else {
            return 1;
        }
        return 0;
    }
}  // namespace meshit

