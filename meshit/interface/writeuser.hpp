#ifndef WRITEUSER_HPP
#define WRITEUSER_HPP

/**************************************************************************/
/* File:    writeuser.hh                                                  */
/* Authors: many                                                          */
/* Date:    10. Dec. 97                                                   */
/**************************************************************************/

#include <string>
#include <iostream>

namespace meshit
{
    class Mesh;

    class SplineGeometry2d;

    void WriteGmsh2Format(const Mesh& mesh, const std::string& filename);
    void WriteGmsh2Format(const Mesh& mesh, std::ostream& os);

    bool WriteUserFormat(const std::string& format, const Mesh& mesh, const std::string& filename);
    bool WriteUserFormat(const std::string& format, const Mesh& mesh, std::ostream& ostream);

}  // namespace meshit
#endif

