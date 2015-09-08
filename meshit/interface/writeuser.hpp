#ifndef WRITEUSER
#define WRITEUSER

/**************************************************************************/
/* File:    writeuser.hh                                                  */
/* Authors: many                                                          */
/* Date:    10. Dec. 97                                                   */
/**************************************************************************/

#include <string>

namespace meshit {

    class Mesh;
    class SplineGeometry2d;
    
    void WriteGmsh2Format(const Mesh & mesh,
            const std::string & filename);

    bool WriteUserFormat(
            const std::string & format,
            const Mesh & mesh,
            const std::string & filename);

}
#endif
