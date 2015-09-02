#ifndef FILE_SPECIALS
#define FILE_SPECIALS

/*

  Very special implementations ..
  
 */
#include "meshclass.hpp"

namespace netgen {

///
void CutOffAndCombine (Mesh & mesh, const Mesh & othermesh);

void HelmholtzMesh (Mesh & mesh);

}

#endif
