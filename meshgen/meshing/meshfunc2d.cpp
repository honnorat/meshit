#include <meshgen.hpp>
#include "improve2.hpp"
#include "meshfunc.hpp"
#include "global.hpp"

namespace netgen
{

  DLL_HEADER void Optimize2d (Mesh & mesh, MeshingParameters & mp)
  {
    mesh.CalcSurfacesOfNode();

    const char * optstr = mp.optimize2d;
    int optsteps = mp.optsteps2d;

    for (int i = 1; i <= optsteps; i++)
      for (size_t j = 1; j <= strlen(optstr); j++)
	{
	  if (multithread.terminate) break;
	  switch (optstr[j-1])
	    {
	    case 's': 
	      {  // topological swap
		MeshOptimize2d meshopt;
		meshopt.SetMetricWeight (0);
		meshopt.EdgeSwapping (mesh, 0);
		break;
	      }
	    case 'S': 
	      {  // metric swap
		MeshOptimize2d meshopt;
		meshopt.SetMetricWeight (0);
		meshopt.EdgeSwapping (mesh, 1);
		break;
	      }
	    case 'm': 
	      {
		MeshOptimize2d meshopt;
		meshopt.SetMetricWeight (1);
		meshopt.ImproveMesh(mesh, mp);
		break;
	      }
	    case 'c': 
	      {
		MeshOptimize2d meshopt;
		meshopt.SetMetricWeight (0.2);
		meshopt.CombineImprove(mesh);
		break;
	      }
	    default:
	      std::cerr << "Optimization code " << optstr[j-1] << " not defined" << std::endl;
	    }  
	}
  }

}
