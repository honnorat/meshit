#include <meshit.hpp>

#include "improve2.hpp"
#include "improve3.hpp"
#include "global.hpp"
#include "../general/ngexception.hpp"
#include "../linalg/opti.hpp"

namespace meshit
{


  void MeshOptimize2d :: ProjectBoundaryPoints(Array<int> & surfaceindex, 
					       const Array<Point<3>* > & from, Array<Point<3>* > & dest)
  {
    for(int i=0; i<surfaceindex.size(); i++)
      {
	if(surfaceindex[i] >= 0)
	  {
	    *dest[i] = *from[i];
	    ProjectPoint(surfaceindex[i],*dest[i]);
	  }
      }
      

  }

  void MeshOptimize2d :: ImproveVolumeMesh (Mesh & mesh)
  {
    
    if (!faceindex)
      {
	PrintMessage (3, "Smoothing");

	for (faceindex = 1; faceindex <= mesh.GetNFD(); faceindex++)
	  {
	    ImproveVolumeMesh (mesh);
	    if (multithread.terminate)
	      throw NgException ("Meshing stopped");
	  }
	faceindex = 0;
	return;
      }

    int i, j, k;
    SurfaceElementIndex sei;

    Array<SurfaceElementIndex> seia;
    mesh.GetSurfaceElementsOfFace (faceindex, seia);

    /*
    bool mixed = 0;
    for (i = 0; i < seia.Size(); i++)
      if (mesh[seia[i]].GetNP() != 3)
	{
	  mixed = 1;
	  break;
	}
    */

    int loci;
    double fact;
    bool moveisok;

    PointGeomInfo ngi;
    Point<3> origp;

    Vector x(3);

    Array<MeshPoint, PointIndex::BASE> savepoints(mesh.GetNP());

    Array<int, PointIndex::BASE> nelementsonpoint(mesh.GetNP());
    nelementsonpoint = 0;

    for (i = 0; i < seia.size(); i++)
      {
	const Element2d & el = mesh[seia[i]];
	for (j = 0; j < el.GetNP(); j++)
	  nelementsonpoint[el[j]]++;
      }


    TABLE<SurfaceElementIndex,PointIndex::BASE> elementsonpoint(nelementsonpoint);
    for (i = 0; i < seia.size(); i++)
      {
	const Element2d & el = mesh[seia[i]];
	for (j = 0; j < el.GetNP(); j++)
	  elementsonpoint.Add (el[j], seia[i]);
      }
    

    JacobianPointFunction pf(mesh.Points(),mesh.VolumeElements());



//     Opti2SurfaceMinFunction surfminf(mesh);
//     Opti2EdgeMinFunction edgeminf(mesh);
//     Opti2SurfaceMinFunctionJacobian surfminfj(mesh);

    OptiParameters par;
    par.maxit_linsearch = 8;
    par.maxit_bfgs = 5;

    int np = mesh.GetNP();
    int ne = mesh.GetNE();

    BitArray badnodes(np);
    badnodes.Clear();

    for (i = 1; i <= ne; i++)
      {
	const Element & el = mesh.VolumeElement(i);
	double bad = el.CalcJacobianBadness (mesh.Points());
	if (bad > 1)
	  for (j = 1; j <= el.GetNP(); j++)
	    badnodes.Set (el.PNum(j));
      }


    bool printeddot = 0;
    char plotchar = '.';
    int modplot = 1;
    if (mesh.GetNP() > 1000)
      {
	plotchar = '+';
	modplot = 10;
      }
    if (mesh.GetNP() > 10000)
      {
	plotchar = 'o';
	modplot = 100;
      }
    int cnt = 0;


    Array<SurfaceElementIndex> locelements(0);
    Array<int> locrots(0);

    for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
      {
	if (mesh[pi].Type() != SURFACEPOINT)
	  continue;

	if (multithread.terminate)
	  throw NgException ("Meshing stopped");
	
	int surfi(-1);

	if(elementsonpoint[pi].size() == 0)
	  continue;

	Element2d & hel = mesh[elementsonpoint[pi][0]];

	if(hel.GetIndex() != faceindex)
	  continue;

	cnt++;
	if (cnt % modplot == 0 && writestatus)
	  {
	    printeddot = 1;
	    PrintDot (plotchar);
	  }

		
	int hpi = 0;
	for (j = 1; j <= hel.GetNP(); j++)
	  if (hel.PNum(j) == pi)
	    {
	      hpi = j;
	      break;
	    }
	PointGeomInfo gi1 = hel.GeomInfoPi(hpi);
	
	locelements.resize(0);
	locrots.resize (0);
	
	for (j = 0; j < elementsonpoint[pi].size(); j++)
	  {
	    sei = elementsonpoint[pi][j];
	    const Element2d & bel = mesh[sei];
	    surfi = mesh.GetFaceDescriptor(bel.GetIndex()).SurfNr();
	    
	    locelements.push_back (sei);
	    
	    for (k = 1; k <= bel.GetNP(); k++)
	      if (bel.PNum(k) == pi)
		{
		  locrots.push_back (k);
		  break;
		}
	  }
	 

	double lh = mesh.GetH(mesh.Point(pi));
	par.typx = lh;

	pf.SetPointIndex(pi);
	
	x = 0;
	bool pok = (pf.Func (x) < 1e10); 
	
	if (pok)
	  {
	    BFGS (x, pf, par);
	    
	    origp = mesh[pi];
	    loci = 1;
	    fact = 1;
	    moveisok = false;
	
	    
	    //optimizer loop (if whole distance is not possible, move only a bit!!!!)
	    while (loci <= 5 && !moveisok)
	      {
		loci ++;
		mesh[pi](0) = origp(0) + x(0)*fact;
		mesh[pi](1) = origp(1) + x(1)*fact;
		mesh[pi](2) = origp(2) + x(2)*fact;
		fact = fact/2.;
	    
	    
		//std::cout << "origp " << origp << " newp " << mesh[pi];
	    
		ngi = gi1;
		moveisok = (ProjectPointGI (surfi, mesh[pi], ngi) != 0);

		//std::cout << " projected " << mesh[pi] <<std::endl;

		// point lies on same chart in stlsurface
		
		if (moveisok)
		  {
		    for (j = 0; j < locelements.size(); j++)
		      mesh[locelements[j]].GeomInfoPi(locrots[j]) = ngi;

		    //std::cout << "moved " << origp << " to " << mesh[pi] <<std::endl;
		  }
		else
		  {
		    mesh[pi] = origp;
		  }
	    
	      }
	  }
	else
	  {
	    std::cout << "el not ok (point " << pi << ": " << mesh[pi] << ")" <<std::endl;
	  }
      }

    if (printeddot)
      PrintDot ('\n');
  
    mesh.SetNextTimeStamp();
  }

  
}
