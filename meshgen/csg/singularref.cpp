#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{

  SingularEdge :: SingularEdge (double abeta, int adomnr, 
                                const CSGeometry & ageom,
                                const Solid * asol1, 
                                const Solid * asol2, double sf,
                                const double maxh_at_initialization)
    : domnr(adomnr), geom(ageom)
  {
  beta = abeta;
  maxhinit = maxh_at_initialization;

  if (beta > 1) 
    {
      beta = 1;
      std::cout << "Warning: beta set to 1" <<std::endl;
    }
  if (beta <= 1e-3)
    {
      beta = 1e-3;
      std::cout << "Warning: beta set to minimal value 0.001" <<std::endl;
    }

  sol1 = asol1;
  sol2 = asol2;
  factor = sf; 
}

void SingularEdge :: FindPointsOnEdge (class Mesh & mesh)
{
  std::cerr << "find points on edge" <<std::endl;
  points.resize(0);
  segms.resize(0);


  Array<int> si1, si2;
  sol1->GetSurfaceIndices (si1);
  sol2->GetSurfaceIndices (si2);

  for (int i = 0; i < si1.size(); i++)
    si1[i] = geom.GetSurfaceClassRepresentant(si1[i]);
  for (int i = 0; i < si2.size(); i++)
    si2[i] = geom.GetSurfaceClassRepresentant(si2[i]);


  for (SegmentIndex si = 0; si < mesh.GetNSeg(); si++)
    {
      INDEX_2 i2 (mesh[si][0], mesh[si][1]);
      /*
      
      bool onedge = 1;
      for (j = 1; j <= 2; j++)
	{
	  const Point<3> p = mesh[ PointIndex (i2.I(j)) ];
	  if (sol1->IsIn (p, 1e-3) && sol2->IsIn(p, 1e-3) &&
	      !sol1->IsStrictIn (p, 1e-3) && !sol2->IsStrictIn(p, 1e-3))
	    {
	      ;
	    }
	  else
	    onedge = 0;
	}
      */

      if (domnr != -1 && domnr != mesh[si].domin && domnr != mesh[si].domout)
	continue;

      /*
      bool onedge = 1;
      for (int j = 0; j < 2; j++)
	{
	  int surfi = (j == 0) ? mesh[si].surfnr1 : mesh[si].surfnr2;
	  surfi = geom.GetSurfaceClassRepresentant(surfi);
	  if (!si1.Contains(surfi) && !si2.Contains(surfi))
	    onedge = 0;
	}
      */
      int surfi1 = geom.GetSurfaceClassRepresentant(mesh[si].surfnr1);
      int surfi2 = geom.GetSurfaceClassRepresentant(mesh[si].surfnr2);

      if ( (si1.Contains(surfi1) && si2.Contains(surfi2)) ||
           (si1.Contains(surfi2) && si2.Contains(surfi1)) )

	// if (onedge)
	{
	  segms.push_back (i2);
	  //	  PrintMessage (5, "sing segment ", i2.I1(), " - ", i2.I2());
	  points.push_back (mesh[ PointIndex (i2.I1())]);
	  points.push_back (mesh[ PointIndex (i2.I2())]);
	  mesh[si].singedge_left = factor;
	  mesh[si].singedge_right = factor;
	}	    
    }
  
  /*
  std::cerr << "Singular edge points:" <<std::endl;
  for (int i = 0; i < points.Size(); i++)
    std::cerr << points[i] <<std::endl;
  */
 
}

void SingularEdge :: SetMeshSize (class Mesh & mesh, double globalh)
{
  double hloc = pow (globalh, 1/beta);
  if(maxhinit > 0 && maxhinit < hloc)
    {
      hloc = maxhinit;
      if(points.size() > 1)
	{
	  for (int i = 0; i < points.size()-1; i++)
	    mesh.RestrictLocalHLine(points[i],points[i+1],hloc);
	}
      else
	{
	  for (int i = 0; i < points.size(); i++)
	    mesh.RestrictLocalH (points[i], hloc);
	}
    }
  else
    {
      for (int i = 0; i < points.size(); i++)
	mesh.RestrictLocalH (points[i], hloc);
    }
}



SingularPoint :: SingularPoint (double abeta, 
				const Solid * asol1, 
				const Solid * asol2,
				const Solid * asol3, double sf)
{
  beta = abeta;
  sol1 = asol1;
  sol2 = asol2;
  sol3 = asol3;
  factor = sf; 
}


void SingularPoint :: FindPoints (class Mesh & mesh)
{
  points.resize(0);
  Array<int> surfk, surf;


  for (PointIndex pi = PointIndex::BASE; 
       pi < mesh.GetNP()+PointIndex::BASE; pi++)
    {
      if (mesh[pi].Type() != FIXEDPOINT) continue;
      const Point<3> p = mesh[pi];

      std::cerr << "check singular point" << p <<std::endl;

      if (sol1->IsIn (p) && sol2->IsIn(p) && sol3->IsIn(p) &&
	  !sol1->IsStrictIn (p) && !sol2->IsStrictIn(p) && !sol3->IsStrictIn(p))
	{
	  surf.resize (0);
	  for (int k = 1; k <= 3; k++)
	    {
	      const Solid * solk(NULL);
	      Solid *tansol;
	      switch (k)
		{
		case 1:  solk = sol1; break;
		case 2:  solk = sol2; break;
		case 3:  solk = sol3; break;
		}

	      solk -> TangentialSolid (p, tansol, surfk, 1e-3);
	      std::cerr << "Tansol = " << *tansol <<std::endl;

	      if (!tansol) continue;

	      ReducePrimitiveIterator rpi(Box<3> (p-Vec<3> (1e-3,1e-3,1e-3),
						  p+Vec<3> (1e-3,1e-3,1e-3)));
	      UnReducePrimitiveIterator urpi;
	      
	      tansol -> IterateSolid (rpi);
	      tansol->GetSurfaceIndices (surfk);
	      tansol -> IterateSolid (urpi);

	      std::cerr << "surfinds = " << surfk <<std::endl;

	      for (int i = 0; i < surfk.size(); i++)
		if (!surf.Contains (surfk[i]))
		  surf.push_back (surfk[i]);
	      
	      delete tansol;
	    }

	  if (surf.size() < 3) continue;

	  points.push_back (p);
	  PrintMessage (5, "Point (", p(0), ", ", p(1), ", ", p(2), ") is singular");
	  mesh[pi].Singularity(factor);
	}
    }  
}


void SingularPoint :: SetMeshSize (class Mesh & mesh, double globalh)
{
  double hloc = pow (globalh, 1/beta);
  for (int i = 1; i <= points.size(); i++)
    mesh.RestrictLocalH (points.Get(i), hloc);  
}
}
