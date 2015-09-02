//
//  Write Abaqus file
//
//

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "writeuser.hpp"




void WriteAbaqusFormat (const Mesh & mesh,
			const string & filename)

{
      
  std::cout << "\nWrite Abaqus Volume Mesh" <<std::endl;

  ofstream outfile (filename.c_str());

  outfile << "*Heading" <<std::endl;
  outfile << " " << filename <<std::endl;

  outfile.precision(8);

  outfile << "*Node" <<std::endl;

  int np = mesh.GetNP();
  int ne = mesh.GetNE();
  int i, j, k;

  for (i = 1; i <= np; i++)
    {
      outfile << i << ", ";
      outfile << mesh.Point(i)(0) << ", ";
      outfile << mesh.Point(i)(1) << ", ";
      outfile << mesh.Point(i)(2) << "\n";
    }

  int elemcnt = 0; //element counter
  int finished = 0;
  int indcnt = 1; //index counter

  while (!finished)
    {
      int actcnt = 0;
      const Element & el1 = mesh.VolumeElement(1);
      int non = el1.GetNP();
      if (non == 4)
	{
	  outfile << "*Element, type=C3D4, ELSET=PART" << indcnt <<std::endl;
	} 
      else if (non == 10)
	{
	  outfile << "*Element, type=C3D10, ELSET=PART" << indcnt <<std::endl;
	} 
      else
	{
	  std::cout << "unsupported Element type!!!" <<std::endl;	  
	}

      for (i = 1; i <= ne; i++)
	{
	  const Element & el = mesh.VolumeElement(i);
	      
	  if (el.GetIndex() == indcnt)
	    {
	      actcnt++;
	      if (el.GetNP() != non) 
		{
		  std::cout << "different element-types in a subdomain are not possible!!!" <<std::endl;
		  continue;
		}
		  
	      elemcnt++;
	      outfile << elemcnt << ", ";
	      if (non == 4)
		{
		  outfile << el.PNum(1) << ", ";
		  outfile << el.PNum(2) << ", ";
		  outfile << el.PNum(4) << ", ";
		  outfile << el.PNum(3) << "\n";
		}
	      else if (non == 10)
		{
		  outfile << el.PNum(1) << ", ";
		  outfile << el.PNum(2) << ", ";
		  outfile << el.PNum(4) << ", ";
		  outfile << el.PNum(3) << ", ";
		  outfile << el.PNum(5) << ", ";
		  outfile << el.PNum(9) << ", ";
		  outfile << el.PNum(7) << ", " << "\n";
		  outfile << el.PNum(6) << ", ";
		  outfile << el.PNum(8) << ", ";
		  outfile << el.PNum(10) << "\n";
		}
	      else
		{
		  std::cout << "unsupported Element type!!!" <<std::endl;
		  for (j = 1; j <= el.GetNP(); j++)
		    {
		      outfile << el.PNum(j);
		      if (j != el.GetNP()) outfile << ", ";
		    }
		  outfile << "\n";
		}
	    }
	}	  
      indcnt++;
      if (elemcnt == ne) {finished = 1; std::cout << "all elements found by Index!" <<std::endl;}
      if (actcnt == 0) {finished = 1;}
    }

  if (mesh.GetIdentifications().GetMaxNr())
    {
      // periodic identification, implementation for
      // Helmut J. Boehm, TU Vienna
	  
      char cfilename[255];
      strcpy (cfilename, filename.c_str());

      char mpcfilename[255];
      strcpy (mpcfilename, cfilename);
      size_t len = strlen (cfilename);
      if (len >= 4 && (strcmp (mpcfilename+len-4, ".inp") == 0))
	strcpy (mpcfilename+len-4, ".mpc");
      else
	strcat (mpcfilename, ".mpc");
	  
      ofstream mpc (mpcfilename);

      int masternode(0);

      Array<INDEX_2> pairs;
      BitArray master(np), help(np);
      master.Set();
      for (i = 1; i <= 3; i++)
	{
	  mesh.GetIdentifications().GetPairs (i, pairs);
	  help.Clear();
	  for (j = 1; j <= pairs.Size(); j++)
	    {
	      help.Set (pairs.Get(j).I1());
	    }
	  master.And (help);
	}
      for (i = 1; i <= np; i++)
	if (master.Test(i))
	  masternode = i;

      std::cout << "masternode = " << masternode << " = "
	   << mesh.Point(masternode) <<std::endl;
      Array<int> slaves(3);
      for (i = 1; i <= 3; i++)
	{
	  mesh.GetIdentifications().GetPairs (i, pairs);
	  for (j = 1; j <= pairs.Size(); j++)
	    {
	      if (pairs.Get(j).I1() == masternode)
		slaves.Elem(i) = pairs.Get(j).I2();
	    }
	  std::cout << "slave(" << i << ") = " << slaves.Get(i)
	       << " = " << mesh.Point(slaves.Get(i)) <<std::endl;
	}
	  
	  
      outfile << "**\n"
	      << "*NSET,NSET=CTENODS\n"
	      << slaves.Get(1) << ", " 
	      << slaves.Get(2) << ", " 
	      << slaves.Get(3) <<std::endl;

	  
      outfile << "**\n"
	      << "**POINT_fixed\n"
	      << "**\n"
	      << "*BOUNDARY, OP=NEW\n";
      for (j = 1; j <= 3; j++)
	outfile << masternode << ", " << j << ",,    0.\n";

      outfile << "**\n"
	      << "*BOUNDARY, OP=NEW\n";
      for (j = 1; j <= 3; j++)
	{
	  Vec3d v(mesh.Point(masternode), mesh.Point(slaves.Get(j)));
	  double vlen = v.Length();
	  int dir = 0;
	  if (fabs (v.X()) > 0.9 * vlen) dir = 2;
	  if (fabs (v.Y()) > 0.9 * vlen) dir = 3;
	  if (fabs (v.Z()) > 0.9 * vlen) dir = 1;
	  if (!dir)
	    std::cout << "ERROR: Problem with rigid body constraints" <<std::endl;
	  outfile << slaves.Get(j) << ", " << dir << ",,    0.\n";
	}

      outfile << "**\n"
	      << "*EQUATION, INPUT=" << mpcfilename <<std::endl;
	  

      BitArray eliminated(np);
      eliminated.Clear();
      for (i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++)
	{
	  mesh.GetIdentifications().GetPairs (i, pairs);
	  if (!pairs.Size())
	    continue;
	      
	  for (j = 1; j <= pairs.Size(); j++)
	    if (pairs.Get(j).I1() != masternode && 
		!eliminated.Test(pairs.Get(j).I2()))
	      {
		eliminated.Set (pairs.Get(j).I2());
		for (k = 1; k <= 3; k++)
		  {
		    mpc << "4" << "\n";
		    mpc << pairs.Get(j).I2() << "," << k << ", -1.0, ";
		    mpc << pairs.Get(j).I1() << "," << k << ", 1.0, ";
		    mpc << slaves.Get(i) << "," << k << ", 1.0, ";
		    mpc << masternode << "," << k << ", -1.0 \n";
		  }
	      }
	}
    }


  std::cout << "done" <<std::endl;
}

}
