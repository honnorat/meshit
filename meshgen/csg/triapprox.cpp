#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <csg.hpp>


namespace netgen
{

  TriangleApproximation :: TriangleApproximation ()
  {
    ;
  }

  int TriangleApproximation :: 
  AddTriangle (const TATriangle & tri, bool invert)
  { 
    trigs.push_back (tri);
    if (invert)
      {
	trigs.Last()[1] = tri[2];
	trigs.Last()[2] = tri[1];
      }
    return trigs.size()-1;
  }


  void TriangleApproximation :: RemoveUnusedPoints ()
  {
    BitArray used(GetNP());
    Array<int> map (GetNP());
    int i, j;
    int cnt = 0;

    used.Clear();
    for (i = 0; i < GetNT(); i++)
      for (j = 0; j < 3; j++)
	used.Set (GetTriangle (i)[j]);

    for (i = 0; i < GetNP(); i++)
      if (used.Test(i))
	map[i] = cnt++;
  
    for (i = 0; i < GetNT(); i++)
      for (j = 0; j < 3; j++)
	trigs[i][j] = map[trigs[i][j]];

    for (i = 0; i < GetNP(); i++)
      if (used.Test(i))
	{
	  points[map[i]] = points[i];
	  normals[map[i]] = normals[i];
	}

    points.resize (cnt);
    normals.resize (cnt);
  }
}
