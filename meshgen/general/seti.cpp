#include <meshgen.hpp>
#include "seti.hpp"

namespace netgen
{
  //using namespace netgen;

  IndexSet :: IndexSet (int maxind)
  {
    SetMaxIndex (maxind);
  }

  IndexSet :: ~IndexSet ()
  {
    Clear();
  }


  void IndexSet :: SetMaxIndex (int maxind)
  {
    if (maxind > flags.Size())
      {
	flags.SetSize (2 * maxind);
	flags.Clear();
      }
  }

  /*
    int IndexSet :: IsIn (int ind) const
    {
    return flags.Test (ind);
    }
  */

  /*
    void IndexSet :: Add (int ind)
    {
    if (ind > flags.Size())
    {
    std::cerr << "out of range" <<std::endl;
    exit (1);
    }

    if (!flags.Test(ind))
    {
    set.Append (ind);
    flags.Set (ind);
    }
    }
  */

  void IndexSet :: Del (int ind)
  {
    for (int i = 1; i <= set.size(); i++)
      if (set.Get(i) == ind)
	{
	  set.DeleteElement (ind);
	  break;
	}
    flags.Clear (ind);
  }

  void IndexSet :: Clear ()
  {
    for (int i = 1; i <= set.size(); i++)
      flags.Clear (set.Get(i));
    set.resize (0);
  }
}
