#ifndef FILE_SETI
#define FILE_SETI


/**************************************************************************/
/* File:   seti.hh                                                        */
/* Author: Joachim Schoeberl                                              */
/* Date:   20. Mar. 98                                                    */
/**************************************************************************/

#include "bitarray.hpp"

namespace meshit
{

/**
  Set of Integers
  */
class IndexSet
{
  Array<int> set;
  BitArray flags;
public:
  IndexSet (int maxind);
  
  ~IndexSet ();
  /// increase range to maxind
  void SetMaxIndex (int maxind);
  int IsIn (int ind) const
  { 
    return flags.Test (ind); 
  }

  void Add (int ind)
  {
    if (!flags.Test(ind))
      {
	set.push_back (ind);
	flags.Set (ind);
      }
  }

  void Del (int ind);
  void Clear ();
  
  const Array<int> & GetArray() { return set; }
};

}

#endif
