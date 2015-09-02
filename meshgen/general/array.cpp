#ifndef FILE_NGSTD_ArrayCPP
#define FILE_NGSTD_ArrayCPP
// necessary for SGI ????

/**************************************************************************/
/* File:   array.cpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type Array
*/

#include <meshgen.hpp>
#include <assert.h>
#include "array.hpp"


namespace netgen
{
  //using namespace netgen;

#ifdef NONE  
  void BASE_Array :: ReSize (int minsize, int elementsize)
  {
    std::cout << "resize, minsize = " << minsize <<std::endl;

    if (inc == -1)
      throw Exception ("Try to resize fixed size array");

    
    void * p;
    int nsize = (inc) ? allocsize + inc : 2 * allocsize;
    if (nsize < minsize) nsize = minsize;

    if (data)
      {
	p = new char [nsize * elementsize];
	
	int mins = (nsize < actsize) ? nsize : actsize; 
	memcpy (p, data, mins * elementsize);
	
	delete [] static_cast<char*> (data);
	data = p;
      }
    else
      {
	data = new char[nsize * elementsize];
      }
    
    allocsize = nsize;
    std::cout << "resize done" <<std::endl;
  }
  
  
  
  void BASE_Array :: RangeCheck (int i) const
  {
    if (i < 0 || i >= actsize)
      throw ArrayRangeException ();
  }
  
  void BASE_Array :: CheckNonEmpty () const
  {
    if (!actsize)
      {
	throw Exception ("Array should not be empty");
	//      std::cerr << "Array souldn't be empty";
      }
  }
#endif
}
#endif

