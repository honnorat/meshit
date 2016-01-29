/**************************************************************************/
/* File:   optmem.cpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

/* 
   Abstract data type Array
*/

#include "../meshit.hpp"
#include "optmem.hpp"

namespace meshit
{

  BlockAllocator :: BlockAllocator (unsigned asize, unsigned ablocks)
    : bablocks (0)
  {
    if (asize < sizeof(void*))
      asize = sizeof(void*);
    size = asize;
    blocks = ablocks;
    freelist = NULL;
  }

  BlockAllocator :: ~BlockAllocator ()
  {
    for (int i = 0; i < bablocks.size(); i++)
      delete [] bablocks[i];
  }

  void * BlockAllocator :: Alloc ()
  {
    //  return new char[size];
    if (!freelist)
      {
	// std::cout << "freelist = " << freelist <<std::endl;
	// std::cout << "BlockAlloc: " << size*blocks <<std::endl;
	char * hcp = new char [size * blocks];
	bablocks.push_back (hcp);
	bablocks.Last() = hcp;
	for (unsigned i = 0; i < blocks-1; i++)
	  *(void**)&(hcp[i * size]) = &(hcp[ (i+1) * size]);
	*(void**)&(hcp[(blocks-1)*size]) = NULL;
	freelist = hcp;
      }

    void * p = freelist;
    freelist = *(void**)freelist;
    return p;
  }

  /*
  void BlockAllocator :: Free (void * p)
  {
    *(void**)p = freelist;
    freelist = p;
  }
  */
}
