#ifndef FILE_OPTMEM_HPP
#define FILE_OPTMEM_HPP

/**************************************************************************/
/* File:   optmem.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

#include "../general/array.hpp"

namespace meshit {

/** 
    Optimized Memory allocation classes
*/

    class BlockAllocator
    {
     private:
        unsigned size, blocks;
        void* freelist;
        Array<char*> bablocks;

     public:
        explicit BlockAllocator(unsigned asize, unsigned ablocks = 100);
        ~BlockAllocator();

        void* Alloc();

        void Free(void* p)
        {
            *(void**) p = freelist;
            freelist = p;
        }
    };

}  // namespace meshit

#endif
