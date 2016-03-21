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
        size_t blocks;
        size_t size;
        void* freelist;
        std::vector<char*> bablocks;

     public:
        explicit BlockAllocator(size_t asize, size_t ablocks = 100);
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
