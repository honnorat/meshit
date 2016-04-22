#ifndef FILE_OPTMEM_HPP
#define FILE_OPTMEM_HPP

/**************************************************************************/
/* File:   optmem.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

#include <cstdlib>
#include <vector>

namespace meshit {

/**
    Optimized Memory allocation classes
*/
class BlockAllocator
{
 public:
    explicit BlockAllocator(size_t asize, size_t ablocks = 512);
    ~BlockAllocator();

    void* Alloc();

    void Free(void* p)
    {
        *reinterpret_cast<void**>(p) = freelist;
        freelist = p;
    }

 private:
    size_t blocks;
    size_t size;
    void* freelist;
    std::vector<char*> bablocks;
};

}  // namespace meshit

#endif
