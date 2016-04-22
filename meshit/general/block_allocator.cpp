/**************************************************************************/
/* File:   optmem.cpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

/*
   Abstract data type Array
*/

#include "block_allocator.hpp"

namespace meshit {

BlockAllocator::BlockAllocator(size_t asize, size_t ablocks)
    : blocks{ablocks}, freelist{nullptr}
{
    if (asize < sizeof(void*)) asize = sizeof(void*);
    size = asize;
}

BlockAllocator::~BlockAllocator()
{
    for (size_t i = 0; i < bablocks.size(); i++) {
        delete[] bablocks[i];
    }
}

void* BlockAllocator::Alloc()
{
    if (!freelist) {
        char* hcp = new char[size * blocks];
        bablocks.push_back(hcp);
        for (size_t i = 0; i < blocks - 1; i++) {
            *reinterpret_cast<void**>(&(hcp[i * size])) = &(hcp[(i + 1) * size]);
        }
        *reinterpret_cast<void**>(&(hcp[(blocks - 1) * size])) = nullptr;
        freelist = hcp;
    }

    void* p = freelist;
    freelist = *static_cast<void**>(freelist);
    return p;
}

}  // namespace meshit
