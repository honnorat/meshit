/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2015-2016 Marc Honnorat <marc.honnorat@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library in the file LICENSE.LGPL; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

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
