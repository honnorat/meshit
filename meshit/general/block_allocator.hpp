#ifndef MESHIT_BLOCK_ALLOC_HPP
#define MESHIT_BLOCK_ALLOC_HPP
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
