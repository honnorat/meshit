#ifndef MESHIT_TEMPLATE_HPP
#define MESHIT_TEMPLATE_HPP
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

#include <iostream>
#include <unordered_map>

namespace meshit {

typedef uint_fast32_t GenericIndex;

class IndexPair
{
 public:
    IndexPair() { }

    IndexPair(GenericIndex i1, GenericIndex i2)
        : idx_{i1, i2} { }

    IndexPair(const IndexPair& other)
    {
        idx_[0] = other.idx_[0];
        idx_[1] = other.idx_[1];
    }

    bool operator==(const IndexPair& other) const
    {
        return idx_[0] == other.idx_[0] && idx_[1] == other.idx_[1];
    }

    IndexPair Sort()
    {
        if (idx_[0] > idx_[1]) {
            std::swap(idx_[0], idx_[1]);
        }
        return *this;
    }

    GenericIndex& I1() { return idx_[0]; }
    GenericIndex& I2() { return idx_[1]; }
    const GenericIndex& I1() const { return idx_[0]; }
    const GenericIndex& I2() const { return idx_[1]; }

    GenericIndex& operator[](size_t j) { return idx_[j]; }
    const GenericIndex& operator[](size_t j) const { return idx_[j]; }

 protected:
    GenericIndex idx_[2];
};

template<class T> using IndexPair_map = std::unordered_map<IndexPair, T>;

}  // namespace meshit

namespace std {

template<>
struct hash<::meshit::IndexPair>
{
    size_t operator()(const ::meshit::IndexPair& idx) const { return idx.I1() ^ (idx.I2() << 16); }
};

}  // namespace std

#endif
