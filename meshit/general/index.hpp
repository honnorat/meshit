#ifndef FILE_TEMPLATE_HPP
#define FILE_TEMPLATE_HPP

/**************************************************************************/
/* File:   template.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

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
