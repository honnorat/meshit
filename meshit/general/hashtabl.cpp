/**************************************************************************/
/* File:   hashtabl.cpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type HASHTABLE
*/

#include <algorithm>
#include <stdexcept>

#include "hashtabl.hpp"

namespace meshit {

    std::ostream& operator<<(std::ostream& s, const INDEX_2& i2)
    {
        return s << i2.I1() << ", " << i2.I2();
    }

    std::ostream& operator<<(std::ostream& s, const INDEX_3& i3)
    {
        return s << i3.I1() << ", " << i3.I2() << ", " << i3.I3();
    }

    BASE_INDEX_2_CLOSED_HASHTABLE::
    BASE_INDEX_2_CLOSED_HASHTABLE(size_t size)
            : hash(size)
    {
        invalid = -1;
        for (size_t i = 0; i < size; i++) {
            hash[i].I1() = invalid;
        }
    }

    void BASE_INDEX_2_CLOSED_HASHTABLE::
    BaseSetSize(size_t size)
    {
        hash.resize(size);
        for (size_t i = 0; i < size; i++) {
            hash[i].I1() = invalid;
        }
    }

    bool BASE_INDEX_2_CLOSED_HASHTABLE::
    PositionCreate2(const INDEX_2& ind, size_t& apos)
    {
        size_t i = HashValue(ind);
        size_t startpos = i;
        while (true) {
            i++;
            if (i > hash.size()) i = 1;
            if (hash[i - 1] == ind) {
                apos = i;
                return false;
            }
            if (hash[i - 1].I1() == invalid) {
                hash[i - 1] = ind;
                apos = i;
                return true;
            }
            if (i == startpos) {
                throw std::runtime_error("Try to set new element in full closed hashtable");
            }
        }
    }

    size_t BASE_INDEX_2_CLOSED_HASHTABLE::UsedElements() const
    {
        size_t n = hash.size();
        size_t cnt = 0;
        for (size_t i = 0; i < n; i++) {
            if (hash[i].I1() != invalid) {
                cnt++;
            }
        }
        return cnt;
    }

    void BASE_INDEX_3_CLOSED_HASHTABLE::
    BaseSetSize(size_t size)
    {
        hash.resize(size);
        for (size_t i = 0; i < size; i++) {
            hash[i].I1() = invalid;
        }
    }

    bool BASE_INDEX_3_CLOSED_HASHTABLE::
    PositionCreate2(const INDEX_3& ind, size_t& apos)
    {
        size_t i = HashValue(ind);
        size_t startpos = i;
        while (true) {
            i = (i + 1) % hash.size();
            if (hash[i] == ind) {
                apos = i;
                return false;
            }
            if (hash[i].I1() == invalid) {
                hash[i] = ind;
                apos = i;
                return true;
            }
            if (i == startpos)
                throw std::runtime_error("Try to set new element in full closed hashtable");
        }
    }

}  // namespace meshit
