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
    //using namespace netgen;

    void INDEX_4::Sort()
    {
        if (i[0] > i[1]) std::swap(i[0], i[1]);
        if (i[2] > i[3]) std::swap(i[2], i[3]);
        if (i[0] > i[2]) std::swap(i[0], i[2]);
        if (i[1] > i[3]) std::swap(i[1], i[3]);
        if (i[1] > i[2]) std::swap(i[1], i[2]);
    }


    void INDEX_4Q::Sort()
    {
        if (std::min(i[1], i[2]) < std::min(i[0], i[3])) {
            std::swap(i[0], i[1]);
            std::swap(i[2], i[3]);
        }
        if (i[3] < i[0]) {
            std::swap(i[0], i[3]);
            std::swap(i[1], i[2]);
        }
        if (i[3] < i[1]) { std::swap(i[1], i[3]); }
    }


    std::ostream& operator<<(std::ostream& s, const INDEX_2& i2)
    {
        return s << i2.I1() << ", " << i2.I2();
    }

    std::ostream& operator<<(std::ostream& s, const INDEX_3& i3)
    {
        return s << i3.I1() << ", " << i3.I2() << ", " << i3.I3();
    }

    std::ostream& operator<<(std::ostream& s, const INDEX_4& i4)
    {
        return s << i4.I1() << ", " << i4.I2() << ", " << i4.I3() << ", " << i4.I4();
    }

    std::ostream& operator<<(std::ostream& s, const INDEX_4Q& i4)
    {
        return s << i4.I1() << ", " << i4.I2() << ", " << i4.I3() << ", " << i4.I4();
    }

    BASE_INDEX_2_CLOSED_HASHTABLE::
    BASE_INDEX_2_CLOSED_HASHTABLE(int size)
            : hash(size)
    {
        // hash.SetName ("i2-hashtable, hash");

        invalid = -1;
        for (int i = 1; i <= size; i++)
            hash.Elem(i).I1() = invalid;
    }

    void BASE_INDEX_2_CLOSED_HASHTABLE::
    BaseSetSize(int size)
    {
        hash.resize(size);
        for (int i = 1; i <= size; i++)
            hash.Elem(i).I1() = invalid;
    }

    int BASE_INDEX_2_CLOSED_HASHTABLE::
    PositionCreate2(const INDEX_2& ind, int& apos)
    {
        int i = HashValue(ind);
        int startpos = i;
        while (1) {
            i++;
            if (i > hash.size()) i = 1;
            if (hash.Get(i) == ind) {
                apos = i;
                return 0;
            }
            if (hash.Get(i).I1() == invalid) {
                hash.Elem(i) = ind;
                apos = i;
                return 1;
            }
            if (i == startpos)
                throw std::runtime_error("Try to set new element in full closed hashtable");
        }
    }

    int BASE_INDEX_2_CLOSED_HASHTABLE::UsedElements() const
    {
        int n = hash.size();
        int cnt = 0;
        for (int i = 1; i <= n; i++)
            if (hash.Get(i).I1() != invalid)
                cnt++;
        return cnt;
    }


    void BASE_INDEX_3_CLOSED_HASHTABLE::
    BaseSetSize(int size)
    {
        hash.resize(size);
        for (int i = 0; i < size; i++)
            hash[i].I1() = invalid;
    }

    bool BASE_INDEX_3_CLOSED_HASHTABLE::
    PositionCreate2(const INDEX_3& ind, int& apos)
    {
        int i = HashValue(ind);
        int startpos = i;
        while (1) {
            /*
        i++;
        if (i >= hash.Size()) i = 0;
            */
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
}

