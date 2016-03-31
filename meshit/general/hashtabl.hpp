#ifndef FILE_HASHTABL_HPP
#define FILE_HASHTABL_HPP

/**************************************************************************/
/* File:   hashtabl.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <iostream>

#include "table.hpp"
#include "index.hpp"

namespace meshit
{
    class BASE_INDEX_2_HASHTABLE
    {
     protected:
        TABLE<INDEX_2> hash;

     public:
        explicit BASE_INDEX_2_HASHTABLE(size_t size)
            : hash(size) { }

        int HashValue(const INDEX_2& ind) const
        {
            return (ind.I1() + ind.I2()) % hash.Size();
        }

        size_t Position(int bnr, const INDEX_2& ind) const
        {
            for (size_t i = 0; i < hash.EntrySize(bnr); i++) {
                if (hash.Get(bnr, i) == ind) {
                    return i + 1;
                }
            }
            return 0;
        }
    };

    template<class T>
    class INDEX_2_HASHTABLE : public BASE_INDEX_2_HASHTABLE
    {
        TABLE<T> cont;

     public:
        explicit INDEX_2_HASHTABLE(size_t size)
            : BASE_INDEX_2_HASHTABLE(size), cont(size) { }

        void Set(const INDEX_2& ahash, const T& acont)
        {
            int bnr = HashValue(ahash);
            int pos = Position(bnr, ahash);
            if (pos) {
                cont.Set(bnr, pos - 1, acont);
            } else {
                hash.Add(bnr, ahash);
                cont.Add(bnr, acont);
            }
        }

        const T& Get(const INDEX_2& ahash) const
        {
            int bnr = HashValue(ahash);
            int pos = Position(bnr, ahash);
            return cont.Get(bnr, pos - 1);
        }

        bool Used(const INDEX_2& ahash) const
        {
            return Position(HashValue(ahash), ahash) > 0;
        }
    };


    class BASE_INDEX_2_CLOSED_HASHTABLE
    {
     protected:
        std::vector<INDEX_2> hash;
        INDEX invalid;

     public:
        explicit BASE_INDEX_2_CLOSED_HASHTABLE(size_t size);

        size_t Size() const
        {
            return hash.size();
        }

        size_t UsedPos(size_t pos) const
        {
            return hash[pos - 1].I1() != invalid;
        }

        size_t HashValue(const INDEX_2& ind) const
        {
            return (ind.I1() + 71 * ind.I2()) % hash.size() + 1;
        }

        size_t Position(const INDEX_2& ind) const
        {
            size_t i = HashValue(ind);
            while (true) {
                if (hash[i - 1] == ind) return i;
                if (hash[i - 1].I1() == invalid) return 0;
                i++;
                if (i > hash.size()) i = 1;
            }
        }

        // returns 1, if new postion is created
        bool PositionCreate(const INDEX_2& ind, size_t& apos)
        {
            size_t i = HashValue(ind);
            if (hash[i - 1] == ind) {
                apos = i;
                return false;
            }
            if (hash[i - 1].I1() == invalid) {
                hash[i - 1] = ind;
                apos = i;
                return true;
            }
            return PositionCreate2(ind, apos);
        }

     protected:
        bool PositionCreate2(const INDEX_2& ind, size_t& apos);
    };

    template<class T>
    class INDEX_2_CLOSED_HASHTABLE : public BASE_INDEX_2_CLOSED_HASHTABLE
    {
        std::vector<T> cont;

     public:
        explicit INDEX_2_CLOSED_HASHTABLE(int size)
            : BASE_INDEX_2_CLOSED_HASHTABLE(size), cont(size) { }

        inline void Set(const INDEX_2& ahash, const T& acont);
        inline const T& Get(const INDEX_2& ahash) const;
        inline bool Used(const INDEX_2& ahash) const;
        inline void GetData(int pos, INDEX_2& ahash, T& acont) const;
    };

    template<typename T>
    inline std::ostream& operator<<(std::ostream& ost, const INDEX_2_CLOSED_HASHTABLE<T>& ht)
    {
        for (int i = 0; i < ht.Size(); i++)
            if (ht.UsedPos(i)) {
                INDEX_2 hash;
                T data;
                ht.GetData(i, hash, data);
                ost << "hash = " << hash << ", data = " << data << std::endl;
            }
        return ost;
    }

    /* *********** Closed Hashing ************************* */

    template<class T>
    inline void INDEX_2_CLOSED_HASHTABLE<T>::
    Set(const INDEX_2& ahash, const T& acont)
    {
        size_t pos;
        PositionCreate(ahash, pos);
        hash[pos - 1] = ahash;
        cont[pos - 1] = acont;
    }

    template<class T>
    inline const T& INDEX_2_CLOSED_HASHTABLE<T>::
    Get(const INDEX_2& ahash) const
    {
        int pos = Position(ahash);
        return cont[pos - 1];
    }

    template<class T>
    inline bool INDEX_2_CLOSED_HASHTABLE<T>::
    Used(const INDEX_2& ahash) const
    {
        int pos = Position(ahash);
        return (pos != 0);
    }

    template<class T>
    inline void INDEX_2_CLOSED_HASHTABLE<T>::
    GetData(int pos, INDEX_2& ahash, T& acont) const
    {
        ahash = hash[pos - 1];
        acont = cont[pos - 1];
    }

}  // namespace meshit


#endif
