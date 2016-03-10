#ifndef FILE_HASHTABL_HPP
#define FILE_HASHTABL_HPP

/**************************************************************************/
/* File:   hashtabl.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <iostream>

#include "table.hpp"
#include "template.hpp"

namespace meshit {

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

        size_t GetNBags() const
        {
            return cont.Size();
        }

        size_t GetBagSize(size_t bnr) const
        {
            return cont.EntrySize(bnr);
        }

        void GetData(size_t bnr, size_t colnr, INDEX_2& ahash, T& acont) const
        {
            ahash = hash.Get(bnr, colnr);
            acont = cont.Get(bnr, colnr);
        }

        void SetData(size_t bnr, size_t colnr, const INDEX_2& ahash, const T& acont)
        {
            hash.Set(bnr - 1, colnr - 1, ahash);
            cont.Set(bnr - 1, colnr - 1, acont);
        }
    };


    class BASE_INDEX_3_HASHTABLE
    {
     protected:
        TABLE<INDEX_3> hash;

     public:
        explicit BASE_INDEX_3_HASHTABLE(int size)
                : hash(size) { }

     protected:
        int HashValue(const INDEX_3& ind) const
        {
            return (ind.I1() + ind.I2() + ind.I3()) % hash.Size();
        }

        int Position(int bnr, const INDEX_3& ind) const
        {
            const INDEX_3* pi = &hash.Get(bnr, 0);
            int n = hash.EntrySize(bnr);
            for (int i = 1; i <= n; ++i, ++pi) {
                if (*pi == ind)
                    return i;
            }
            return 0;
        }
    };

    template<class T>
    class INDEX_3_HASHTABLE : private BASE_INDEX_3_HASHTABLE
    {
        TABLE<T> cont;

     public:
        explicit INDEX_3_HASHTABLE(size_t size)
                : BASE_INDEX_3_HASHTABLE(size), cont(size) { }

        inline void Set(const INDEX_3& ahash, const T& acont);
        inline const T& Get(const INDEX_3& ahash) const;
        inline bool Used(const INDEX_3& ahash) const;
        inline size_t GetNBags() const;
        inline size_t GetBagSize(size_t bnr) const;
        inline void GetData(size_t bnr, size_t colnr, INDEX_3& ahash, T& acont) const;

        class Iterator
        {
         protected:
            const INDEX_3_HASHTABLE& ht;
            int bagnr, pos;

         public:
            Iterator(const INDEX_3_HASHTABLE& aht, int abagnr, int apos)
                    : ht(aht), bagnr(abagnr), pos(apos) { }

            int BagNr() const
            {
                return bagnr;
            }

            int Pos() const
            {
                return pos;
            }

            void operator++(int)
            {
                pos++;
                while (bagnr < ht.GetNBags() &&
                       pos == ht.GetBagSize(bagnr)) {
                    pos = 0;
                    bagnr++;
                }
            }

            bool operator!=(int i) const
            {
                return bagnr != i;
            }
        };

        Iterator Begin() const
        {
            Iterator it(*this, 0, -1);
            it++;
            return it;
        }

        int End() const
        {
            return GetNBags();
        }

        const INDEX_3& GetHash(const Iterator& it) const
        {
            return hash[it.BagNr()][it.Pos()];
        }

        const T& GetData(const Iterator& it) const
        {
            return cont[it.BagNr()][it.Pos()];
        }
    };

    template<typename T>
    inline std::ostream& operator<<(std::ostream& ost, const INDEX_3_HASHTABLE<T>& ht)
    {
        for (typename INDEX_3_HASHTABLE<T>::Iterator it = ht.Begin();
             it != ht.End(); it++) {
            ost << ht.GetHash(it) << ": " << ht.GetData(it) << std::endl;
        }

        return ost;
    }


    class BASE_INDEX_2_CLOSED_HASHTABLE
    {
     protected:
        Array<INDEX_2> hash;
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

        size_t UsedElements() const;

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
        void BaseSetSize(size_t asize);
    };

    template<class T>
    class INDEX_2_CLOSED_HASHTABLE : public BASE_INDEX_2_CLOSED_HASHTABLE
    {
        Array<T> cont;

     public:
        explicit INDEX_2_CLOSED_HASHTABLE(int size)
                : BASE_INDEX_2_CLOSED_HASHTABLE(size), cont(size) { }

        inline void Set(const INDEX_2& ahash, const T& acont);
        inline const T& Get(const INDEX_2& ahash) const;
        inline bool Used(const INDEX_2& ahash) const;
        inline void GetData(int pos, INDEX_2& ahash, T& acont) const;
        inline void SetSize(int size);
        inline void PrintMemInfo(std::ostream& ost) const;
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

    class BASE_INDEX_3_CLOSED_HASHTABLE
    {
     protected:
        Array<INDEX_3> hash;
        INDEX invalid;

     protected:
        explicit BASE_INDEX_3_CLOSED_HASHTABLE(size_t size)
                : hash(size)
        {
            invalid = -1;
            for (size_t i = 0; i < size; i++) {
                hash[i].I1() = invalid;
            }
        }

     public:
        size_t Size() const
        {
            return hash.size();
        }

        bool UsedPos(size_t pos) const
        {
            return hash[pos].I1() != invalid;
        }

        size_t HashValue(const INDEX_3& ind) const
        {
            return (ind.I1() + 15 * ind.I2() + 41 * ind.I3()) % hash.size();
        }

        int Position(const INDEX_3& ind) const
        {
            size_t i = HashValue(ind);
            while (true) {
                if (hash[i] == ind) return i;
                if (hash[i].I1() == invalid) return -1;
                i = (i + 1) % hash.size();
            }
        }

        bool PositionCreate(const INDEX_3& ind, size_t& apos)
        {
            size_t i = HashValue(ind);
            if (hash[i] == ind) {
                apos = i;
                return false;
            }
            if (hash[i].I1() == invalid) {
                hash[i] = ind;
                apos = i;
                return true;
            }
            return PositionCreate2(ind, apos);
        }

     protected:
        bool PositionCreate2(const INDEX_3& ind, size_t& apos);
        void BaseSetSize(size_t asize);
    };

    template<class T>
    class INDEX_3_CLOSED_HASHTABLE : public BASE_INDEX_3_CLOSED_HASHTABLE
    {
        Array<T> cont;

     public:
        explicit INDEX_3_CLOSED_HASHTABLE(int size)
                : BASE_INDEX_3_CLOSED_HASHTABLE(size), cont(size) { }

        void Set(const INDEX_3& ahash, const T& acont)
        {
            size_t pos;
            PositionCreate(ahash, pos);
            hash[pos] = ahash;
            cont[pos] = acont;
        }

        bool Used(const INDEX_3& ahash) const
        {
            return (Position(ahash) != -1);
        }

        void GetData(int pos, INDEX_3& ahash, T& acont) const
        {
            ahash = hash[pos];
            acont = cont[pos];
        }

        void SetSize(size_t size)
        {
            BaseSetSize(size);
            cont.resize(size);
        }

        void PrintMemInfo(std::ostream& ost) const
        {
            std::cout << "Hashtable: " << Size()
            << " entries of size " << sizeof(INDEX_3) << " + " << sizeof(T)
            << " = " << Size() * (sizeof(INDEX_3) + sizeof(T)) << " bytes" << std::endl;
        }
    };

    template<typename T>
    inline std::ostream& operator<<(std::ostream& ost, const INDEX_3_CLOSED_HASHTABLE<T>& ht)
    {
        for (int i = 0; i < ht.Size(); i++) {
            if (ht.UsedPos(i)) {
                INDEX_3 hash;
                T data;
                ht.GetData(i, hash, data);
                ost << "hash = " << hash << ", data = " << data << std::endl;
            }
        }
        return ost;
    }

    template<class T>
    inline void INDEX_3_HASHTABLE<T>::Set(const INDEX_3& ahash, const T& acont)
    {
        int bnr = HashValue(ahash);
        int pos = Position(bnr, ahash);
        if (pos > 0) {
            cont.Set(bnr, pos - 1, acont);
        } else {
            hash.Add(bnr, ahash);
            cont.Add(bnr, acont);
        }
    }

    template<class T>
    inline const T& INDEX_3_HASHTABLE<T>::Get(const INDEX_3& ahash) const
    {
        int bnr = HashValue(ahash);
        int pos = Position(bnr, ahash);
        return cont.Get(bnr, pos - 1);
    }

    template<class T>
    inline bool INDEX_3_HASHTABLE<T>::Used(const INDEX_3& ahash) const
    {
        return (Position(HashValue(ahash), ahash) > 0);
    }

    template<class T>
    inline size_t INDEX_3_HASHTABLE<T>::GetNBags() const
    {
        return cont.Size();
    }

    template<class T>
    inline size_t INDEX_3_HASHTABLE<T>::GetBagSize(size_t bnr) const
    {
        return cont.EntrySize(bnr);
    }

    template<class T>
    inline void INDEX_3_HASHTABLE<T>::GetData(size_t bnr, size_t colnr, INDEX_3& ahash, T& acont) const
    {
        ahash = hash.Get(bnr, colnr);
        acont = cont.Get(bnr, colnr);
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

    template<class T>
    inline void INDEX_2_CLOSED_HASHTABLE<T>::
    SetSize(int size)
    {
        BaseSetSize(size);
        cont.resize(size);
    }

    template<class T>
    inline void INDEX_2_CLOSED_HASHTABLE<T>::
    PrintMemInfo(std::ostream& ost) const
    {
        std::cout << "Hashtable: " << Size()
        << " entries of size " << sizeof(INDEX_2) << " + " << sizeof(T)
        << " = " << Size() * (sizeof(INDEX_2) + sizeof(T)) << " bytes."
        << " Used els: " << UsedElements()
        << std::endl;
    }

}  // namespace meshit


#endif
