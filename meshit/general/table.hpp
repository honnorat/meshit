#ifndef FILE_TABLE_HPP
#define FILE_TABLE_HPP

/**************************************************************************/
/* File:   table.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include "array.hpp"

namespace meshit {


    /// Base class to generic class TABLE.

    class BASE_TABLE
    {
     protected:
        class linestruct
        {
         public:
            size_t size;
            size_t maxsize;
            void* col;
        };

        Array<linestruct> data;
        char* oneblock;

     public:
        explicit BASE_TABLE(size_t size);
        BASE_TABLE(const FlatArray<int>& entrysizes, size_t elemsize);
        ~BASE_TABLE();

        void SetSize(size_t size);
        void ChangeSize(size_t size);

        /// increment size of entry i by one, i is 0-based

        void IncSize(size_t i, size_t elsize)
        {
            if (data[i].size < data[i].maxsize) {
                data[i].size++;
            } else {
                IncSize2(i, elsize);
            }
        }

        void IncSize2(size_t i, size_t elsize);

        void AllocateElementsOneBlock(size_t elemsize);

        size_t AllocatedElements() const;
        size_t UsedElements() const;
    };

    /** 
       Abstract data type TABLE.
   
       To an integer i in the range from 1 to size a set of elements of the
       generic type T is associated. 
     */
    template<class T>
    class TABLE : public BASE_TABLE
    {
     public:
        /// Creates table of size size
        explicit TABLE(size_t size = 0) : BASE_TABLE(size) { }

        /// Creates fixed maximal element size table
        explicit TABLE(const FlatArray<int>& entrysizes)
                : BASE_TABLE(FlatArray<int>(entrysizes.size(), const_cast<int*> (&entrysizes[0])), sizeof(T)) { }

        /// Changes Size of table to size, deletes data
        inline void SetSize(size_t size)
        {
            BASE_TABLE::SetSize(size);
        }

        /// Changes Size of table to size, keep data
        inline void ChangeSize(size_t size)
        {
            BASE_TABLE::ChangeSize(size);
        }

        /// Inserts element acont into row i. Does not test if already used.
        inline void Add(size_t i, const T& acont)
        {
            IncSize(i, sizeof(T));
            ((T*) data[i].col)[data[i].size - 1] = acont;
        }

        /// Inserts element acont into row i, 1-based. Does not test if already used.
        inline void Add1(size_t i, const T& acont)
        {
            IncSize(i - 1, sizeof(T));
            ((T*) data[i - 1].col)[data[i - 1].size - 1] = acont;
        }

        void IncSizePrepare(size_t i)
        {
            data[i].maxsize++;
        }

        /// Inserts element acont into row i. Does not test if already used, assumes to have enough memory
        inline void AddSave(size_t i, const T& acont)
        {
            ((T*) data[i].col)[data[i].size] = acont;
            data[i].size++;
        }

        /** Set the nr-th element in the i-th row to acont.
          Does not check for overflow. */
        inline void Set(size_t i, size_t nr, const T& acont)
        {
            ((T*) data[i - 1].col)[nr - 1] = acont;
        }

        /** Returns the nr-th element in the i-th row.
          Does not check for overflow. */
        inline const T& Get(size_t i, size_t nr) const
        {
            return ((T*) data[i - 1].col)[nr - 1];
        }

        /// Returns size of the table.
        inline size_t Size() const
        {
            return data.size();
        }

        /// Returns size of the i-th row.
        inline size_t EntrySize(int i) const
        {
            return data[i - 1].size;
        }

        void AllocateElementsOneBlock()
        {
            BASE_TABLE::AllocateElementsOneBlock(sizeof(T));
        }

        inline void PrintMemInfo(std::ostream& ost) const
        {
            int els = AllocatedElements();
            ost << "table: allocaed " << els
            << " a " << sizeof(T) << " Byts = "
            << els * sizeof(T)
            << " bytes in " << Size() << " bags."
            << " used: " << UsedElements()
            << std::endl;
        }

        /// Access entry.
        FlatArray<T> operator[](int i) const
        {
#ifdef DEBUG
            if (i < 0 || i >= data.size())
                std::cout << "table out of range, i = " << i << ", s = " << data.size() << std::endl;
#endif
            return FlatArray<T>(data[i].size, (T*) data[i].col);
        }
    };

    template<class T>
    inline std::ostream& operator<<(std::ostream& ost, const TABLE<T>& table)
    {
        for (int i = 0; i < table.Size(); i++) {
            ost << i << ": ";
            FlatArray<T> row = table[i];
            ost << "(" << row.size() << ") ";
            for (int j = 0; j < row.size(); j++)
                ost << row[j] << " ";
            ost << std::endl;
        }
        return ost;
    }

}  // namespace meshit

#endif
