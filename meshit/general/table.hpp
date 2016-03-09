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
            int size;
            int maxsize;
            void* col;
        };

        Array<linestruct> data;
        char* oneblock;

     public:
        explicit BASE_TABLE(int size);
        BASE_TABLE(const FlatArray<int>& entrysizes, int elemsize);
        ~BASE_TABLE();

        void SetSize(int size);
        void ChangeSize(int size);

        /// increment size of entry i by one, i is 0-based

        void IncSize(int i, int elsize)
        {
            if (data[i].size < data[i].maxsize)
                data[i].size++;
            else
                IncSize2(i, elsize);
        }

        void IncSize2(int i, int elsize);
        void SetEntrySize2(int i, int newsize, int elsize);

        void AllocateElementsOneBlock(int elemsize);

        int AllocatedElements() const;
        int UsedElements() const;
    };

    /** 
       Abstract data type TABLE.
   
       To an integer i in the range from 1 to size a set of elements of the
       generic type T is associated. 
     */
    template<class T, int BASE = 0>
    class TABLE : public BASE_TABLE
    {
     public:
        /// Creates table.
        TABLE() : BASE_TABLE(0) { }

        /// Creates table of size size
        explicit TABLE(int size) : BASE_TABLE(size) { }

        /// Creates fixed maximal element size table
        explicit TABLE(const FlatArray<int, BASE>& entrysizes)
                : BASE_TABLE(FlatArray<int>(entrysizes.size(), const_cast<int*> (&entrysizes[BASE])), sizeof(T))
        {
        }

        /// Changes Size of table to size, deletes data
        inline void SetSize(int size)
        {
            BASE_TABLE::SetSize(size);
        }

        /// Changes Size of table to size, keep data
        inline void ChangeSize(int size)
        {
            BASE_TABLE::ChangeSize(size);
        }

        /// Inserts element acont into row i, BASE-based. Does not test if already used.
        inline void Add(int i, const T& acont)
        {
            IncSize(i - BASE, sizeof(T));
            ((T*) data[i - BASE].col)[data[i - BASE].size - 1] = acont;
        }

        /// Inserts element acont into row i, 1-based. Does not test if already used.
        inline void Add1(int i, const T& acont)
        {
            IncSize(i - 1, sizeof(T));
            ((T*) data[i - 1].col)[data[i - 1].size - 1] = acont;
        }

        void IncSizePrepare(int i)
        {
            data[i - BASE].maxsize++;
        }

        /// Inserts element acont into row i. BASE-based. Does not test if already used, assumes to have enough memory
        inline void AddSave(int i, const T& acont)
        {
            ((T*) data[i - BASE].col)[data[i - BASE].size] = acont;
            data[i - BASE].size++;
        }

        /// Inserts element acont into row i. Does not test if already used.
        inline void AddEmpty(int i)
        {
            IncSize(i - BASE, sizeof(T));
        }

        /** Set the nr-th element in the i-th row to acont.
          Does not check for overflow. */
        inline void Set(int i, int nr, const T& acont)
        {
            ((T*) data[i - 1].col)[nr - 1] = acont;
        }

        /** Returns the nr-th element in the i-th row.
          Does not check for overflow. */
        inline const T& Get(int i, int nr) const
        {
            return ((T*) data[i - 1].col)[nr - 1];
        }

        /// Returns size of the table.
        inline int Size() const
        {
            return data.size();
        }

        /// Returns size of the i-th row.
        inline int EntrySize(int i) const
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
            if (i - BASE < 0 || i - BASE >= data.size())
                std::cout << "table out of range, i = " << i << ", s = " << data.size() << std::endl;
#endif
            return FlatArray<T>(data[i - BASE].size, (T*) data[i - BASE].col);
        }
    };

    template<class T, int BASE>
    inline std::ostream& operator<<(std::ostream& ost, const TABLE<T, BASE>& table)
    {
        for (int i = BASE; i < table.Size() + BASE; i++) {
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
