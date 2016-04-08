#ifndef FILE_TABLE_HPP
#define FILE_TABLE_HPP

/**************************************************************************/
/* File:   table.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <iostream>
#include <vector>

#include <cstring>

namespace meshit
{
    /**
       Abstract data type TABLE.
   
       To an integer i in the range from 1 to size a set of elements of the
       generic type T is associated. 
     */
    template<class T>
    class TABLE
    {
     public:
        /// Creates table of size size
        explicit TABLE(size_t size = 0);
        explicit TABLE(const std::vector<int>& entrysizes);
        ~TABLE();

        /// Returns size of the table.
        inline size_t Size() const
        {
            return data.size();
        }

        void SetSize(size_t size);

        /// increment size of entry i by one, i is 0-based
        void IncSize(size_t i);

        /// Inserts element acont into row i. Does not test if already used.
        inline void Add(size_t i, const T& acont)
        {
            IncSize(i);
            static_cast<T*>(data[i].col)[data[i].size - 1] = acont;
        }

        /// Access entry.
        std::vector<T> operator[](size_t i) const
        {
            std::vector<T> vec(data[i].size);
            for (size_t j = 0; j < vec.size(); j++) {
                vec[j] = static_cast<T*>(data[i].col)[j];
            }
            return vec;
        }

     protected:
        class linestruct
        {
         public:
            size_t size;
            size_t maxsize;
            void* col;
        };

        std::vector<linestruct> data;
        char* oneblock;

    };

    template<class T>
    TABLE<T>::TABLE(size_t size)
        : data(size)
    {
        for (size_t i = 0; i < size; i++) {
            data[i].maxsize = 0;
            data[i].size = 0;
            data[i].col = NULL;
        }
        oneblock = NULL;
    }

    template<class T>
    TABLE<T>::TABLE(const std::vector<int>& entrysizes)
        : data(entrysizes.size())
    {
        constexpr size_t elemsize = sizeof(T);
        size_t cnt = 0;
        size_t n = entrysizes.size();

        for (size_t i = 0; i < n; i++) {
            cnt += entrysizes[i];
        }
        oneblock = new char[elemsize * cnt];

        cnt = 0;
        for (size_t i = 0; i < n; i++) {
            data[i].maxsize = entrysizes[i];
            data[i].size = 0;
            data[i].col = &oneblock[elemsize * cnt];
            cnt += entrysizes[i];
        }
    }

    template<class T>
    TABLE<T>::~TABLE()
    {
        if (oneblock) {
            delete[] oneblock;
        } else {
            for (size_t i = 0; i < data.size(); i++) {
                delete[] reinterpret_cast<char*>(data[i].col);
            }
        }
    }

    template<class T>
    void TABLE<T>::SetSize(size_t size)
    {
        for (size_t i = 0; i < data.size(); i++) {
            delete[] reinterpret_cast<char*>(data[i].col);
        }
        data.resize(size);
        for (size_t i = 0; i < size; i++) {
            data[i].maxsize = 0;
            data[i].size = 0;
            data[i].col = NULL;
        }
    }

    template<class T>
    void TABLE<T>::IncSize(size_t i)
    {
#ifdef DEBUG
        if (i >= data.size())
        {
            MyError("BASE_TABLE::Inc: Out of range");
            return;
        }
#endif
        if (data[i].size < data[i].maxsize) {
            data[i].size++;
            return;
        }

        constexpr size_t elsize = sizeof(T);
        linestruct& line = data[i];
        if (line.size == line.maxsize) {
            size_t new_size = 2 * std::max(static_cast<size_t>(1), line.maxsize);
            void* new_col = new char[new_size * elsize];

            memcpy(new_col, line.col, line.maxsize * elsize);
            delete[] reinterpret_cast<char*>(line.col);

            line.col = new_col;
            line.maxsize = new_size;
        }

        line.size++;
    }

}  // namespace meshit

#endif
