#ifndef FILE_ARRAY_HPP
#define FILE_ARRAY_HPP

#include <algorithm>
#include <iostream>
#include <cstring>

namespace meshit {

    /**
       A simple array container.
       Array represented by size and data-pointer.
       No memory allocation and deallocation, must be provided by user.
       Helper functions for printing. 
       Optional range check by macro RANGE_CHECK
     */

    template<typename T>
    class FlatArray
    {
     protected:
        size_t _size;
        T* _data;

     public:
        typedef T TELEM;

        FlatArray(size_t asize, T* adata)
                : _size(asize), _data(adata) { }

        size_t size() const { return _size; }

        size_t begin() const { return 0; }

        T& back() { return _data[_size - 1]; }

        const T& back() const { return _data[_size - 1]; }

        /// Access array.
        T& operator[](size_t i)
        {
            return _data[i];
        }

        const T& operator[](size_t i) const
        {
            return _data[i];
        }

        /// Fill array with value val

        FlatArray& operator=(const T& val)
        {
            for (size_t i = 0; i < _size; i++) {
                _data[i] = val;
            }
            return *this;
        }

        /// takes range starting from position start of end-start elements
        const FlatArray<T> Range(size_t start, size_t end)
        {
            return FlatArray<T>(end - start, _data + start);
        }
    };

}  // namespace meshit

#endif
