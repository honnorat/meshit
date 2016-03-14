#ifndef FILE_Array
#define FILE_Array

#include <algorithm>
#include <iostream>
#include <cstring>

namespace meshit {

    template<class TA1, class TA2>
    class IndirectArray;

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

        size_t size() const
        {
            return _size;
        }

        /// Access array.
        T& operator[](size_t i)
        {
            return _data[i];
        }

        const T& operator[](size_t i) const
        {
            return _data[i];
        }

        template<typename T2>
        IndirectArray<FlatArray, FlatArray<T2> > operator[](const FlatArray<T2>& ia) const
        {
            return IndirectArray<FlatArray, FlatArray<T2> >(*this, ia);
        }

        /// access last element. check by macro CHECK_RANGE

        T& Last() const
        {
            return _data[_size - 1];
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

    /** 
        Dynamic array container.
   
        Array<T> is an automatically increasing array container.
        The allocated memory doubles on overflow. 
        Either the container takes care of memory allocation and deallocation,
        or the user provides one block of data.
     */
    template<class T>
    class Array : public FlatArray<T>
    {
     protected:
        using FlatArray<T>::_size;
        using FlatArray<T>::_data;

        size_t _allocsize;  // physical size of array
        bool _ownmem;       // memory is responsibility of container

     public:
        Array() : FlatArray<T>(0, nullptr)
        {
            _allocsize = 0;
            _ownmem = 1;
        }

        explicit Array(size_t asize)
                : FlatArray<T>(asize, new T[asize])
        {
            _allocsize = asize;
            _ownmem = 1;
        }

        /// Generate array in user data

        Array(size_t asize, T* adata)
                : FlatArray<T>(asize, adata)
        {
            _allocsize = asize;
            _ownmem = 0;
        }

        /// array copy

        explicit Array(const Array<T>& a2)
                : FlatArray<T>(a2.size(), (a2.size() ? new T[a2.size()] : 0))
        {
            _allocsize = _size;
            _ownmem = 1;
            for (size_t i = 0; i < _size; i++) {
                (*this)[i] = a2[i];
            }
        }

        ~Array()
        {
            if (_ownmem)
                delete[] _data;
        }

        /// Change logical size. If necessary, do reallocation. Keeps contents.

        void resize(size_t nsize)
        {
            if (nsize > _allocsize)
                ReSize(nsize);
            _size = nsize;
        }

        /// Change physical size. Keeps logical size. Keeps contents.

        void reserve(size_t nallocsize)
        {
            if (nallocsize > _allocsize)
                ReSize(nallocsize);
        }


        /// Add element at end of array. reallocation if necessary.
        size_t push_back(const T& el)
        {
            if (_size == _allocsize)
                ReSize(_size + 1);
            _data[_size] = el;
            _size++;
            return _size;
        }

        template<typename T2>
        void Append(FlatArray<T2>& a2)
        {
            if (_size + a2.size() > _allocsize)
                ReSize(_size + a2.size());
            for (size_t i = 0; i < a2.size(); i++)
                _data[_size + i] = a2[i];
            _size += a2.size();
        }

        void
        assign(size_t nsize, const T& el)
        {
            if (_data && _ownmem) {
                delete[] _data;
            }
            _ownmem = 1;
            _data = new T[nsize];
            _size = _allocsize = nsize;

            for (size_t i = 0; i < _size; i++) {
                _data[i] = el;
            }
        }

        /// Delete element i (0-based). Move last element to position i.
        void Delete(size_t i)
        {
            _data[i] = _data[_size - 1];
            _size--;
        }

        /// Delete last element.
        void pop_back()
        {
            _size--;
        }

        /// Deallocate memory
        void DeleteAll()
        {
            if (_ownmem)
                delete[] _data;
            _data = 0;
            _size = _allocsize = 0;
        }

        /// Fill array with val
        Array& operator=(const T& val)
        {
            FlatArray<T>::operator=(val);
            return *this;
        }

        /// array copy
        Array& operator=(const Array& a2)
        {
            resize(a2.size());
            for (size_t i = 0; i < _size; i++)
                (*this)[i] = a2[i];
            return *this;
        }

        /// array copy
        Array& operator=(const FlatArray<T>& a2)
        {
            resize(a2.size());
            for (size_t i = 0; i < _size; i++)
                (*this)[i] = a2[i];
            return *this;
        }

     private:
        void ReSize(size_t minsize)
        {
            // resize array, at least to size minsize. copy contents
            size_t nsize = 2 * _allocsize;
            if (nsize < minsize) nsize = minsize;

            if (_data) {
                T* p = new T[nsize];

                memcpy(p, _data, std::min(nsize, _size) * sizeof(T));

                if (_ownmem) {
                    delete[] _data;
                }
                _ownmem = 1;
                _data = p;
            } else {
                _data = new T[nsize];
                _ownmem = 1;
            }
            _allocsize = nsize;
        }
    };

    template<class T, size_t S>
    class ArrayMem : public Array<T>
    {
        using Array<T>::_size;
        using Array<T>::_data;
        using Array<T>::_ownmem;

        // T mem[S];     // Intel C++ calls dummy constructor
        // char mem[S*sizeof(T)];
        double mem[(S * sizeof(T) + 7) / 8];

     public:
        /// Generate array of logical and physical size asize
        explicit ArrayMem(size_t asize = 0)
                : Array<T>(S, reinterpret_cast<T*>(mem))
        {
            _size = asize;
            if (asize > S) {
                _data = new T[asize];
                _ownmem = 1;
            }
        }

        ArrayMem& operator=(const T& val)
        {
            Array<T>::operator=(val);
            return *this;
        }

        /// array copy
        ArrayMem& operator=(const FlatArray<T>& a2)
        {
            this->resize(a2.size());
            for (size_t i = 0; i < _size; i++)
                (*this)[i] = a2[i];
            return *this;
        }
    };

    template<class TA1, class TA2>
    class IndirectArray
    {
        const TA1& array;
        const TA2& ia;

     public:
        IndirectArray(const TA1& aa, const TA2& aia)
                : array(aa), ia(aia) { }

        size_t size() const
        {
            return ia.size();
        }

        const typename TA1::TELEM& operator[](size_t i) const
        {
            return array[ia[i]];
        }
    };

}  // namespace meshit

#endif
