#ifndef FILE_Array
#define FILE_Array

#include <iostream>
#include <cstring>

namespace netgen {

    template <class TA1, class TA2> class IndirectArray;

    /**
       A simple array container.
       Array represented by size and data-pointer.
       No memory allocation and deallocation, must be provided by user.
       Helper functions for printing. 
       Optional range check by macro RANGE_CHECK
     */

    template <typename T, int BASE = 0, typename TIND = int>
    class FlatArray
    {
      protected:

        size_t _size;
        T * _data;

      public:

        typedef T TELEM;

        FlatArray(int asize, T * adata)
            : _size(asize), _data(adata) { }

        int size() const
        {
            return _size;
        }

        TIND Begin() const
        {
            return TIND(BASE);
        }

        TIND End() const
        {
            return TIND(_size + BASE);
        }

        /// Access array. BASE-based

        T & operator[](TIND i) const
        {
            return _data[i - BASE];
        }

        template <typename T2, int B2>
        IndirectArray<FlatArray, FlatArray<T2, B2> > operator[](const FlatArray<T2, B2> & ia) const
        {
            return IndirectArray<FlatArray, FlatArray<T2, B2> > (*this, ia);
        }

        /// Access array, one-based  (old fashioned)

        T & Elem(int i)
        {
            return ((T*) _data)[i - 1];
        }

        /// Access array, one-based  (old fashioned)

        const T & Get(int i) const
        {
            return ((const T*) _data)[i - 1];
        }

        /// Access array, one-based  (old fashioned)

        void Set(int i, const T & el)
        {
            ((T*) _data)[i - 1] = el;
        }

        /// access first element

        T & First() const
        {
            return _data[0];
        }

        /// access last element. check by macro CHECK_RANGE

        T & Last() const
        {
            return _data[_size - 1];
        }

        /// Fill array with value val

        FlatArray & operator=(const T & val)
        {
            for (size_t i = 0; i < _size; i++)
                _data[i] = val;
            return *this;
        }

        /// takes range starting from position start of end-start elements

        const FlatArray<T> Range(TIND start, TIND end)
        {
            return FlatArray<T> (end - start, _data + start);
        }

        /// first position of element elem, returns -1 if element not contained in array 

        TIND Pos(const T & elem) const
        {
            TIND pos = -1;
            for (TIND i = 0; pos == -1 && i < (TIND) this->_size; i++)
                if (elem == _data[i]) pos = i;
            return pos;
        }

        /// does the array contain element elem ?

        bool Contains(const T & elem) const
        {
            return ( Pos(elem) >= 0 );
        }
    };

    // print array

    template <typename T, int BASE, typename TIND>
    inline std::ostream & operator<<(std::ostream & s, const FlatArray<T, BASE, TIND> & a)
    {
        for (TIND i = a.Begin(); i < a.End(); i++)
            s << i << ": " << a[i] << std::endl;
        return s;
    }

    /** 
        Dynamic array container.
   
        Array<T> is an automatically increasing array container.
        The allocated memory doubles on overflow. 
        Either the container takes care of memory allocation and deallocation,
        or the user provides one block of data.
     */
    template <class T, int BASE = 0, typename TIND = int>
    class Array : public FlatArray<T, BASE, TIND>
    {
      protected:
        using FlatArray<T, BASE, TIND>::_size;
        using FlatArray<T, BASE, TIND>::_data;

        size_t _allocsize;  // physical size of array
        bool _ownmem;       // memory is responsibility of container

      public:

        /// Generate array of logical and physical size asize

        explicit Array()
            : FlatArray<T, BASE, TIND> (0, NULL)
        {
            _allocsize = 0;
            _ownmem = 1;
        }

        explicit Array(int asize)
            : FlatArray<T, BASE, TIND> (asize, new T[asize])
        {
            _allocsize = asize;
            _ownmem = 1;
        }

        /// Generate array in user data

        Array(TIND asize, T* adata)
            : FlatArray<T, BASE, TIND> (asize, adata)
        {
            _allocsize = asize;
            _ownmem = 0;
        }

        /// array copy 

        explicit Array(const Array<T, BASE, TIND> & a2)
            : FlatArray<T, BASE, TIND> (a2.size(), a2.size() ? new T[a2.size()] : 0)
        {
            _allocsize = _size;
            _ownmem = 1;
            for (TIND i = BASE; i < (TIND) _size + BASE; i++)
                (*this)[i] = a2[i];
        }

        ~Array()
        {
            if (_ownmem)
                delete [] _data;
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

        int push_back(const T & el)
        {
            if (_size == _allocsize)
                ReSize(_size + 1);
            _data[_size] = el;
            _size++;
            return _size;
        }

        template <typename T2, int B2>
        void Append(FlatArray<T2, B2> a2)
        {
            if (_size + a2.size() > _allocsize)
                ReSize(_size + a2.size());
            for (int i = 0; i < a2.size(); i++)
                _data[_size + i] = a2[i + B2];
            _size += a2.size();
        }

        void
        assign(size_t nsize, const T & el)
        {
            if (_data && _ownmem) {
                delete [] _data;
            }
            _ownmem = 1;
            _data = new T[nsize];
            _size = _allocsize = nsize;

            for (size_t i = 0; i < _size; i++) {
                _data[i] = el;
            }
        }

        /// Delete element i (0-based). Move last element to position i.

        void Delete(TIND i)
        {
            _data[i] = _data[_size - 1];
            _size--;
        }

        /// Delete element i (1-based). Move last element to position i.

        void DeleteElement(TIND i)
        {
            _data[i - 1] = _data[_size - 1];
            _size--;
        }

        /// Delete last element. 

        void DeleteLast()
        {
            _size--;
        }

        /// Deallocate memory

        void DeleteAll()
        {
            if (_ownmem)
                delete [] _data;
            _data = 0;
            _size = _allocsize = 0;
        }

        /// Fill array with val

        Array & operator=(const T & val)
        {
            FlatArray<T, BASE, TIND>::operator=(val);
            return *this;
        }

        /// array copy

        Array & operator=(const Array & a2)
        {
            resize(a2.size());
            for (size_t i = BASE; i < _size + BASE; i++)
                (*this)[i] = a2[i];
            return *this;
        }

        /// array copy

        Array & operator=(const FlatArray<T> & a2)
        {
            resize(a2.size());
            for (size_t i = BASE; i < _size + BASE; i++)
                (*this)[i] = a2[i];
            return *this;
        }


      private:

        /// resize array, at least to size minsize. copy contents

        void ReSize(size_t minsize)
        {
            size_t nsize = 2 * _allocsize;
            if (nsize < minsize) nsize = minsize;

            if (_data) {
                T * p = new T[nsize];

                memcpy(p, _data, std::min(nsize, _size) * sizeof (T));

                if (_ownmem)
                    delete [] _data;
                _ownmem = 1;
                _data = p;
            }
            else {
                _data = new T[nsize];
                _ownmem = 1;
            }
            _allocsize = nsize;
        }
    };

    template <class T, int S>
    class ArrayMem : public Array<T>
    {
        using Array<T>::_size;
        using Array<T>::_data;
        using Array<T>::_ownmem;

        // T mem[S];     // Intel C++ calls dummy constructor
        // char mem[S*sizeof(T)];
        double mem[(S * sizeof (T) + 7) / 8];
      public:
        /// Generate array of logical and physical size asize

        explicit ArrayMem(int asize = 0)
            : Array<T> (S, static_cast<T*> (static_cast<void*> (&mem[0])))
        {
            _size = asize;
            if (asize > S) {
                _data = new T[asize];
                _ownmem = 1;
            }
            // SetSize (asize);
        }

        ArrayMem & operator=(const T & val)
        {
            Array<T>::operator=(val);
            return *this;
        }

        /// array copy

        ArrayMem & operator=(const FlatArray<T> & a2)
        {
            this->resize(a2.size());
            for (size_t i = 0; i < _size; i++)
                (*this)[i] = a2[i];
            return *this;
        }

    };

    /*
    template <class T, int B1, int B2>
    class IndirectArray
    {
      const FlatArray<T, B1> & array;
      const FlatArray<int, B2> & ia; 
    
    public:
      IndirectArray (const FlatArray<T,B1> & aa, const FlatArray<int, B2> & aia)
      : array(aa), ia(aia) { ; }
      int Size() const { return ia.Size(); }
      const T & operator[] (int i) const { return array[ia[i]]; }
    };
     */

    template <class TA1, class TA2>
    class IndirectArray
    {
        const TA1 & array;
        const TA2 & ia;

      public:

        IndirectArray(const TA1 & aa, const TA2 & aia)
            : array(aa), ia(aia)
        {
            ;
        }

        int Size() const
        {
            return ia.Size();
        }

        int Begin() const
        {
            return ia.Begin();
        }

        int End() const
        {
            return ia.End();
        }

        const typename TA1::TELEM & operator[](int i) const
        {
            return array[ia[i]];
        }
    };

    template <typename T1, typename T2>
    inline std::ostream & operator<<(std::ostream & s, const IndirectArray<T1, T2> & ia)
    {
        for (int i = ia.Begin(); i < ia.End(); i++)
            s << i << ": " << ia[i] << std::endl;
        return s;
    }



    /*

    ///
    template <class T, int BASE = 0> 
    class MoveableArray 
    {
      int size;
      int allocsize;
      DynamicMem<T> data;

    public:

      MoveableArray()
      { 
        size = allocsize = 0; 
        data.SetName ("MoveableArray");
      }

      MoveableArray(int asize)
        : size(asize), allocsize(asize), data(asize)
      { ; }
  
      ~MoveableArray () { ; }

      int Size() const { return size; }

      void SetSize(int nsize)
      {
        if (nsize > allocsize) 
      {
        data.ReAlloc (nsize);
        allocsize = nsize;
      }
        size = nsize;
      }

      void SetAllocSize (int nallocsize)
      {
        data.ReAlloc (nallocsize);
        allocsize = nallocsize;
      }

      ///
      T & operator[] (int i)
      { return ((T*)data)[i-BASE]; }

      ///
      const T & operator[] (int i) const
      { return ((const T*)data)[i-BASE]; }

      ///
      T & Elem (int i)
      { return ((T*)data)[i-1]; }
  
      ///
      const T & Get (int i) const 
      { return ((const T*)data)[i-1]; }

      ///
      void Set (int i, const T & el)
      { ((T*)data)[i-1] = el; }

      ///
      T & Last ()
      { return ((T*)data)[size-1]; }
  
      ///
      const T & Last () const
      { return ((const T*)data)[size-1]; }
  
      ///
      int Append (const T & el)
      {
        if (size == allocsize) 
      {
        SetAllocSize (2*allocsize+1);
      }
        ((T*)data)[size] = el;
        size++;
        return size;
      }
  
      ///
      void Delete (int i)
      {
        DeleteElement (i+1);
      }

      ///
      void DeleteElement (int i)
      {
        ((T*)data)[i-1] = ((T*)data)[size-1];
        size--;
      }
  
      ///
      void DeleteLast ()
      { size--; }

      ///
      void DeleteAll ()
      {
        size = allocsize = 0;
        data.Free();
      }

      ///
      void PrintMemInfo (ostream & ost) const
      {
        ost << Size() << " elements of size " << sizeof(T) << " = " 
        << Size() * sizeof(T) <<std::endl;
      }

      MoveableArray & operator= (const T & el)
      {
        for (int i = 0; i < size; i++)
      ((T*)data)[i] = el;
        return *this;
      }


      MoveableArray & Copy (const MoveableArray & a2)
      {
        SetSize (a2.Size());
        for (int i = 0; i < this->size; i++)
      data[i] = a2.data[i];
        return *this;
      }

      /// array copy
      MoveableArray & operator= (const MoveableArray & a2)
      {
        return Copy(a2);
      }


      void SetName (const char * aname)
      {
        data.SetName(aname);
      }
    private:
      ///
      //MoveableArray & operator= (MoveableArray &); //???
      ///
      //MoveableArray (const MoveableArray &); //???
    };


    template <class T>
    inline ostream & operator<< (ostream & ost, MoveableArray<T> & a)
    {
      for (int i = 0; i < a.Size(); i++)
        ost << i << ": " << a[i] <<std::endl;
      return ost;
    }
     */


    /// bubble sort array

    template <class T>
    inline void BubbleSort(const FlatArray<T> & data)
    {
        for (int i = 0; i < data.size(); i++)
            for (int j = i + 1; j < data.size(); j++)
                if (data[i] > data[j]) {
                    T hv = data[i];
                    data[i] = data[j];
                    data[j] = hv;
                }
    }

    /// bubble sort array

    template <class T, class S>
    inline void BubbleSort(FlatArray<T> & data, FlatArray<S> & slave)
    {
        for (int i = 0; i < data.size(); i++)
            for (int j = i + 1; j < data.size(); j++)
                if (data[i] > data[j]) {
                    T hv = data[i];
                    data[i] = data[j];
                    data[j] = hv;

                    S hvs = slave[i];
                    slave[i] = slave[j];
                    slave[j] = hvs;
                }
    }

    template <class T, class S>
    void QuickSortRec(FlatArray<T> & data,
            FlatArray<S> & slave,
            int left, int right)
    {
        int i = left;
        int j = right;
        T midval = data[(left + right) / 2];

        do {
            while (data[i] < midval) i++;
            while (midval < data[j]) j--;

            if (i <= j) {
                std::swap(data[i], data[j]);
                std::swap(slave[i], slave[j]);
                i++;
                j--;
            }
        } while (i <= j);
        if (left < j) QuickSortRec(data, slave, left, j);
        if (i < right) QuickSortRec(data, slave, i, right);
    }

    template <class T, class S>
    void QuickSort(FlatArray<T> & data, FlatArray<S> & slave)
    {
        QuickSortRec(data, slave, 0, data.size() - 1);
    }

    template <class T>
    void Intersection(const FlatArray<T> & in1, const FlatArray<T> & in2,
            Array<T> & out)
    {
        out.resize(0);
        for (int i = 0; i < in1.size(); i++)
            if (in2.Contains(in1[i]))
                out.push_back(in1[i]);
    }

    template <class T>
    void Intersection(const FlatArray<T> & in1, const FlatArray<T> & in2, const FlatArray<T> & in3,
            Array<T> & out)
    {
        out.resize(0);
        for (int i = 0; i < in1.size(); i++)
            if (in2.Contains(in1[i]) && in3.Contains(in1[i]))
                out.push_back(in1[i]);
    }


}

#endif

