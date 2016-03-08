#ifndef FILE_BitArray_HPP
#define FILE_BitArray_HPP

/**************************************************************************/
/* File:   bitarray.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

#include <climits>
#include <iostream>

#include "array.hpp"
#include "template.hpp"

namespace meshit {

    /**
       data type BitArray
   
       BitArray is a compressed array of Boolean information. By Set and Clear
       the whole array or one bit can be set or reset, respectively. 
       Test returns the state of the accoring bit.
       No range checking is done.

       index ranges from 0 to size-1
     */
    class BitArray
    {
        INDEX size;
        unsigned char* data;

     public:
        BitArray();
        explicit BitArray(INDEX asize);
        ~BitArray();

        void SetSize(INDEX asize);

        INDEX Size() const
        {
            return size;
        }

        void Set();

        void Set(INDEX i)
        {
            data[Addr(i)] |= Mask(i);
        }

        void Clear();

        bool Test(INDEX i) const
        {
            return data[i / CHAR_BIT] & (static_cast<char>(1) << (i % CHAR_BIT));
        }

     private:
        inline unsigned char Mask(INDEX i) const
        {
            return static_cast<char>(1) << (i % CHAR_BIT);
        }

        inline INDEX Addr(INDEX i) const
        {
            return (i / CHAR_BIT);
        }

        BitArray& operator=(BitArray&) { return *this; }

        BitArray(const BitArray&) { }
    };

    // print bitarray

    inline std::ostream& operator<<(std::ostream& s, const BitArray& a)
    {
        for (int i = 1; i <= a.Size(); i++) {
            s << int(a.Test(i));
            if (i % 40 == 0) s << "\n";
        }
        if (a.Size() % 40 != 0) s << "\n";
        return s;
    }

    /**
       data type BitArrayChar
   
       BitArray is an array of Boolean information. By Set and Clear
       the whole array or one bit can be set or reset, respectively. 
       Test returns the state of the accoring bit.
       No range checking is done.
     */
    template<int BASE = 1>
    class BitArrayChar
    {
        Array<char, BASE> data;

     public:
        BitArrayChar() { }

        explicit BitArrayChar(int asize)
                : data(asize) { }

        ~BitArrayChar() { }

        inline int Size() const
        {
            return data.size();
        }

        void Set();

        inline void Set(int i)
        {
            data[i] = 1;
        }

        void Clear();

        inline int Test(int i) const
        {
            return data[i];
        }

     private:
        ///  copy bitarray is not supported

        BitArrayChar& operator=(BitArrayChar&)
        {
            return *this;
        }
        ///  copy bitarray is not supported

        BitArrayChar(const BitArrayChar&) { }
    };

    template<int BASE>
    inline std::ostream& operator<<(std::ostream& s, const BitArrayChar<BASE>& a)
    {
        for (int i = BASE; i < a.Size() + BASE; i++) {
            s << a.Test(i);
            if ((i - BASE) % 40 == 39) s << "\n";
        }
        if (a.Size() % 40 != 0) s << "\n";
        return s;
    }

}  // namespace meshit

#endif
