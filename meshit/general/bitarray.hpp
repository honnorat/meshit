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
       data type BitArrayChar
   
       BitArray is an array of Boolean information. By Set and Clear
       the whole array or one bit can be set or reset, respectively. 
       Test returns the state of the accoring bit.
       No range checking is done.
     */
    class BitArrayChar
    {
        Array<char> data;

     public:
        BitArrayChar() { }

        explicit BitArrayChar(size_t asize)
                : data(asize) { }

        ~BitArrayChar() { }

        inline size_t Size() const
        {
            return data.size();
        }

        inline void Set(size_t i)
        {
            data[i] = 1;
        }

        void Clear()
        {
            data = 0;
        }

        inline int Test(size_t i) const
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

    inline std::ostream& operator<<(std::ostream& s, const BitArrayChar& a)
    {
        for (size_t i = 0; i < a.Size(); i++) {
            s << a.Test(i);
            if (i % 40 == 39) s << "\n";
        }
        if (a.Size() % 40 != 0) s << "\n";
        return s;
    }

}  // namespace meshit

#endif
