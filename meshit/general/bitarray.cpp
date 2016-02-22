/**************************************************************************/
/* File:   bitarray.cc                                                    */
/* Autho: Joachim Schoeberl                                               */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   data type BitArray
*/

#include "bitarray.hpp"

namespace meshit {

    BitArray::BitArray()
    {
        size = 0;
        data = NULL;
    }

    BitArray::BitArray(int asize)
    {
        size = 0;
        data = NULL;
        SetSize(asize);
    }

    BitArray::~BitArray()
    {
        delete[] data;
    }

    void BitArray::SetSize(int asize)
    {
        if (size == asize) return;
        delete[] data;

        size = asize;
        data = new unsigned char[Addr(size) + 1];
    }

    void BitArray::Set()
    {
        if (!size) return;
        for (int i = 0; i <= Addr(size); i++)
            data[i] = UCHAR_MAX;
    }

    void BitArray::Clear()
    {
        if (!size) return;
        for (int i = 0; i <= Addr(size); i++)
            data[i] = 0;
    }

    template<int BASE>
    void BitArrayChar<BASE>::Set()
    {
        data = 1;
    }

    template<int BASE>
    void BitArrayChar<BASE>::Clear()
    {
        data = 0;
    }

    template
    class BitArrayChar<0>;

    template
    class BitArrayChar<1>;

}  // namespace meshit
