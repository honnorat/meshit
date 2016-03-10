/**************************************************************************/
/* File:   table.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type TABLE
*/

#include "table.hpp"

namespace meshit {

    BASE_TABLE::BASE_TABLE(size_t size)
            : data(size)
    {
        for (size_t i = 0; i < size; i++) {
            data[i].maxsize = 0;
            data[i].size = 0;
            data[i].col = NULL;
        }
        oneblock = NULL;
    }

    BASE_TABLE::BASE_TABLE(const FlatArray<int>& entrysizes, size_t elemsize)
            : data(entrysizes.size())
    {
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

    BASE_TABLE::~BASE_TABLE()
    {
        if (oneblock) {
            delete[] oneblock;
        } else {
            for (size_t i = 0; i < data.size(); i++) {
                delete[] (char*) data[i].col;
            }
        }
    }

    void BASE_TABLE::SetSize(size_t size)
    {
        for (size_t i = 0; i < data.size(); i++) {
            delete[] (char*) data[i].col;
        }
        data.resize(size);
        for (size_t i = 0; i < size; i++) {
            data[i].maxsize = 0;
            data[i].size = 0;
            data[i].col = NULL;
        }
    }

    void BASE_TABLE::ChangeSize(size_t size)
    {
        size_t oldsize = data.size();
        if (size == oldsize)
            return;

        if (size < oldsize) {
            for (size_t i = size; i < oldsize; i++) {
                delete[] (char*) data[i].col;
            }
        }

        data.resize(size);

        for (size_t i = oldsize; i < size; i++) {
            data[i].maxsize = 0;
            data[i].size = 0;
            data[i].col = NULL;
        }
    }

    void BASE_TABLE::IncSize2(size_t i, size_t elsize)
    {
#ifdef DEBUG
        if (i < 0 || i >= data.size())
        {
            MyError("BASE_TABLE::Inc: Out of range");
            return;
        }
#endif
        linestruct& line = data[i];
        if (line.size == line.maxsize) {
            void* p = new char[(line.maxsize + 5) * elsize];

            memcpy(p, line.col, line.maxsize * elsize);
            delete[] (char*) line.col;

            line.col = p;
            line.maxsize += 5;
        }

        line.size++;
    }

    void BASE_TABLE::AllocateElementsOneBlock(size_t elemsize)
    {
        size_t cnt = 0;
        size_t n = data.size();

        for (size_t i = 0; i < n; i++) {
            cnt += data[i].maxsize;
        }
        oneblock = new char[elemsize * cnt];

        cnt = 0;
        for (size_t i = 0; i < n; i++) {
            data[i].size = 0;
            data[i].col = &oneblock[elemsize * cnt];
            cnt += data[i].maxsize;
        }
    }

    size_t BASE_TABLE::AllocatedElements() const
    {
        size_t els = 0;
        for (size_t i = 0; i < data.size(); i++) {
            els += data[i].maxsize;
        }
        return els;
    }

    size_t BASE_TABLE::UsedElements() const
    {
        size_t els = 0;
        for (size_t i = 0; i < data.size(); i++) {
            els += data[i].size;
        }
        return els;
    }

}  // namespace meshit
