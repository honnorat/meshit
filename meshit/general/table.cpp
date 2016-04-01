/**************************************************************************/
/* File:   table.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type TABLE
*/

#include "table.hpp"

#include <cstring>

namespace meshit
{
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

    BASE_TABLE::BASE_TABLE(const std::vector<int>& entrysizes, size_t elemsize)
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
                delete[] reinterpret_cast<char*>(data[i].col);
            }
        }
    }

    void BASE_TABLE::SetSize(size_t size)
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
            size_t new_size = 2 * std::max(static_cast<size_t>(1), line.maxsize);
            void* new_col = new char[new_size * elsize];

            memcpy(new_col, line.col, line.maxsize * elsize);
            delete[] reinterpret_cast<char*>(line.col);

            line.col = new_col;
            line.maxsize = new_size;
        }

        line.size++;
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
