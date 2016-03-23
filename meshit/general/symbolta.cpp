/**************************************************************************/
/* File:   symbolta.cc                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type Symbol Table
*/

#include "symbolta.hpp"

namespace meshit
{
    BASE_SYMBOLTABLE::~BASE_SYMBOLTABLE()
    {
        DelNames();
    }

    void BASE_SYMBOLTABLE::DelNames()
    {
        for (size_t i = 0; i < names.size(); i++) {
            delete[] names[i];
        }
        names.resize(0);
    }

    size_t BASE_SYMBOLTABLE::Index(const char* name) const
    {
        if (name) {
            for (size_t i = 0; i < names.size(); i++) {
                if (strcmp(names[i], name) == 0) {
                    return i + 1;
                }
            }
        }
        return 0;
    }
}  // namespace meshit
