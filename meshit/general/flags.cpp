/**************************************************************************/
/* File:   flags.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

/* 
   Datatype Flags
*/

#include "flags.hpp"

#include <sstream>
#include <fstream>

namespace meshit
{
    void Flags::DeleteFlags()
    {
        for (size_t i = 0; i < strflags.Size(); i++) {
            delete[] strflags[i];
        }
        strflags.DeleteAll();
        numflags.DeleteAll();
        defflags.DeleteAll();
    }

    void Flags::SetFlag(const char* name, const char* val)
    {
        char* hval = new char[strlen(val) + 1];
        strcpy(hval, val);
        strflags.Set(name, hval);
    }

    void Flags::SetFlag(const char* name, double val)
    {
        numflags.Set(name, val);
    }

    void Flags::SetFlag(const char* name)
    {
        defflags.Set(name, 1);
    }

    double Flags::GetNumFlag(const char* name, double def) const
    {
        if (numflags.Used(name))
            return numflags.Get(name);
        else
            return def;
    }


    void Flags::SetCommandLineFlag(const char* st)
    {
        if (st[0] != '-') {
            std::cerr << "flag must start with '-'" << std::endl;
            return;
        }

        const char* pos = strchr(st, '=');

        if (!pos) {
            SetFlag(st + 1);
        } else {
            char name[100];
            strncpy(name, st + 1, (pos - st) - 1);
            name[pos - st - 1] = 0;
            pos++;
            char* endptr = NULL;

            double val = strtod(pos, &endptr);

            if (endptr == pos) {
                SetFlag(name, pos);
            } else {
                SetFlag(name, val);
            }
        }
    }

}  // namespace meshit
