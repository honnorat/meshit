/**************************************************************************/
/* File:   flags.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

/* 
   Datatype Flags
*/

#include <sstream>
#include <fstream>
#include "flags.hpp"

namespace meshit {

    Flags::Flags() { }

    Flags::~Flags()
    {
        DeleteFlags();
    }

    void Flags::DeleteFlags()
    {
        for (int i = 0; i < strflags.Size(); i++)
            delete[] strflags[i];
        for (int i = 0; i < numlistflags.Size(); i++)
            delete numlistflags[i];
        strflags.DeleteAll();
        numflags.DeleteAll();
        defflags.DeleteAll();
        strlistflags.DeleteAll();
        numlistflags.DeleteAll();
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

    const char*
    Flags::GetStringFlag(const char* name, const char* def) const
    {
        if (strflags.Used(name))
            return strflags.Get(name);
        else
            return def;
    }

    double Flags::GetNumFlag(const char* name, double def) const
    {
        if (numflags.Used(name))
            return numflags.Get(name);
        else
            return def;
    }

    bool Flags::GetDefineFlag(const char* name) const
    {
        return defflags.Used(name);
    }

    bool Flags::StringFlagDefined(const char* name) const
    {
        return strflags.Used(name);
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
