#ifndef FILE_FLAGS
#define FILE_FLAGS


/**************************************************************************/
/* File:   flags.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                   */
/**************************************************************************/

#include "../meshit.hpp"
#include "symbolta.hpp"

#include <iostream>

namespace meshit {

/** 
   Flag - Table.
   A flag table maintains string variables, numerical 
   variables and boolean flags.
*/
    class Flags
    {
     protected:
        SYMBOLTABLE<char*> strflags;
        SYMBOLTABLE<double> numflags;
        SYMBOLTABLE<int> defflags;
        SYMBOLTABLE<Array<char*>*> strlistflags;
        SYMBOLTABLE<Array<double>*> numlistflags;
     public:
        DLL_HEADER Flags();
        DLL_HEADER ~Flags();

        /// Deletes all flags
        DLL_HEADER void DeleteFlags();

        /// Sets string flag, overwrite if exists
        DLL_HEADER void SetFlag(const char* name, const char* val);
        DLL_HEADER void SetFlag(const char* name);
        DLL_HEADER void SetFlag(const char* name, double val);

        /// set flag of form -name=hello -val=0.5 -defined
        DLL_HEADER void SetCommandLineFlag(const char* st);

        /// Returns string flag, default value if not exists
        DLL_HEADER const char* GetStringFlag(const char* name, const char* def) const;
        /// Returns numerical flag, default value if not exists
        DLL_HEADER double GetNumFlag(const char* name, double def) const;
        /// Returns boolean flag
        DLL_HEADER bool GetDefineFlag(const char* name) const;


        /// Test, if string flag is defined
        DLL_HEADER bool StringFlagDefined(const char* name) const;
    };

}

#endif

