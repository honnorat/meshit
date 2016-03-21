#ifndef FILE_FLAGS
#define FILE_FLAGS


/**************************************************************************/
/* File:   flags.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                   */
/**************************************************************************/

#include <iostream>

#include "../meshit.hpp"
#include "symbolta.hpp"

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

     public:
        Flags();
        ~Flags();

        /// Deletes all flags
        void DeleteFlags();

        /// Sets string flag, overwrite if exists
        void SetFlag(const char* name, const char* val);
        void SetFlag(const char* name);
        void SetFlag(const char* name, double val);

        /// set flag of form -name=hello -val=0.5 -defined
        void SetCommandLineFlag(const char* st);

        /// Returns numerical flag, default value if not exists
        double GetNumFlag(const char* name, double def) const;
    };

}  // namespace meshit

#endif

