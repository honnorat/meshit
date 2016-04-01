#ifndef FILE_FLAGS
#define FILE_FLAGS


/**************************************************************************/
/* File:   flags.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                   */
/**************************************************************************/

#include <iostream>
#include <unordered_map>

#include <cstring>

namespace meshit
{
/**
   Flag - Table.
   A flag table maintains string variables, numerical 
   variables and boolean flags.
*/
    class Flags
    {
     protected:
        std::unordered_map<std::string, double> num_flags;

     public:
        Flags() { }

        ~Flags() { }

        void SetCommandLineFlag(const std::string& st);

        /// Returns numerical flag, default value if not exists
        double GetNumFlag(const std::string& name, double default_value);
    };

}  // namespace meshit

#endif

