#ifndef FILE_FLAGS_HPP
#define FILE_FLAGS_HPP

/**************************************************************************/
/* File:   flags.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                   */
/**************************************************************************/

#include <iostream>
#include <unordered_map>

#include <cstring>

namespace meshit {
/**
   Flag - Table.
   A flag table maintains string variables, numerical
   variables and boolean flags.
*/
class Flags
{
 public:
    Flags() { }

    ~Flags() { }

    void SetCommandLineFlag(const std::string& st);

    // Returns numerical flag, default value if not exists
    double GetNumFlag(const std::string& name, double default_value);
    int GetIntFlag(const std::string& name, int default_value);

 protected:
    std::unordered_map<std::string, double> num_flags;
};

}  // namespace meshit

#endif
