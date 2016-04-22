/**************************************************************************/
/* File:   flags.cc                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

/*
   Datatype Flags
*/

#include "flags.hpp"

#include <string>

#include "logging.hpp"

namespace meshit {

double Flags::GetNumFlag(const std::string& name, double default_value)
{
    if (num_flags.count(name)) {
        return num_flags[name];
    } else {
        return default_value;
    }
}

int Flags::GetIntFlag(const std::string& name, int default_value)
{
    if (num_flags.count(name)) {
        return static_cast<int>(num_flags[name]);
    } else {
        return default_value;
    }
}

void Flags::SetCommandLineFlag(const std::string& st)
{
    size_t eq_pos = st.find('=');
    if (eq_pos == std::string::npos) {
        MESHIT_LOG_WARNING("invalid flag: " << st);
        return;
    }

    double value;
    try {
        value = std::stod(st.substr(eq_pos + 1));
    } catch (const std::invalid_argument& e) {
        MESHIT_LOG_WARNING("Invalid argument for flag: '" << st << "' (" << e.what() << ").");
        return;
    }

    std::string name = st.substr(0, eq_pos);
    num_flags[name] = value;
}

}  // namespace meshit
