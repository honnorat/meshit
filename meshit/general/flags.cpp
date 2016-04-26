/**
 * meshit - a 2d mesh generator
 *
 * Copyright © 1995-2015 Joachim Schoeberl <joachim.schoeberl@tuwien.ac.at>
 * Copyright © 2015-2016 Marc Honnorat <marc.honnorat@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library in the file LICENSE.LGPL; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

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
