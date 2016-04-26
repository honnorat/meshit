#ifndef MESHIT_FLAGS_HPP
#define MESHIT_FLAGS_HPP
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
