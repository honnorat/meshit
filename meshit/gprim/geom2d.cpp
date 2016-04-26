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

#define _USE_MATH_DEFINES 1

#include "geom2d.hpp"
#include "geom3d.hpp"

namespace meshit {

inline const Box2d& Box2d::operator+=(const Box2d& b)
{
    pmin.X() = std::min(pmin.X(), b.pmin.X());
    pmin.Y() = std::min(pmin.Y(), b.pmin.Y());
    pmax.X() = std::max(pmax.X(), b.pmax.X());
    pmax.Y() = std::max(pmax.Y(), b.pmax.Y());
    return *this;
}

std::ostream& operator<<(std::ostream& s, const Point2d& p)
{
    return s << "(" << p.px << ", " << p.py << ")";
}

std::ostream& operator<<(std::ostream& s, const Vec2d& v)
{
    return s << "(" << v.vx << ", " << v.vy << ")";
}

}  // namespace meshit
