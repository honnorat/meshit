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

#include "geom3d.hpp"

namespace meshit {

inline Vec3d& Vec3d::operator+=(const Vec3d& v2)
{
    x[0] += v2.X();
    x[1] += v2.Y();
    x[2] += v2.Z();
    return *this;
}

inline Vec3d& Vec3d::operator-=(const Vec3d& v2)
{
    x[0] -= v2.X();
    x[1] -= v2.Y();
    x[2] -= v2.Z();
    return *this;
}

inline Vec3d& Vec3d::operator*=(double s)
{
    x[0] *= s;
    x[1] *= s;
    x[2] *= s;
    return *this;
}

inline Vec3d& Vec3d::operator/=(double s)
{
    if (s != 0) {
        x[0] /= s;
        x[1] /= s;
        x[2] /= s;
    }
    return *this;
}

inline Vec3d operator-(const Vec3d& v1, const Vec3d& v2)
{
    return Vec3d(v1.x[0] - v2.x[0], v1.x[1] - v2.x[1], v1.x[2] - v2.x[2]);
}

inline Vec3d operator+(const Vec3d& v1, const Vec3d& v2)
{
    return Vec3d(v1.x[0] + v2.x[0], v1.x[1] + v2.x[1], v1.x[2] + v2.x[2]);
}

inline Vec3d operator*(double scal, const Vec3d& v)
{
    return Vec3d(scal * v.x[0], scal * v.x[1], scal * v.x[2]);
}

inline double operator*(const Vec3d& v1, const Vec3d& v2)
{
    return v1.x[0] * v2.x[0] + v1.x[1] * v2.x[1] + v1.x[2] * v2.x[2];
}

Vec3d operator*(const Mat3x3& m, const Vec3d& v)
{
    Vec3d res;
    res.X() = m(0, 0) * v.X() + m(0, 1) * v.Y() + m(0, 2) * v.Z();
    res.Y() = m(1, 0) * v.X() + m(1, 1) * v.Y() + m(1, 2) * v.Z();
    res.Z() = m(2, 0) * v.X() + m(2, 1) * v.Y() + m(2, 2) * v.Z();
    return res;
}

std::ostream& operator<<(std::ostream& s, const Vec3d& v)
{
    return s << "(" << v.x[0] << ", " << v.x[1] << ", " << v.x[2] << ")";
}

}  // namespace meshit
