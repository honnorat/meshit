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

#include "geomtest3d.hpp"

#include "../linalg/densemat.hpp"

namespace meshit {

bool IntersectTriangleLine(const Point2d** tri, const Point2d** line)
{
    Vec2d vl(*line[0], *line[1]);
    Vec2d vt1(*tri[0], *tri[1]);
    Vec2d vt2(*tri[0], *tri[2]);
    Vec2d vrs(*tri[0], *line[0]);

    Mat3x3 a, ainv;
    Vec3d rs = Vec3d(vrs);
    Vec3d lami;

    a = 0;
    a(0, 0) = -vl.X();
    a(0, 1) = vt1.X();
    a(0, 2) = vt2.X();
    a(1, 0) = -vl.Y();
    a(1, 1) = vt1.Y();
    a(1, 2) = vt2.Y();

    double det = a.Det();
    double arel = vl.Length() * vt1.Length() * vt2.Length();

    if (fabs(det) <= 1e-10 * arel) {
        return false;
    }

    a.CalcInverse(ainv);
    lami = ainv * rs;

    return lami.X() >= 0 && lami.X() <= 1 && lami.Y() >= 0 && lami.Z() >= 0 && lami.Y() + lami.Z() <= 1;
}

bool IntersectTriangleTriangle(const Point2d** tri1, const Point2d** tri2)
{
    constexpr double epsrel = 1e-8;
    const double diam = Dist(*tri1[0], *tri1[1]);
    const double eps = diam * epsrel;
    const double eps2 = eps * eps;

    uint32_t cnt = 0;
    for (uint32_t i = 0; i <= 2; i++) {
        for (uint32_t j = 0; j <= 2; j++) {
            if (Dist2(*tri1[j], *tri2[i]) < eps2) {
                cnt++;
                break;
            }
        }
    }
    switch (cnt) {
        case 0: {
            const Point2d* line[2];
            for (size_t i = 0; i <= 2; i++) {
                line[0] = tri2[i];
                line[1] = tri2[(i + 1) % 3];

                if (IntersectTriangleLine(tri1, &line[0])) {
                    std::cerr << "int1, line = " << *line[0] << " - " << *line[1] << std::endl;
                    return true;
                }
            }
            for (size_t i = 0; i <= 2; i++) {
                line[0] = tri1[i];
                line[1] = tri1[(i + 1) % 3];

                if (IntersectTriangleLine(tri2, &line[0])) {
                    std::cerr << "int2, line = " << *line[0] << " - " << *line[1] << std::endl;
                    return true;
                }
            }
            break;
        }
        default:
            return false;
    }
    return false;
}

}  // namespace meshit
