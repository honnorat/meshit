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

#include "mesh_class.hpp"

#include "../gprim/geomtest3d.hpp"

namespace meshit {

int Mesh::CheckSurface()
{
    MESHIT_LOG_DEBUG("Check Surface mesh");

    size_t nf = elements.size();
    IndexPair_map<int> edges(nf + 2);
    IndexPair i2;
    int cnt1 = 0, cnt2 = 0;

    for (size_t i = 0; i < nf; i++) {
        for (size_t j = 0; j < 3; j++) {
            i2.I1() = elements[i].PointID(j);
            i2.I2() = elements[i].PointID((j + 1) % 3);
            if (edges.count(i2)) {
                int hi = edges[i2];
                if (hi != 1) MESHIT_LOG_ERROR("CheckSurface, hi = " << hi);
                edges[i2] = 2;
                cnt2++;
            } else {
                std::swap(i2.I1(), i2.I2());
                edges[i2] = 1;
                cnt1++;
            }
        }
    }

    if (cnt1 != cnt2) {
        MESHIT_LOG_ERROR("Surface mesh not consistent : cnt1 = " << cnt1 << " / cnt2 = " << cnt2);
        return 0;
    }
    return 1;
}

static double TriangleQualityInst(const MeshPoint& p1, const MeshPoint& p2, const MeshPoint& p3)
{
    // quality 0 (worst) .. 1 (optimal)

    Vec2d v1, v2, v3;
    double s1, s2, s3;
    double an1, an2, an3;

    v1 = p2 - p1;
    v2 = p3 - p1;
    v3 = p3 - p2;

    an1 = Angle(v1, v2);
    v1 *= -1;
    an2 = Angle(v1, v3);
    an3 = Angle(v2, v3);

    s1 = sin(an1 / 2);
    s2 = sin(an2 / 2);
    s3 = sin(an3 / 2);

    return 8 * s1 * s2 * s3;
}


void Mesh::PrintQuality()
{
    size_t ncl = 20, cl;
    std::vector<size_t> incl(ncl, 0);
    double qual;

    for (size_t sei = 0; sei < elements.size(); sei++) {
        qual = TriangleQualityInst(points[elements[sei].PointID(0)],
                                   points[elements[sei].PointID(1)],
                                   points[elements[sei].PointID(2)]);

        cl = static_cast<size_t>((ncl - 1e-3) * qual);
        incl[cl]++;
    }

    MESHIT_LOG_INFO("\n\n");
    MESHIT_LOG_INFO("Points:           " << GetNbPoints());
    MESHIT_LOG_INFO("Surface Elements: " << GetNbElements());
    MESHIT_LOG_INFO("\nElements in qualityclasses:");
    for (size_t i = 0; i < ncl; i++) {
        MESHIT_LOG_INFO(std::fixed << std::setprecision(2) << std::setw(4) << static_cast<double>(i) / ncl
                        << " - " << std::setw(4) << static_cast<double>(i + 1) / ncl << ": " << incl[i]);
    }
}


}  // namespace meshit
