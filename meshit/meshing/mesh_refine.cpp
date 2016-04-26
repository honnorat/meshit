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

namespace meshit {

void Mesh::Refine()
{
    // reduce 2nd order
    size_t nb_vertices = ComputeNVertices();

    points.resize(nb_vertices);
    IndexPair_map<PointIndex> between(nb_vertices + 5);

    size_t oldns = segments.size();
    for (size_t si = 0; si < oldns; si++) {
        const Segment& seg = segments[si];
        const MeshPoint& p1 = points[seg[0]];
        const MeshPoint& p2 = points[seg[1]];

        IndexPair i2(seg[0], seg[1]);
        i2.Sort();
        PointIndex pi_new;

        if (between.count(i2) == 1) {
            pi_new = between[i2];
        } else {
            Point2d pnew;
            pnew.X() = 0.5 * (p1.X() + p2.X());
            pnew.Y() = 0.5 * (p1.Y() + p2.Y());
            pi_new = AddPoint(pnew);
            between[i2] = pi_new;
        }

        Segment ns1 = seg;
        Segment ns2 = seg;
        ns1[1] = pi_new;
        ns2[0] = pi_new;

        segments[si] = ns1;
        AddSegment(ns2);
    }

    // refine surface elements
    size_t oldnf = elements.size();
    for (size_t sei = 0; sei < oldnf; sei++) {
        int j, k;
        const Element2d& el = elements[sei];

        PointIndex pnums[6];

        static int betw[3][3] = {{1, 2, 3},
                                 {0, 2, 4},
                                 {0, 1, 5}};

        for (j = 0; j < 3; j++) {
            pnums[j] = el.PointID(j);
        }

        for (j = 0; j < 3; j++) {
            PointIndex pi1 = pnums[betw[j][0]];
            PointIndex pi2 = pnums[betw[j][1]];

            IndexPair i2(pi1, pi2);
            i2.Sort();

            if (between.count(i2) == 0) {
                const MeshPoint& p1 = points[pi1];
                const MeshPoint& p2 = points[pi2];
                Point2d pnew;
                pnew.X() = 0.5 * (p1.X() + p2.X());
                pnew.Y() = 0.5 * (p1.Y() + p2.Y());
                between[i2] = AddPoint(pnew);
            }
            pnums[3 + j] = between[i2];
        }

        static int reftab[4][3] = {{0, 5, 4},
                                   {1, 3, 5},
                                   {2, 4, 3},
                                   {5, 3, 4}};

        DomainIndex ind = el.FaceID();
        for (j = 0; j < 4; j++) {
            Element2d new_element;
            for (k = 0; k < 3; k++) {
                new_element.PointID(k) = pnums[reftab[j][k]];
            }
            new_element.SetFaceID(ind);

            if (j == 0) {
                elements[sei] = new_element;
            } else {
                AddElement(new_element);
            }
        }
    }

    ComputeNVertices();
    RebuildSurfaceElementLists();
}

}  // namespace meshit
