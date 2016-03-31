#include "refine.hpp"

namespace meshit
{
    void Refinement::PointBetween(const Point3d& p1, const Point3d& p2, double secpoint, Point3d& newp) const
    {
        newp[0] = p1[0] + secpoint * (p2[0] - p1[0]);
        newp[1] = p1[1] + secpoint * (p2[1] - p1[1]);
        newp[2] = p1[2] + secpoint * (p2[2] - p1[2]);
    }

    void Refinement::Refine(Mesh& mesh)
    {
        // reduce 2nd order
        mesh.ComputeNVertices();
        mesh.SetNP(mesh.GetNV());
        INDEX_2_HASHTABLE<PointIndex> between(mesh.GetNP() + 5);

        // refine edges
        std::vector<EdgePointGeomInfo> epgi;

        size_t oldns = mesh.GetNSeg();
        for (size_t si = 0; si < oldns; si++) {
            const Segment& el = mesh.LineSegment(si);

            INDEX_2 i2 = INDEX_2::Sort(el[0], el[1]);
            PointIndex pinew;
            EdgePointGeomInfo ngi;

            if (between.Used(i2)) {
                pinew = between.Get(i2);
                ngi = epgi[pinew];
            } else {
                Point3d pnew;
                PointBetween(mesh.Point(el[0]),
                             mesh.Point(el[1]), 0.5, pnew);
                pinew = mesh.AddPoint(pnew);
                between.Set(i2, pinew);

                if (pinew >= static_cast<PointIndex>(epgi.size())) {
                    epgi.resize(pinew + 1);
                }
                epgi[pinew] = ngi;
            }

            Segment ns1 = el;
            Segment ns2 = el;
            ns1[1] = pinew;
            ns1.epgeominfo[1] = ngi;
            ns2[0] = pinew;
            ns2.epgeominfo[0] = ngi;

            mesh.LineSegment(si) = ns1;
            mesh.AddSegment(ns2);
        }

        // refine surface elements
        size_t oldnf = mesh.GetNSE();
        for (size_t sei = 0; sei < oldnf; sei++) {
            int j, k;
            const Element2d& el = mesh.SurfaceElement(sei);

            PointIndex pnums[6];

            static int betw[3][3] = {
                {1, 2, 3},
                {0, 2, 4},
                {0, 1, 5}
            };

            for (j = 1; j <= 3; j++) {
                pnums[j - 1] = el.PNum(j);
            }

            for (j = 0; j < 3; j++) {
                PointIndex pi1 = pnums[betw[j][0]];
                PointIndex pi2 = pnums[betw[j][1]];

                INDEX_2 i2(pi1, pi2);
                i2.Sort();

                Point3d pb;
                PointBetween(mesh.Point(pi1), mesh.Point(pi2), 0.5, pb);

                if (between.Used(i2)) {
                    pnums[3 + j] = between.Get(i2);
                } else {
                    pnums[3 + j] = mesh.AddPoint(pb);
                    between.Set(i2, pnums[3 + j]);
                }
            }

            static int reftab[4][3] = {
                {0, 5, 4},
                {1, 3, 5},
                {2, 4, 3},
                {5, 3, 4}
            };

            int ind = el.GetIndex();
            for (j = 0; j < 4; j++) {
                Element2d nel;
                for (k = 1; k <= 3; k++) {
                    nel.PNum(k) = pnums[reftab[j][k - 1]];
                }
                nel.SetIndex(ind);

                if (j == 0)
                    mesh.SurfaceElement(sei) = nel;
                else
                    mesh.AddSurfaceElement(nel);
            }
        }

        mesh.ComputeNVertices();
        mesh.RebuildSurfaceElementLists();
        return;
    }
}  // namespace meshit
