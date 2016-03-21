#include "refine.hpp"

namespace meshit {

    void Refinement::PointBetween(
            const Point<3>& p1, const Point<3>& p2, double secpoint, Point<3>& newp) const
    {
        newp = p1 + secpoint * (p2 - p1);
    }

    void Refinement::PointBetween(const Point<3>& p1, const Point<3>& p2, double secpoint,
                                  const EdgePointGeomInfo& ap1,
                                  const EdgePointGeomInfo& ap2,
                                  Point<3>& newp, EdgePointGeomInfo& newgi) const
    {
        newp = p1 + secpoint * (p2 - p1);
    }

    void Refinement::Refine(Mesh& mesh)
    {
        // reduce 2nd order
        mesh.ComputeNVertices();
        mesh.SetNP(mesh.GetNV());
        INDEX_2_HASHTABLE<PointIndex> between(mesh.GetNP() + 5);

        // refine edges
        Array<EdgePointGeomInfo> epgi;

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
                Point<3> pnew;
                PointBetween(mesh.Point(el[0]),
                             mesh.Point(el[1]), 0.5,
                             el.epgeominfo[0], el.epgeominfo[1],
                             pnew, ngi);

                pinew = mesh.AddPoint(pnew);
                between.Set(i2, pinew);

                if (pinew >= epgi.size()) {
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

            ArrayMem<PointIndex, 6> pnums(6);

            static int betw[3][3] = {
                    {2, 3, 4},
                    {1, 3, 5},
                    {1, 2, 6}
            };

            for (j = 1; j <= 3; j++) {
                pnums[j - 1] = el.PNum(j);
            }

            for (j = 0; j < 3; j++) {
                PointIndex pi1 = pnums[betw[j][0] - 1];
                PointIndex pi2 = pnums[betw[j][1] - 1];

                INDEX_2 i2(pi1, pi2);
                i2.Sort();

                Point<3> pb;
                PointBetween(mesh.Point(pi1), mesh.Point(pi2), 0.5, pb);

                if (between.Used(i2)) {
                    pnums[3 + j] = between.Get(i2);
                } else {
                    pnums[3 + j] = mesh.AddPoint(pb);
                    between.Set(i2, pnums[3 + j]);
                }
            }

            static int reftab[4][3] = {
                    {1, 6, 5},
                    {2, 4, 6},
                    {3, 5, 4},
                    {6, 4, 5}
            };

            int ind = el.GetIndex();
            for (j = 0; j < 4; j++) {
                Element2d nel;
                for (k = 1; k <= 3; k++) {
                    nel.PNum(k) = pnums[reftab[j][k - 1] - 1];
                }
                nel.SetIndex(ind);

                if (j == 0)
                    mesh.SurfaceElement(sei) = nel;
                else
                    mesh.AddSurfaceElement(nel);
            }
        }

        // update identification tables
        for (size_t i = 0; i < mesh.GetIdentifications().GetMaxNr(); i++) {
            Array<int> identmap;
            mesh.GetIdentifications().GetMap(i + 1, identmap);

            for (size_t j = 0; j < between.GetNBags(); j++)
                for (size_t k = 0; k < between.GetBagSize(j); k++) {
                    INDEX_2 i2;
                    PointIndex newpi;
                    between.GetData(j, k, i2, newpi);
                    INDEX_2 oi2(identmap[i2.I1() - 1], identmap[i2.I2() - 1]);
                    oi2.Sort();
                    if (between.Used(oi2)) {
                        PointIndex onewpi = between.Get(oi2);
                        mesh.GetIdentifications().Add(newpi, onewpi, i + 1);
                    }
                }

        }

        mesh.ComputeNVertices();
        mesh.RebuildSurfaceElementLists();
        return;
    }
}  // namespace meshit
