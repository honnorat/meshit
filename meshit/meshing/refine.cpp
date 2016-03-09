#include "refine.hpp"

namespace meshit {

    void Refinement::PointBetween(
            const Point<3>& p1, const Point<3>& p2, double secpoint,
            const PointGeomInfo& gi1,
            const PointGeomInfo& gi2,
            Point<3>& newp, PointGeomInfo& newgi) const
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

        int oldns, oldnf;

        // refine edges

        Array<EdgePointGeomInfo> epgi;

        oldns = mesh.GetNSeg();
        for (SegmentIndex si = 0; si < oldns; si++) {
            const Segment& el = mesh.LineSegment(si);

            INDEX_2 i2 = INDEX_2::Sort(el[0], el[1]);
            PointIndex pinew;
            EdgePointGeomInfo ngi;

            if (between.Used(i2)) {
                pinew = between.Get(i2);
                ngi = epgi[pinew];
            }
            else {
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
        Array<PointGeomInfo> surfgi(8 * mesh.GetNP());
        for (int i = 0; i < surfgi.size(); i++) {
            surfgi[i].trignum = -1;
        }

        oldnf = mesh.GetNSE();
        for (SurfaceElementIndex sei = 0; sei < oldnf; sei++) {
            int j, k;
            const Element2d& el = mesh.SurfaceElement(sei);

            switch (el.GetType()) {
                case TRIG:
                case TRIG6: {
                    ArrayMem<PointIndex, 6> pnums(6);
                    ArrayMem<PointGeomInfo, 6> pgis(6);

                    static int betw[3][3] = {
                            {2, 3, 4},
                            {1, 3, 5},
                            {1, 2, 6}
                    };

                    for (j = 1; j <= 3; j++) {
                        pnums[j - 1] = el.PNum(j);
                        pgis[j - 1] = el.GeomInfoPi(j);
                    }

                    for (j = 0; j < 3; j++) {
                        PointIndex pi1 = pnums[betw[j][0] - 1];
                        PointIndex pi2 = pnums[betw[j][1] - 1];

                        INDEX_2 i2(pi1, pi2);
                        i2.Sort();

                        Point<3> pb;
                        PointGeomInfo pgi;
                        PointBetween(mesh.Point(pi1),
                                     mesh.Point(pi2), 0.5,
                                     el.GeomInfoPi(betw[j][0]),
                                     el.GeomInfoPi(betw[j][1]),
                                     pb, pgi);

                        pgis[4 + j - 1] = pgi;
                        if (between.Used(i2))
                            pnums[4 + j - 1] = between.Get(i2);
                        else {
                            pnums[4 + j - 1] = mesh.AddPoint(pb);
                            between.Set(i2, pnums.Get(4 + j));
                        }

                        if (surfgi.size() < pnums[4 + j - 1])
                            surfgi.resize(pnums[4 + j - 1]);
                        surfgi[pnums[4 + j - 1] - 1] = pgis[4 + j - 1];
                    }

                    static int reftab[4][3] = {
                            {1, 6, 5},
                            {2, 4, 6},
                            {3, 5, 4},
                            {6, 4, 5}
                    };

                    int ind = el.GetIndex();
                    for (j = 0; j < 4; j++) {
                        Element2d nel(TRIG);
                        for (k = 1; k <= 3; k++) {
                            nel.PNum(k) = pnums.Get(reftab[j][k - 1]);
                            nel.GeomInfoPi(k) = pgis.Get(reftab[j][k - 1]);
                        }
                        nel.SetIndex(ind);

                        if (j == 0)
                            mesh.SurfaceElement(sei) = nel;
                        else
                            mesh.AddSurfaceElement(nel);
                    }
                    break;
                }
                case QUAD:
                case QUAD6:
                case QUAD8: {
                    ArrayMem<PointIndex, 9> pnums(9);
                    ArrayMem<PointGeomInfo, 9> pgis(9);

                    static int betw[5][3] = {
                            {1, 2, 5},
                            {2, 3, 6},
                            {3, 4, 7},
                            {1, 4, 8},
                            {5, 7, 9}
                    };

                    for (j = 1; j <= 4; j++) {
                        pnums[j - 1] = el.PNum(j);
                        pgis[j - 1] = el.GeomInfoPi(j);
                    }

                    for (j = 0; j < 5; j++) {
                        int pi1 = pnums[betw[j][0] - 1];
                        int pi2 = pnums[betw[j][1] - 1];

                        INDEX_2 i2(pi1, pi2);
                        i2.Sort();

                        if (between.Used(i2)) {
                            pnums[4 + j] = between.Get(i2);
                            pgis[4 + j] = surfgi.Get(pnums[4 + j - 1]);
                        }
                        else {
                            Point<3> pb;
                            PointBetween(mesh.Point(pi1),
                                         mesh.Point(pi2), 0.5,
                                         el.GeomInfoPi(betw[j][0]),
                                         el.GeomInfoPi(betw[j][1]),
                                         pb, pgis[4 + j]);

                            pnums[4 + j] = mesh.AddPoint(pb);

                            between.Set(i2, pnums.Get(5 + j));

                            if (surfgi.size() < pnums[4 + j])
                                surfgi.resize(pnums[4 + j]);
                            surfgi[pnums[4 + j] - 1] = pgis[4 + j];
                        }
                    }

                    static int reftab[4][4] = {
                            {1, 5, 9, 8},
                            {5, 2, 6, 9},
                            {8, 9, 7, 4},
                            {9, 6, 3, 7}
                    };

                    int ind = el.GetIndex();
                    for (j = 0; j < 4; j++) {
                        Element2d nel(QUAD);
                        for (k = 1; k <= 4; k++) {
                            nel.PNum(k) = pnums.Get(reftab[j][k - 1]);
                            nel.GeomInfoPi(k) = pgis.Get(reftab[j][k - 1]);
                        }
                        nel.SetIndex(ind);

                        if (j == 0)
                            mesh.SurfaceElement(sei) = nel;
                        else
                            mesh.AddSurfaceElement(nel);
                    }
                    break;
                }
                default:
                    MESHIT_LOG_ERROR("Refine: undefined surface element type " << int(el.GetType()));
            }
        }

        // update identification tables
        for (int i = 1; i <= mesh.GetIdentifications().GetMaxNr(); i++) {
            Array<int> identmap;
            mesh.GetIdentifications().GetMap(i, identmap);

            for (int j = 1; j <= between.GetNBags(); j++)
                for (int k = 1; k <= between.GetBagSize(j); k++) {
                    INDEX_2 i2;
                    PointIndex newpi;
                    between.GetData(j, k, i2, newpi);
                    INDEX_2 oi2(identmap.Get(i2.I1()),
                                identmap.Get(i2.I2()));
                    oi2.Sort();
                    if (between.Used(oi2)) {
                        PointIndex onewpi = between.Get(oi2);
                        mesh.GetIdentifications().Add(newpi, onewpi, i);
                    }
                }

        }

        mesh.ComputeNVertices();
        mesh.RebuildSurfaceElementLists();
        return;
    }
}
