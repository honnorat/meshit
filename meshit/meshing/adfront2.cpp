/*
  Advancing front class for surfaces
 */
#include <iostream>
#include "adfront2.hpp"

namespace meshit {

    FrontPoint2::FrontPoint2(const Point3d& ap, PointIndex agi,
                             MultiPointGeomInfo* amgi, bool aonsurface)
            : globalindex(agi)
    {
        p = ap;
        nlinetopoint = 0;
        frontnr = INT_MAX - 10;
        onsurface = aonsurface;

        if (amgi) {
            mgi = new MultiPointGeomInfo(*amgi);
            for (int i = 1; i <= mgi->GetNPGI(); i++)
                if (mgi->GetPGI(i).trignum <= 0)
                    std::cout << "Add FrontPoint2, illegal geominfo = " << mgi->GetPGI(i).trignum << std::endl;
        }
        else
            mgi = NULL;
    }

    AdFront2::AdFront2(const Box3d& aboundingbox)
            : boundingbox(aboundingbox),
              linesearchtree(boundingbox.PMin(), boundingbox.PMax()),
              pointsearchtree(boundingbox.PMin(), boundingbox.PMax()),
              cpointsearchtree(boundingbox.PMin(), boundingbox.PMax())
    {
        nfl = 0;
        allflines = 0;

        minval = 0;
        starti = 0;
    }

    AdFront2::~AdFront2()
    {
        delete allflines;
    }

    void AdFront2::PrintOpenSegments(std::ostream& ost) const
    {
        if (nfl > 0) {
            ost << nfl << " open front segments left:" << std::endl;
            for (size_t i = 0; i < lines.size(); i++)
                if (lines[i].Valid())
                    ost << i << ": "
                    << GetGlobalIndex(lines[i].L().I1()) << "-"
                    << GetGlobalIndex(lines[i].L().I2()) << std::endl;
        }
    }

    int AdFront2::AddPoint(
            const Point3d& p, PointIndex globind,
            MultiPointGeomInfo* mgi,
            bool pointonsurface)
    {
        // inserts at empty position or resizes array
        int pi;

        if (delpointl.size() != 0) {
            pi = delpointl.Last();
            delpointl.DeleteLast();

            points[pi] = FrontPoint2(p, globind, mgi, pointonsurface);
        } else {
            points.push_back(FrontPoint2(p, globind, mgi, pointonsurface));
            pi = points.size() - 1;
        }

        if (mgi)
            cpointsearchtree.Insert(p, pi);

        if (pointonsurface)
            pointsearchtree.Insert(p, pi);

        return pi;
    }

    int AdFront2::AddLine(int pi1, int pi2,
                          const PointGeomInfo& gi1, const PointGeomInfo& gi2)
    {
        int minfn;
        int li;

        FrontPoint2& p1 = points[pi1];
        FrontPoint2& p2 = points[pi2];

        nfl++;

        p1.AddLine();
        p2.AddLine();

        minfn = std::min(p1.FrontNr(), p2.FrontNr());
        p1.DecFrontNr(minfn + 1);
        p2.DecFrontNr(minfn + 1);

        if (dellinel.size() != 0) {
            li = dellinel.Last();
            dellinel.DeleteLast();
            lines[li] = FrontLine(INDEX_2(pi1, pi2));
        } else {
            lines.push_back(FrontLine(INDEX_2(pi1, pi2)));
            li = lines.size() - 1;
        }

        if (!gi1.trignum || !gi2.trignum) {
            std::cout << "ERROR: in AdFront::AddLine, illegal geominfo" << std::endl;
        }

        lines[li].SetGeomInfo(gi1, gi2);

        Box3d lbox;
        lbox.SetPoint(p1.P());
        lbox.AddPoint(p2.P());

        linesearchtree.Insert(lbox.PMin(), lbox.PMax(), li);

        if (allflines) {
            if (allflines->Used(INDEX_2(GetGlobalIndex(pi1),
                                        GetGlobalIndex(pi2)))) {
                std::cerr << "ERROR Adfront2::AddLine: line exists" << std::endl;
            }

            allflines->Set(INDEX_2(GetGlobalIndex(pi1),
                                   GetGlobalIndex(pi2)), 1);
        }

        return li;
    }

    void AdFront2::DeleteLine(int li)
    {
        int pi;

        nfl--;

        for (int i = 1; i <= 2; i++) {
            pi = lines[li].L().I(i);
            points[pi].RemoveLine();

            if (!points[pi].Valid()) {
                delpointl.push_back(pi);
                if (points[pi].mgi) {
                    cpointsearchtree.DeleteElement(pi);
                    delete points[pi].mgi;
                    points[pi].mgi = NULL;
                }

                pointsearchtree.DeleteElement(pi);
            }
        }

        if (allflines) {
            allflines->Set(INDEX_2(GetGlobalIndex(lines[li].L().I1()),
                                   GetGlobalIndex(lines[li].L().I2())), 2);
        }

        lines[li].Invalidate();
        linesearchtree.DeleteElement(li);

        dellinel.push_back(li);
    }

    int AdFront2::ExistsLine(int pi1, int pi2)
    {
        if (!allflines)
            return 0;
        if (allflines->Used(INDEX_2(pi1, pi2)))
            return allflines->Get(INDEX_2(pi1, pi2));
        else
            return 0;
    }

    int AdFront2::SelectBaseLine(Point<3>& p1, Point<3>& p2,
                                 const PointGeomInfo*& geominfo1,
                                 const PointGeomInfo*& geominfo2,
                                 int& qualclass)
    {
        int baselineindex = -1;

        for (size_t i = starti; i < lines.size(); i++) {
            if (lines[i].Valid()) {
                int hi = lines[i].LineClass() +
                         points[lines[i].L().I1()].FrontNr() +
                         points[lines[i].L().I2()].FrontNr();

                if (hi <= minval) {
                    minval = hi;
                    baselineindex = i;
                    break;
                }
            }
        }

        if (baselineindex == -1) {
            minval = INT_MAX;
            for (size_t i = 0; i < lines.size(); i++)
                if (lines[i].Valid()) {
                    int hi = lines[i].LineClass() +
                             points[lines[i].L().I1()].FrontNr() +
                             points[lines[i].L().I2()].FrontNr();

                    if (hi < minval) {
                        minval = hi;
                        baselineindex = i;
                    }
                }
        }
        starti = baselineindex + 1;

        p1 = points[lines[baselineindex].L().I1()].P();
        p2 = points[lines[baselineindex].L().I2()].P();
        geominfo1 = &lines[baselineindex].GetGeomInfo(1);
        geominfo2 = &lines[baselineindex].GetGeomInfo(2);

        qualclass = lines[baselineindex].LineClass();

        return baselineindex;
    }

    int AdFront2::GetLocals(
            int baselineindex,
            Array<Point3d>& locpoints,
            Array<MultiPointGeomInfo>& pgeominfo,
            Array<INDEX_2>& loclines,  // local index
            Array<INDEX>& pindex,
            Array<INDEX>& lindex,
            double xh)
    {
        int pstind;
        Point<3> p0;

        pstind = lines[baselineindex].L().I1();
        p0 = points[pstind].P();

        loclines.push_back(lines[baselineindex].L());
        lindex.push_back(baselineindex);

        ArrayMem<int, 1000> nearlines(0);
        ArrayMem<int, 1000> nearpoints(0);

        // dominating costs !!
        linesearchtree.GetIntersecting(
                p0 - Vec3d(xh, xh, xh),
                p0 + Vec3d(xh, xh, xh),
                nearlines);

        pointsearchtree.GetIntersecting(
                p0 - Vec3d(xh, xh, xh),
                p0 + Vec3d(xh, xh, xh),
                nearpoints);

        for (int ii = 0; ii < nearlines.size(); ii++) {
            int i = nearlines[ii];
            if (lines[i].Valid() && i != baselineindex) {
                loclines.push_back(lines[i].L());
                lindex.push_back(i);
            }
        }

        // static Array<int> invpindex;
        invpindex.resize(points.size());
        // invpindex = -1;
        for (int i = 0; i < nearpoints.size(); i++) {
            invpindex[nearpoints[i]] = -1;
        }

        for (int i = 0; i < loclines.size(); i++) {
            invpindex[loclines[i].I1()] = 0;
            invpindex[loclines[i].I2()] = 0;
        }

        for (int i = 0; i < loclines.size(); i++) {
            for (int j = 0; j < 2; j++) {
                int pi = loclines[i][j];
                if (invpindex[pi] == 0) {
                    pindex.push_back(pi);
                    invpindex[pi] = pindex.size();
                    locpoints.push_back(points[pi].P());
                    loclines[i][j] = locpoints.size();
                }
                else
                    loclines[i][j] = invpindex[pi];
            }
        }

        // double xh2 = xh*xh;
        for (int ii = 0; ii < nearpoints.size(); ii++) {
            int i = nearpoints[ii];
            if (points[i].Valid() && points[i].OnSurface() && invpindex[i] <= 0) {
                locpoints.push_back(points[i].P());
                invpindex[i] = locpoints.size();
                pindex.push_back(i);
            }
        }

        pgeominfo.resize(locpoints.size());
        for (int i = 0; i < pgeominfo.size(); i++) {
            pgeominfo[i].Init();
        }

        for (int i = 0; i < loclines.size(); i++) {
            for (int j = 0; j < 2; j++) {
                int lpi = loclines[i][j];

                const PointGeomInfo& gi =
                        lines[lindex[i]].GetGeomInfo(j + 1);
                pgeominfo.Elem(lpi).AddPointGeomInfo(gi);
            }
        }

        for (int i = 0; i < locpoints.size(); i++) {
            int pi = pindex[i];

            if (points[pi].mgi) {
                for (int j = 1; j <= points[pi].mgi->GetNPGI(); j++) {
                    pgeominfo[i].AddPointGeomInfo(points[pi].mgi->GetPGI(j));
                }
            }
        }
        return lines[baselineindex].LineClass();
    }

    void AdFront2::SetStartFront()
    {
        for (size_t i = 0; i < lines.size(); i++) {
            if (lines[i].Valid()) {
                for (int j = 1; j <= 2; j++) {
                    points[lines[i].L().I(j)].DecFrontNr(0);
                }
            }
        }
    }

    bool AdFront2::Inside(const Point<2>& p) const
    {
        int cnt;
        Vec<2> n;
        Vec<3> v1;
        DenseMatrix a(2), ainv(2);
        Vector b(2), u(2);

        // quasi-random numbers:
        n(0) = 0.123871;
        n(1) = 0.15432;

        cnt = 0;
        for (size_t i = 0; i < lines.size(); i++) {
            if (lines[i].Valid()) {
                const Point<3>& p1 = points[lines[i].L().I1()].P();
                const Point<3>& p2 = points[lines[i].L().I2()].P();

                v1 = p2 - p1;

                a(0, 0) = v1(0);
                a(1, 0) = v1(1);

                a(0, 1) = -n(0);
                a(1, 1) = -n(1);

                b(0) = p(0) - p1(0);
                b(1) = p(1) - p1(1);

                CalcInverse(a, ainv);
                ainv.Mult(b, u);

                if (u(0) >= 0 && u(0) <= 1 && u(1) > 0)
                    cnt++;
            }
        }
        return ((cnt % 2) != 0);
    }
}  // namespace meshit
