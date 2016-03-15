/*
  Advancing front class for surfaces
 */
#include <iostream>
#include "adfront2.hpp"

namespace meshit {

    FrontPoint2::FrontPoint2(const Point3d& ap, PointIndex agi)
            : globalindex(agi)
    {
        p = ap;
        nlinetopoint = 0;
        frontnr = INT_MAX - 10;
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

    int AdFront2::AddPoint(const Point3d& p, PointIndex globind)
    {
        // inserts at empty position or resizes array
        int pi;

        if (delpointl.size() != 0) {
            pi = delpointl.Last();
            delpointl.pop_back();
            points[pi] = FrontPoint2(p, globind);
        } else {
            points.push_back(FrontPoint2(p, globind));
            pi = points.size() - 1;
        }
        pointsearchtree.Insert(p, pi);

        return pi;
    }

    int AdFront2::AddLine(int pi1, int pi2)
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
            dellinel.pop_back();
            lines[li] = FrontLine(INDEX_2(pi1, pi2));
        } else {
            lines.push_back(FrontLine(INDEX_2(pi1, pi2)));
            li = lines.size() - 1;
        }

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

        for (size_t i = 1; i <= 2; i++) {
            pi = lines[li].L().I(i);
            points[pi].RemoveLine();

            if (!points[pi].Valid()) {
                delpointl.push_back(pi);
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

    int AdFront2::SelectBaseLine(Point3d& p1, Point3d& p2, int& qualclass)
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
            for (size_t i = 0; i < lines.size(); i++) {
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
        }
        starti = baselineindex + 1;

        p1 = points[lines[baselineindex].L().I1()].P();
        p2 = points[lines[baselineindex].L().I2()].P();
        qualclass = lines[baselineindex].LineClass();

        return baselineindex;
    }

    int AdFront2::GetLocals(int baselineindex,
                            Array<Point3d>& locpoints,
                            Array<INDEX_2>& loclines,  // local index
                            Array<INDEX>& pindex,
                            Array<INDEX>& lindex,
                            double xh)
    {
        int pstind;
        Point3d p0;

        pstind = lines[baselineindex].L().I1();
        p0 = points[pstind].P();

        loclines.push_back(lines[baselineindex].L());
        lindex.push_back(baselineindex);

        std::vector<size_t> nearlines;
        std::vector<size_t> nearpoints;

        nearlines.reserve(1000);
        nearpoints.reserve(1000);

        // dominating costs !!
        const Point3d pmin(p0.X() - xh, p0.Y() - xh);
        const Point3d pmax(p0.X() + xh, p0.Y() + xh);
        linesearchtree.GetIntersecting(pmin, pmax, nearlines);
        pointsearchtree.GetIntersecting(pmin, pmax, nearpoints);

        for (size_t ii = 0; ii < nearlines.size(); ii++) {
            int i = nearlines[ii];
            if (lines[i].Valid() && i != baselineindex) {
                loclines.push_back(lines[i].L());
                lindex.push_back(i);
            }
        }

        invpindex.resize(points.size());
        for (size_t i = 0; i < nearpoints.size(); i++) {
            invpindex[nearpoints[i]] = -1;
        }

        for (size_t i = 0; i < loclines.size(); i++) {
            invpindex[loclines[i].I1()] = 0;
            invpindex[loclines[i].I2()] = 0;
        }

        for (size_t i = 0; i < loclines.size(); i++) {
            for (size_t j = 0; j < 2; j++) {
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

        for (size_t ii = 0; ii < nearpoints.size(); ii++) {
            size_t i = nearpoints[ii];
            if (points[i].Valid() && invpindex[i] <= 0) {
                locpoints.push_back(points[i].P());
                invpindex[i] = locpoints.size();
                pindex.push_back(i);
            }
        }
        return lines[baselineindex].LineClass();
    }

    void AdFront2::SetStartFront()
    {
        for (size_t i = 0; i < lines.size(); i++) {
            if (lines[i].Valid()) {
                for (size_t j = 1; j <= 2; j++) {
                    points[lines[i].L().I(j)].DecFrontNr(0);
                }
            }
        }
    }

}  // namespace meshit
