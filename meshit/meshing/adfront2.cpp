/*
  Advancing front class for surfaces
 */
#include <iostream>
#include "adfront2.hpp"

namespace meshit
{
    FrontPoint2::FrontPoint2(const Point2d& ap, PointIndex agi)
        : globalindex(agi)
    {
        p = ap;
        nlinetopoint = 0;
        frontnr = INT_MAX - 10;
    }

    AdFront2::AdFront2(const Box2d& aboundingbox)
        : boundingbox(aboundingbox),
          linesearchtree(boundingbox.PMin(), boundingbox.PMax()),
          pointsearchtree(boundingbox.PMin(), boundingbox.PMax()),
          cpointsearchtree(boundingbox.PMin(), boundingbox.PMax())
    {
        nfl = 0;
        minval = 0;
        starti = 0;
        nearlines.reserve(1000);
        nearpoints.reserve(1000);
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

    int AdFront2::AddPoint(const Point2d& p, PointIndex globind)
    {
        // inserts at empty position or resizes array
        int pi;

        if (delpointl.size() != 0) {
            pi = delpointl[delpointl.size() - 1];
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
            li = dellinel[dellinel.size() - 1];
            dellinel.pop_back();
            lines[li] = FrontLine(INDEX_2(pi1, pi2));
        } else {
            lines.push_back(FrontLine(INDEX_2(pi1, pi2)));
            li = lines.size() - 1;
        }

        Box2d lbox;
        lbox.SetPoint(p1.P());
        lbox.AddPoint(p2.P());

        linesearchtree.Insert(lbox.PMin(), lbox.PMax(), li);

        INDEX_2 globline(GetGlobalIndex(pi1), GetGlobalIndex(pi2));
        if (allflines.count(globline)) {
            MESHIT_LOG_ERROR("Adfront2::AddLine: line exists");
        }
        allflines[globline] = 1;

        return li;
    }

    void AdFront2::DeleteLine(int li)
    {
        nfl--;

        const INDEX_2& lidx = lines[li].L();

        INDEX pi1 = lidx.I1();
        points[pi1].RemoveLine();
        if (!points[pi1].Valid()) {
            delpointl.push_back(pi1);
            pointsearchtree.DeleteElement(pi1);
        }
        INDEX pi2 = lidx.I2();
        points[pi2].RemoveLine();
        if (!points[pi2].Valid()) {
            delpointl.push_back(pi2);
            pointsearchtree.DeleteElement(pi2);
        }

        allflines[INDEX_2(GetGlobalIndex(pi1),
                          GetGlobalIndex(pi2))] = 2;

        lines[li].Invalidate();
        linesearchtree.DeleteElement(li);

        dellinel.push_back(li);
    }

    int AdFront2::ExistsLine(int pi1, int pi2)
    {
        INDEX_2 line(pi1, pi2);
        if (allflines.count(line) == 1) {
            return allflines[line];
        } else {
            return 0;
        }
    }

    int AdFront2::SelectBaseLine(Point2d& p1, Point2d& p2, int& qualclass)
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
        starti = static_cast<size_t>(baselineindex + 1);

        p1 = points[lines[baselineindex].L().I1()].P();
        p2 = points[lines[baselineindex].L().I2()].P();
        qualclass = lines[baselineindex].LineClass();

        return baselineindex;
    }

    int AdFront2::GetLocals(int baselineindex,
                            std::vector<Point2d>& locpoints,
                            std::vector<INDEX_2>& loclines,  // local index
                            std::vector<INDEX>& pindex,
                            std::vector<INDEX>& lindex,
                            double xh)
    {
        int pstind = lines[baselineindex].L().I1();
        Point2d p0 = points[pstind].P();

        loclines.push_back(lines[baselineindex].L());
        lindex.push_back(baselineindex);

        const Point2d pmin(p0.X() - xh, p0.Y() - xh);
        const Point2d pmax(p0.X() + xh, p0.Y() + xh);
        // dominating costs !!
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
                } else {
                    loclines[i][j] = invpindex[pi];
                }
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
                points[lines[i].L().I1()].DecFrontNr(0);
                points[lines[i].L().I2()].DecFrontNr(0);
            }
        }
    }

}  // namespace meshit
