#include "localh.hpp"
#include "adfront2.hpp"

namespace meshit
{
    GradingBox::GradingBox(const Point2d& px1, const Point2d& px2)
    {
        SetBox(px1, px2);
    }

    void GradingBox::SetBox(const Point2d& px1, const Point2d& px2)
    {
        hopt = Dist(px2, px1);
        h2 = 0.5 * hopt;

        xmid = Center(px1, px2);

        for (int i = 0; i < 4; i++) {
            childs[i] = nullptr;
        }
    }

    BlockAllocator GradingBox::ball(sizeof(GradingBox));

    void* GradingBox::operator new(size_t)
    {
        return ball.Alloc();
    }

    void GradingBox::operator delete(void* p)
    {
        ball.Free(p);
    }

    void GradingBox::DeleteChilds()
    {
        for (int i = 0; i < 4; i++)
            if (childs[i]) {
                childs[i]->DeleteChilds();
                delete childs[i];
                childs[i] = nullptr;
            }
    }

    LocalH::~LocalH()
    {
        CleanRoot();
    }

    void LocalH::Init(const Point2d& pmin, const Point2d& pmax, double agrading)
    {
        grading = agrading;

        // a small enlargement, non-regular points
        constexpr double alpha = 0.13;
        constexpr double al_p1 = 1.0 + alpha;

        double hmax = std::max(pmax.X() - pmin.X(),
                               pmax.Y() - pmin.Y());

        Point2d x1, x2;

        x1.X() = al_p1 * pmin.X() - alpha * pmax.X();
        x1.Y() = al_p1 * pmin.Y() - alpha * pmax.Y();

        x2.X() = x1.X() + hmax;
        x2.Y() = x1.Y() + hmax;

        CleanRoot();
        root = new GradingBox(x1, x2);
    }

    void LocalH::CleanRoot()
    {
        if (root) {
            root->DeleteChilds();
            delete root;
        }
    }

    void LocalH::SetH(const Point2d& p, double h)
    {
        double p_x = p.X();
        double p_y = p.Y();

        if (fabs(p_x - root->xmid.X()) > root->h2 ||
            fabs(p_y - root->xmid.Y()) > root->h2) {
            return;
        }

        if (GetH(p) <= 1.2 * h) return;

        GradingBox* box = root;
        GradingBox* nbox = root;
        Point2d x1, x2;

        while (nbox) {
            box = nbox;
            int childnr = 0;
            childnr += 1 * static_cast<int>((p_x > box->xmid.X()));
            childnr += 2 * static_cast<int>((p_y > box->xmid.Y()));
            nbox = box->childs[childnr];
        }

        const double h_half = 0.5 * h;

        while (box->h2 > h_half) {
            int childnr = 0;
            childnr += 1 * static_cast<int>((p_x > box->xmid.X()));
            childnr += 2 * static_cast<int>((p_y > box->xmid.Y()));

            double h2 = box->h2;
            if (childnr & 1) {
                x1.X() = box->xmid.X();
                x2.X() = x1.X() + h2;
            } else {
                x2.X() = box->xmid.X();
                x1.X() = x2.X() - h2;
            }
            if (childnr & 2) {
                x1.Y() = box->xmid.Y();
                x2.Y() = x1.Y() + h2;
            } else {
                x2.Y() = box->xmid.Y();
                x1.Y() = x2.Y() - h2;
            }

            box = box->childs[childnr] = new GradingBox(x1, x2);
        }

        box->hopt = h;

        double hbox = 2 * box->h2;  // box->x2[0] - box->x1[0];
        double hnp = h + grading * hbox;

        SetH(Point2d(p_x + hbox, p_y), hnp);
        SetH(Point2d(p_x - hbox, p_y), hnp);

        SetH(Point2d(p_x, p_y + hbox), hnp);
        SetH(Point2d(p_x, p_y - hbox), hnp);
    }

    double LocalH::GetH(const Point2d& x) const
    {
        const GradingBox* box = root;

        while (true) {
            int childnr = 0;
            if (x.X() > box->xmid.X()) childnr += 1;
            if (x.Y() > box->xmid.Y()) childnr += 2;

            if (box->childs[childnr])
                box = box->childs[childnr];
            else
                return box->hopt;
        }
    }

    double LocalH::GetMinH(const Point2d& pmin, const Point2d& pmax) const
    {
        // minimal h in box (pmin, pmax)
        Point2d pmin2, pmax2;

        pmin2.X() = std::min(pmin.X(), pmax.X());
        pmax2.X() = std::max(pmin.X(), pmax.X());
        pmin2.Y() = std::min(pmin.Y(), pmax.Y());
        pmax2.Y() = std::max(pmin.Y(), pmax.Y());

        return GetMinHRec(pmin2, pmax2, root);
    }

    double LocalH::GetMinHRec(const Point2d& pmin, const Point2d& pmax, const GradingBox* box) const
    {
        double h2 = box->h2;
        if (pmax.X() < box->xmid.X() - h2 || pmin.X() > box->xmid.X() + h2 ||
            pmax.Y() < box->xmid.Y() - h2 || pmin.Y() > box->xmid.Y() + h2) {
            return 1e8;
        }

        double hmin = 2.0 * h2;  // box->x2[0] - box->x1[0];

        for (int i = 0; i < 4; i++) {
            if (box->childs[i]) {
                hmin = std::min(hmin, GetMinHRec(pmin, pmax, box->childs[i]));
            }
        }
        return hmin;
    }

}  // namespace meshit
