#include "localh.hpp"
#include "adfront2.hpp"

namespace meshit
{
    GradingBox::GradingBox(const double* ax1, const double* ax2)
    {
        SetBox(ax1, ax2);
    }

    void GradingBox::SetBox(const double* ax1, const double* ax2)
    {
        hopt = ax2[0] - ax1[0];
        h2 = 0.5 * hopt;

        xmid[0] = 0.5 * (ax1[0] + ax2[0]);
        xmid[1] = 0.5 * (ax1[1] + ax2[1]);

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

        double x1[2], x2[2];
        double hmax = std::max(pmax.X() - pmin.X(),
                               pmax.Y() - pmin.Y());

        x1[0] = al_p1 * pmin.X() - alpha * pmax.X();
        x1[1] = al_p1 * pmin.Y() - alpha * pmax.Y();

        x2[0] = x1[0] + hmax;
        x2[1] = x1[1] + hmax;

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

        if (fabs(p_x - root->xmid[0]) > root->h2 ||
            fabs(p_y - root->xmid[1]) > root->h2) {
            return;
        }

        if (GetH(p) <= 1.2 * h) return;

        GradingBox* box = root;
        GradingBox* nbox = root;
        double x1[2], x2[2];

        while (nbox) {
            box = nbox;
            int childnr = 0;
            childnr += 1 * static_cast<int>((p_x > box->xmid[0]));
            childnr += 2 * static_cast<int>((p_y > box->xmid[1]));
            nbox = box->childs[childnr];
        }

        const double h_half = 0.5 * h;

        while (box->h2 > h_half) {
            int childnr = 0;
            childnr += 1 * static_cast<int>((p_x > box->xmid[0]));
            childnr += 2 * static_cast<int>((p_y > box->xmid[1]));

            double h2 = box->h2;
            if (childnr & 1) {
                x1[0] = box->xmid[0];
                x2[0] = x1[0] + h2;
            } else {
                x2[0] = box->xmid[0];
                x1[0] = x2[0] - h2;
            }
            if (childnr & 2) {
                x1[1] = box->xmid[1];
                x2[1] = x1[1] + h2;
            } else {
                x2[1] = box->xmid[1];
                x1[1] = x2[1] - h2;
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
            if (x.X() > box->xmid[0]) childnr += 1;
            if (x.Y() > box->xmid[1]) childnr += 2;

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
        if (pmax.X() < box->xmid[0] - h2 || pmin.X() > box->xmid[0] + h2 ||
            pmax.Y() < box->xmid[1] - h2 || pmin.Y() > box->xmid[1] + h2) {
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
