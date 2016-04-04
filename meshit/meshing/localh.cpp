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
        h2 = 0.5 * (ax2[0] - ax1[0]);
        for (int i = 0; i < 3; i++)
            xmid[i] = 0.5 * (ax1[i] + ax2[i]);

        for (int i = 0; i < 8; i++) {
            childs[i] = NULL;
        }

        hopt = 2.0 * h2;
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
        for (int i = 0; i < 8; i++)
            if (childs[i]) {
                childs[i]->DeleteChilds();
                delete childs[i];
                childs[i] = NULL;
            }
    }

    LocalH::~LocalH()
    {
        CleanRoot();
    }

    void LocalH::Init(const Point3d& pmin, const Point3d& pmax, double agrading)
    {
        grading = agrading;

        // a small enlargement, non-regular points
        constexpr double val = 0.0879;
        constexpr double val2 = 2 * val;
        constexpr double val3 = 3 * val;

        double x1[3], x2[3];
        x1[0] = (1.0 + val) * pmin.X() - val * pmax.X();
        x2[0] = 1.1 * pmax.X() - 0.1 * pmin.X();

        x1[1] = (1.0 + val2) * pmin.Y() - val2 * pmax.Y();
        x2[1] = 1.1 * pmax.Y() - 0.1 * pmin.Y();

        x1[2] = (1.0 + val3) * pmin.Z() - val3 * pmax.Z();
        x2[2] = 1.1 * pmax.Z() - 0.1 * pmin.Z();

        double hmax = x2[0] - x1[0];
        hmax = std::max(hmax, x2[1] - x1[1]);
        hmax = std::max(hmax, x2[2] - x1[2]);

        x2[0] = x1[0] + hmax;
        x2[1] = x1[1] + hmax;
        x2[2] = x1[2] + hmax;

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

    void LocalH::SetH(const Point3d& p, double h)
    {
        double p_x = p.X();
        double p_y = p.Y();
        double p_z = p.Z();

        if (fabs(p_x - root->xmid[0]) > root->h2 ||
            fabs(p_y - root->xmid[1]) > root->h2 ||
            fabs(p_z - root->xmid[2]) > root->h2) {
            return;
        }

        if (GetH(p) <= 1.2 * h) return;

        GradingBox* box = root;
        GradingBox* nbox = root;
        double x1[3], x2[3];

        while (nbox) {
            box = nbox;
            int childnr = 0;
            childnr += 1 * static_cast<int>((p_x > box->xmid[0]));
            childnr += 2 * static_cast<int>((p_y > box->xmid[1]));
            childnr += 4 * static_cast<int>((p_z > box->xmid[2]));
            nbox = box->childs[childnr];
        }

        const double h_half = 0.5 * h;

        while (box->h2 > h_half) {
            int childnr = 0;
            childnr += 1 * static_cast<int>((p_x > box->xmid[0]));
            childnr += 2 * static_cast<int>((p_y > box->xmid[1]));
            childnr += 4 * static_cast<int>((p_z > box->xmid[2]));

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
            if (childnr & 4) {
                x1[2] = box->xmid[2];
                x2[2] = x1[2] + h2;
            } else {
                x2[2] = box->xmid[2];
                x1[2] = x2[2] - h2;
            }

            box = box->childs[childnr] = new GradingBox(x1, x2);
        }

        box->hopt = h;

        double hbox = 2 * box->h2;  // box->x2[0] - box->x1[0];
        double hnp = h + grading * hbox;

        SetH(Point3d(p_x + hbox, p_y, p_z), hnp);
        SetH(Point3d(p_x - hbox, p_y, p_z), hnp);

        SetH(Point3d(p_x, p_y + hbox, p_z), hnp);
        SetH(Point3d(p_x, p_y - hbox, p_z), hnp);

        SetH(Point3d(p_x, p_y, p_z + hbox), hnp);
        SetH(Point3d(p_x, p_y, p_z - hbox), hnp);
    }

    double LocalH::GetH(const Point3d& x) const
    {
        const GradingBox* box = root;

        while (1) {
            int childnr = 0;
            if (x.X() > box->xmid[0]) childnr += 1;
            if (x.Y() > box->xmid[1]) childnr += 2;
            if (x.Z() > box->xmid[2]) childnr += 4;

            if (box->childs[childnr])
                box = box->childs[childnr];
            else
                return box->hopt;
        }
    }

    double LocalH::GetMinH(const Point3d& pmin, const Point3d& pmax) const
    {
        // minimal h in box (pmin, pmax)
        Point3d pmin2, pmax2;
        for (int j = 1; j <= 3; j++) {
            if (pmin.X(j) < pmax.X(j)) {
                pmin2.X(j) = pmin.X(j);
                pmax2.X(j) = pmax.X(j);
            } else {
                pmin2.X(j) = pmax.X(j);
                pmax2.X(j) = pmin.X(j);
            }
        }
        return GetMinHRec(pmin2, pmax2, root);
    }

    double LocalH::GetMinHRec(
        const Point3d& pmin, const Point3d& pmax,
        const GradingBox* box) const
    {
        double h2 = box->h2;
        if (pmax.X() < box->xmid[0] - h2 || pmin.X() > box->xmid[0] + h2 ||
            pmax.Y() < box->xmid[1] - h2 || pmin.Y() > box->xmid[1] + h2 ||
            pmax.Z() < box->xmid[2] - h2 || pmin.Z() > box->xmid[2] + h2)
            return 1e8;

        double hmin = 2 * box->h2;  // box->x2[0] - box->x1[0];

        for (int i = 0; i < 8; i++) {
            if (box->childs[i])
                hmin = std::min(hmin, GetMinHRec(pmin, pmax, box->childs[i]));
        }

        return hmin;
    }
}  // namespace meshit
