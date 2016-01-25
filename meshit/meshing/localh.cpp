#include <meshit/meshit.hpp>
#include "localh.hpp"
#include "adfront2.hpp"
#include "../gprim/geomfuncs.hpp"

namespace meshit {

    GradingBox::GradingBox(const double * ax1, const double * ax2)
    {
        h2 = 0.5 * (ax2[0] - ax1[0]);
        for (int i = 0; i < 3; i++)
            xmid[i] = 0.5 * (ax1[i] + ax2[i]);

        for (int i = 0; i < 8; i++) {
            childs[i] = NULL;
        }
        father = NULL;

        flags.cutboundary = 0;
        flags.isinner = 0;
        flags.oldcell = 0;
        flags.pinner = 0;

        hopt = 2.0 * h2;
    }

    BlockAllocator GradingBox::ball(sizeof (GradingBox));

    void * GradingBox::operator new(size_t)
    {
        return ball.Alloc();
    }

    void GradingBox::operator delete (void * p)
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

    LocalH::LocalH(const Point3d & pmin, const Point3d & pmax, double agrading)
    {
        double x1[3], x2[3];
        double hmax;

        boundingbox = Box3d(pmin, pmax);
        grading = agrading;

        // a small enlargement, non-regular points 
        double val = 0.0879;
        for (int i = 1; i <= 3; i++) {
            x1[i - 1] = (1 + val * i) * pmin.X(i) - val * i * pmax.X(i);
            x2[i - 1] = 1.1 * pmax.X(i) - 0.1 * pmin.X(i);
        }

        hmax = x2[0] - x1[0];
        for (int i = 1; i <= 2; i++)
            if (x2[i] - x1[i] > hmax)
                hmax = x2[i] - x1[i];

        for (int i = 0; i <= 2; i++)
            x2[i] = x1[i] + hmax;

        root = new GradingBox(x1, x2);
        boxes.push_back(root);
    }

    LocalH::LocalH(const Box<3> & box, double agrading)
    {
        Point3d pmin = box.PMin();
        Point3d pmax = box.PMax();

        double x1[3], x2[3];
        double hmax;

        boundingbox = Box3d(pmin, pmax);
        grading = agrading;

        // a small enlargement, non-regular points 
        double val = 0.0879;
        for (int i = 1; i <= 3; i++) {
            x1[i - 1] = (1 + val * i) * pmin.X(i) - val * i * pmax.X(i);
            x2[i - 1] = 1.1 * pmax.X(i) - 0.1 * pmin.X(i);
        }

        hmax = x2[0] - x1[0];
        for (int i = 1; i <= 2; i++)
            if (x2[i] - x1[i] > hmax)
                hmax = x2[i] - x1[i];

        for (int i = 0; i <= 2; i++)
            x2[i] = x1[i] + hmax;

        root = new GradingBox(x1, x2);
        boxes.push_back(root);
    }

    LocalH::~LocalH()
    {
        root->DeleteChilds();
        delete root;
    }

    void LocalH::Delete()
    {
        root->DeleteChilds();
    }

    void LocalH::SetH(const Point3d & p, double h)
    {
        if (
                fabs(p.X() - root->xmid[0]) > root->h2 ||
                fabs(p.Y() - root->xmid[1]) > root->h2 ||
                fabs(p.Z() - root->xmid[2]) > root->h2)
            return;

        if (GetH(p) <= 1.2 * h) return;

        GradingBox * box = root;
        GradingBox * nbox = root;
        GradingBox * ngb;
        int childnr;
        double x1[3], x2[3];

        while (nbox) {
            box = nbox;
            childnr = 0;
            if (p.X() > box->xmid[0]) childnr += 1;
            if (p.Y() > box->xmid[1]) childnr += 2;
            if (p.Z() > box->xmid[2]) childnr += 4;
            nbox = box->childs[childnr];
        };

        const double h_half = 0.5 * h;

        while (box->h2 > h_half) {
            childnr = 0;
            if (p.X() > box->xmid[0]) childnr += 1;
            if (p.Y() > box->xmid[1]) childnr += 2;
            if (p.Z() > box->xmid[2]) childnr += 4;

            double h2 = box->h2;
            if (childnr & 1) {
                x1[0] = box->xmid[0];
                x2[0] = x1[0] + h2; // box->x2[0];
            }
            else {
                x2[0] = box->xmid[0];
                x1[0] = x2[0] - h2; // box->x1[0];
            }

            if (childnr & 2) {
                x1[1] = box->xmid[1];
                x2[1] = x1[1] + h2; // box->x2[1];
            }
            else {
                x2[1] = box->xmid[1];
                x1[1] = x2[1] - h2; // box->x1[1];
            }

            if (childnr & 4) {
                x1[2] = box->xmid[2];
                x2[2] = x1[2] + h2; // box->x2[2];
            }
            else {
                x2[2] = box->xmid[2];
                x1[2] = x2[2] - h2; // box->x1[2];
            }

            ngb = new GradingBox(x1, x2);
            box->childs[childnr] = ngb;
            ngb->father = box;

            boxes.push_back(ngb);
            box = box->childs[childnr];
        }

        box->hopt = h;

        double hbox = 2 * box->h2; // box->x2[0] - box->x1[0];
        double hnp = h + grading * hbox;

        Point3d np;
        for (int i = 1; i <= 3; i++) {
            np = p;
            np.X(i) = p.X(i) + hbox;
            SetH(np, hnp);

            np.X(i) = p.X(i) - hbox;
            SetH(np, hnp);
        }
    }

    double LocalH::GetH(const Point3d & x) const
    {
        const GradingBox * box = root;

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

    double LocalH::GetMinH(const Point3d & pmin, const Point3d & pmax) const
    {
        // minimal h in box (pmin, pmax)

        Point3d pmin2, pmax2;
        for (int j = 1; j <= 3; j++) {
            if (pmin.X(j) < pmax.X(j)) {
                pmin2.X(j) = pmin.X(j);
                pmax2.X(j) = pmax.X(j);
            }
            else {
                pmin2.X(j) = pmax.X(j);
                pmax2.X(j) = pmin.X(j);
            }
        }
        return GetMinHRec(pmin2, pmax2, root);
    }

    double LocalH::GetMinHRec(
            const Point3d & pmin, const Point3d & pmax,
            const GradingBox * box) const
    {
        double h2 = box->h2;
        if (
                pmax.X() < box->xmid[0] - h2 || pmin.X() > box->xmid[0] + h2 ||
                pmax.Y() < box->xmid[1] - h2 || pmin.Y() > box->xmid[1] + h2 ||
                pmax.Z() < box->xmid[2] - h2 || pmin.Z() > box->xmid[2] + h2)
            return 1e8;

        double hmin = 2 * box->h2; // box->x2[0] - box->x1[0];

        for (int i = 0; i < 8; i++) {
            if (box->childs[i])
                hmin = std::min(hmin, GetMinHRec(pmin, pmax, box->childs[i]));
        }

        return hmin;
    }
}
