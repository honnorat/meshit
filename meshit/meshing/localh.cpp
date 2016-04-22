#include "localh.hpp"
#include "adfront2.hpp"

namespace meshit {

GradingBox::GradingBox(const Point2d& px1, const Point2d& px2)
{
    SetBox(px1, px2);
}

void GradingBox::SetBox(const Point2d& px1, const Point2d& px2)
{
    hopt = px2.X() - px1.X();
    h2 = 0.5 * hopt;

    xmid = Center(px1, px2);

    childs[0] = nullptr;
    childs[1] = nullptr;
    childs[2] = nullptr;
    childs[3] = nullptr;
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

void LocalH::Init(const Box2d& bbox, double agrading)
{
    grading = agrading;

    // a small enlargement, non-regular points
    constexpr double fact = 0.1;
    constexpr double alpha = 0.5 * (1.0 + fact);

    double hmax = alpha * bbox.LargestSide();

    Point2d pc = bbox.Center();
    Point2d x1(pc.X() - hmax, pc.Y() - hmax);
    Point2d x2(pc.X() + hmax, pc.Y() + hmax);

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

    if (fabs(p_x - root->xmid.X()) > root->h2 || fabs(p_y - root->xmid.Y()) > root->h2) {
        return;
    }

    if (GetH(p) <= h) return;

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
    return GetMinHRec(pmin, pmax, root);
}

double LocalH::GetMinHRec(const Point2d& pmin, const Point2d& pmax, const GradingBox* box) const
{
    double h2 = box->h2;
    if (pmax.X() < box->xmid.X() - h2 || pmin.X() > box->xmid.X() + h2 || pmax.Y() < box->xmid.Y() - h2 ||
        pmin.Y() > box->xmid.Y() + h2) {
        return 1e8;
    }

    double hmin = 2.0 * h2;
    double hmin_0 = 1e100;
    double hmin_1 = 1e100;
    double hmin_2 = 1e100;
    double hmin_3 = 1e100;

    if (box->childs[0]) hmin_0 = GetMinHRec(pmin, pmax, box->childs[0]);
    if (box->childs[1]) hmin_1 = GetMinHRec(pmin, pmax, box->childs[1]);
    if (box->childs[2]) hmin_2 = GetMinHRec(pmin, pmax, box->childs[2]);
    if (box->childs[3]) hmin_3 = GetMinHRec(pmin, pmax, box->childs[3]);

    hmin = std::min(hmin, hmin_0);
    hmin = std::min(hmin, hmin_1);
    hmin = std::min(hmin, hmin_2);
    hmin = std::min(hmin, hmin_3);

    return hmin;
}

}  // namespace meshit
