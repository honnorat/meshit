#include "adtree.hpp"
#include <cstring>

namespace meshit {

#define ADTREE_MAX_STACK_SIZE 20

/* ******************************* ADTree3 ******************************* */

ADTreeNode3::ADTreeNode3()
{
    pi = -1;

    left = nullptr;
    right = nullptr;
    father = nullptr;
    nchilds = 0;
}

void ADTreeNode3::DeleteChilds()
{
    if (left) {
        left->DeleteChilds();
        delete left;
        left = nullptr;
    }
    if (right) {
        right->DeleteChilds();
        delete right;
        right = nullptr;
    }
}

BlockAllocator ADTreeNode3::ball(sizeof(ADTreeNode3));

void* ADTreeNode3::operator new(size_t /*s*/)
{
    return ball.Alloc();
}

void ADTreeNode3::operator delete(void* p)
{
    ball.Free(p);
}

ADTree3::ADTree3(const Point2d& acmin, const Point2d& acmax)
{
    cmin[0] = acmin.X();
    cmin[1] = acmin.Y();
    cmax[0] = acmax.X();
    cmax[1] = acmax.Y();

    root = new ADTreeNode3;
    root->sep = (cmin[0] + cmax[0]) / 2;
}

ADTree3::~ADTree3()
{
    root->DeleteChilds();
    delete root;
}

void ADTree3::Insert(const Point2d& p, PointIndex pi)
{
    ADTreeNode3* node = nullptr;
    ADTreeNode3* next = root;
    int dir = 0;
    bool lr = false;

    double bmin[3];
    double bmax[3];

    memcpy(bmin, cmin, 3 * sizeof(double));
    memcpy(bmax, cmax, 3 * sizeof(double));

    while (next) {
        node = next;

        if (node->pi == -1) {
            node->data[0] = p[0];
            node->data[1] = p[1];
            node->pi = pi;

            if (ela.size() < pi + 1) {
                ela.resize(pi + 1);
            }
            ela[pi] = node;

            return;
        }

        if (node->sep > p[dir]) {
            next = node->left;
            bmax[dir] = node->sep;
            lr = false;
        } else {
            next = node->right;
            bmin[dir] = node->sep;
            lr = true;
        }

        if (++dir == 2) {
            dir = 0;
        }
    }

    next = new ADTreeNode3;
    next->data[0] = p[0];
    next->data[1] = p[1];
    next->pi = pi;
    next->sep = (bmin[dir] + bmax[dir]) / 2;

    if (ela.size() < static_cast<size_t>(pi + 1)) {
        ela.resize(pi + 1);
    }
    ela[pi] = next;

    if (lr)
        node->right = next;
    else
        node->left = next;
    next->father = node;

    while (node) {
        node->nchilds++;
        node = node->father;
    }
}

void ADTree3::DeleteElement(PointIndex pi)
{
    ADTreeNode3* node = ela[pi];

    node->pi = CONST<PointIndex>::undefined;

    node = node->father;
    while (node) {
        node->nchilds--;
        node = node->father;
    }
}

struct inttn3
{
    ADTreeNode3* node;
    int dir;
};

void ADTree3::GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const
{
    pis.clear();

    static inttn3 stack[ADTREE_MAX_STACK_SIZE];
    stack[0].node = root;
    stack[0].dir = 0;
    size_t stacks = 1;

    while (stacks) {
        stacks--;
        ADTreeNode3* node = stack[stacks].node;

        if (node->pi != -1 &&
            !(node->data[0] < pmin.X() || node->data[0] > pmax.X() || node->data[1] < pmin.Y() ||
              node->data[1] > pmax.Y())) {
            pis.push_back(node->pi);
        }

        int dir = stack[stacks].dir;

        if (node->left && pmin[dir] <= node->sep) {
            stack[stacks].node = node->left;
            stack[stacks].dir = (dir + 1) % 2;
            stacks++;
        }
        if (node->right && pmax[dir] >= node->sep) {
            stack[stacks].node = node->right;
            stack[stacks].dir = (dir + 1) % 2;
            stacks++;
        }
    }
}

void ADTree3::PrintRec(std::ostream& ost, const ADTreeNode3* node) const
{
    if (node->data) {
        ost << node->pi << ": ";
        ost << node->nchilds << " childs, ";
        ost << node->data[0] << " ";
        ost << node->data[1] << " ";
        ost << std::endl;
    }
    if (node->left) PrintRec(ost, node->left);
    if (node->right) PrintRec(ost, node->right);
}

/* ******************************* ADTree6 ******************************* */

ADTreeNode6::ADTreeNode6()
{
    pi_ = CONST<PointIndex>::undefined;

    left = nullptr;
    right = nullptr;
    father = nullptr;
    nchilds = 0;
}

void ADTreeNode6::DeleteChilds()
{
    if (left) {
        left->DeleteChilds();
        delete left;
        left = nullptr;
    }
    if (right) {
        right->DeleteChilds();
        delete right;
        right = nullptr;
    }
}

inline void ADTreeNode6::SetData(const Point2d& pmin, const Point2d& pmax, PointIndex pi)
{
    pmin_ = pmin;
    pmax_ = pmax;
    pi_ = pi;
}

BlockAllocator ADTreeNode6::ball(sizeof(ADTreeNode6));

void* ADTreeNode6::operator new(size_t /*s*/)
{
    return ball.Alloc();
}

void ADTreeNode6::operator delete(void* p)
{
    ball.Free(p);
}

ADTree6::ADTree6(const Point2d& pmin, const Point2d& pmax)
    : pmin_{pmin}, pmax_{pmax}
{
    root = new ADTreeNode6;
    root->sep_ = 0.5 * (pmin_.X() + pmax_.X());
}

ADTree6::~ADTree6()
{
    root->DeleteChilds();
    delete root;
}

void ADTree6::Insert(const Point2d& bmin_, const Point2d& bmax_, PointIndex pi)
{
    ADTreeNode6* node = nullptr;
    ADTreeNode6* next = root;
    int dir = 0;
    bool lr = false;

    double bmin[4] = {pmin_.X(), pmin_.Y(), pmin_.X(), pmin_.Y()};
    double bmax[4] = {pmax_.X(), pmax_.Y(), pmax_.X(), pmax_.Y()};
    double p[4] = {bmin_.X(), bmin_.Y(), bmax_.X(), bmax_.Y()};

    while (next) {
        node = next;

        if (node->pi_ == CONST<PointIndex>::undefined) {
            node->SetData(bmin_, bmax_, pi);

            if (pi >= ela.size()) {
                ela.push_back(node);
            } else {
                ela[pi] = node;
            }
            return;
        }

        if (node->sep_ > p[dir]) {
            next = node->left;
            bmax[dir] = node->sep_;
            lr = false;
        } else {
            next = node->right;
            bmin[dir] = node->sep_;
            lr = true;
        }

        if (++dir == 4) dir = 0;
    }

    next = new ADTreeNode6;
    next->SetData(bmin_, bmax_, pi);
    next->sep_ = (bmin[dir] + bmax[dir]) / 2;

    if (pi >= ela.size()) {
        ela.push_back(next);
    } else {
        ela[pi] = next;
    }

    if (lr)
        node->right = next;
    else
        node->left = next;
    next->father = node;

    while (node) {
        node->nchilds++;
        node = node->father;
    }
}

void ADTree6::DeleteElement(PointIndex pi)
{
    ADTreeNode6* node = ela[pi];

    node->pi_ = CONST<PointIndex>::undefined;

    node = node->father;
    while (node) {
        node->nchilds--;
        node = node->father;
    }
}

struct inttn6
{
    ADTreeNode6* node;
    int dir;
};

void ADTree6::GetIntersecting(const Point2d& bmin, const Point2d& bmax,
                              const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const
{
    pis.clear();

    const double abmin[4] = {bmin.X(), bmin.Y(), pmin.X(), pmin.Y()};
    const double abmax[4] = {pmax.X(), pmax.Y(), bmax.X(), bmax.Y()};

    static inttn6 stack[ADTREE_MAX_STACK_SIZE];
    stack[0].node = root;
    stack[0].dir = 0;
    size_t stacks = 1;

    while (stacks) {
        stacks--;
        ADTreeNode6* node = stack[stacks].node;
        if (node->pi_ != CONST<PointIndex>::undefined &&
            !(node->pmin_.X() > pmax.X() || node->pmax_.X() < pmin.X() ||
              node->pmin_.Y() > pmax.Y() || node->pmax_.Y() < pmin.Y())) {
            pis.push_back(node->pi_);
        }

        int dir = stack[stacks].dir;

        if (node->left && abmin[dir] <= node->sep_) {
            stack[stacks].node = node->left;
            stack[stacks].dir = (dir + 1) % 4;
            stacks++;
        }
        if (node->right && abmax[dir] >= node->sep_) {
            stack[stacks].node = node->right;
            stack[stacks].dir = (dir + 1) % 4;
            stacks++;
        }
    }
}

/* ************************************* Point3dTree ********************** */

Point3dTree::Point3dTree(const Point2d& pmin, const Point2d& pmax)
{
    tree = new ADTree3(pmin, pmax);
}

Point3dTree::~Point3dTree()
{
    delete tree;
}

void Point3dTree::Insert(const Point2d& p, PointIndex pi)
{
    tree->Insert(p, pi);
}

void Point3dTree::GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const
{
    tree->GetIntersecting(pmin, pmax, pis);
}

Box3dTree::Box3dTree(const Point2d& apmin, const Point2d& apmax)
    : boxpmin{apmin}, boxpmax{apmax}
{
    tree = new ADTree6(apmin, apmax);
}

Box3dTree::~Box3dTree()
{
    delete tree;
}

void Box3dTree::Insert(const Point2d& bmin, const Point2d& bmax, PointIndex pi)
{
    tree->Insert(bmin, bmax, pi);
}

void Box3dTree::GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<size_t>& pis) const
{
    tree->GetIntersecting(boxpmin, boxpmax, pmin, pmax, pis);
}

}  // namespace meshit
