#include "adtree.hpp"
#include <cstring>

namespace meshit {

#define ADTREE_MAX_STACK_SIZE 20

/* ******************************* ADTree3 ******************************* */

void ADTreeNode3::DeleteChilds()
{
    if (left_) {
        left_->DeleteChilds();
        delete left_;
        left_ = nullptr;
    }
    if (right_) {
        right_->DeleteChilds();
        delete right_;
        right_ = nullptr;
    }
}

void ADTreeNode3::SetData(const Point2d& point, GenericIndex pi)
{
    data_[0] = point.X();
    data_[1] = point.Y();
    pi_ = pi;
}

BlockAllocator ADTreeNode3::ball_(sizeof(ADTreeNode3));

void* ADTreeNode3::operator new(size_t /*s*/)
{
    return ball_.Alloc();
}

void ADTreeNode3::operator delete(void* p)
{
    ball_.Free(p);
}

ADTree3::ADTree3(const Point2d& cmin, const Point2d& cmax)
    : cmin_{cmin}, cmax_{cmax}
{
    root_ = new ADTreeNode3;
    root_->sep_ = 0.5 * (cmin_.X() + cmax_.X());
}

ADTree3::~ADTree3()
{
    root_->DeleteChilds();
    delete root_;
}

void ADTree3::Insert(const Point2d& p, GenericIndex pi)
{
    ADTreeNode3* node = nullptr;
    ADTreeNode3* next = root_;
    int dir = 0;
    bool lr = false;

    double bmin[2] = {cmin_.X(), cmin_.Y()};
    double bmax[2] = {cmax_.X(), cmax_.Y()};

    while (next) {
        node = next;

        if (node->pi_ == CONST<GenericIndex>::undefined) {
            node->SetData(p, pi);

            if (pi >= nodes_.size()) {
                nodes_.push_back(node);
            } else {
                nodes_[pi] = node;
            }
            return;
        }
        if (node->sep_ > p[dir]) {
            next = node->left_;
            bmax[dir] = node->sep_;
            lr = false;
        } else {
            next = node->right_;
            bmin[dir] = node->sep_;
            lr = true;
        }
        if (++dir == 2) dir = 0;
    }

    next = new ADTreeNode3;
    next->SetData(p, pi);
    next->sep_ = 0.5 * (bmin[dir] + bmax[dir]);

    if (pi >= nodes_.size()) {
        nodes_.push_back(next);
    } else {
        nodes_[pi] = next;
    }

    if (lr) {
        node->right_ = next;
    } else {
        node->left_ = next;
    }
}

void ADTree3::DeleteElement(GenericIndex pi)
{
    nodes_[pi]->pi_ = CONST<PointIndex>::undefined;
}

struct inttn3
{
    ADTreeNode3* node;
    int dir;
};

void ADTree3::GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<GenericIndex>& pis) const
{
    static inttn3 stack[ADTREE_MAX_STACK_SIZE];
    stack[0].node = root_;
    stack[0].dir = 0;
    size_t stacks = 1;

    pis.clear();
    while (stacks > 0) {
        stacks--;
        ADTreeNode3* node = stack[stacks].node;

        if (node->pi_ != CONST<GenericIndex>::undefined &&
            !(node->data_[0] < pmin.X() || node->data_[0] > pmax.X() ||
              node->data_[1] < pmin.Y() || node->data_[1] > pmax.Y())) {
            pis.push_back(node->pi_);
        }

        int dir = stack[stacks].dir;

        if (node->left_ && pmin[dir] <= node->sep_) {
            stack[stacks].node = node->left_;
            stack[stacks].dir = (dir + 1) % 2;
            stacks++;
        }
        if (node->right_ && pmax[dir] >= node->sep_) {
            stack[stacks].node = node->right_;
            stack[stacks].dir = (dir + 1) % 2;
            stacks++;
        }
    }
}

/* ******************************* ADTree6 ******************************* */

void ADTreeNode6::DeleteChilds()
{
    if (left_) {
        left_->DeleteChilds();
        delete left_;
        left_ = nullptr;
    }
    if (right_) {
        right_->DeleteChilds();
        delete right_;
        right_ = nullptr;
    }
}

inline void ADTreeNode6::SetData(const Point2d& pmin, const Point2d& pmax, GenericIndex pi)
{
    pmin_ = pmin;
    pmax_ = pmax;
    pi_ = pi;
}

BlockAllocator ADTreeNode6::ball_(sizeof(ADTreeNode6));

void* ADTreeNode6::operator new(size_t /*s*/)
{
    return ball_.Alloc();
}

void ADTreeNode6::operator delete(void* p)
{
    ball_.Free(p);
}

ADTree6::ADTree6(const Point2d& pmin, const Point2d& pmax)
    : pmin_{pmin}, pmax_{pmax}
{
    root_ = new ADTreeNode6;
    root_->sep_ = 0.5 * (pmin_.X() + pmax_.X());
}

ADTree6::~ADTree6()
{
    root_->DeleteChilds();
    delete root_;
}

void ADTree6::Insert(const Point2d& bmin_, const Point2d& bmax_, GenericIndex pi)
{
    ADTreeNode6* node = nullptr;
    ADTreeNode6* next = root_;
    int dir = 0;
    bool lr = false;

    double bmin[4] = {pmin_.X(), pmin_.Y(), pmin_.X(), pmin_.Y()};
    double bmax[4] = {pmax_.X(), pmax_.Y(), pmax_.X(), pmax_.Y()};
    double p[4] = {bmin_.X(), bmin_.Y(), bmax_.X(), bmax_.Y()};

    while (next) {
        node = next;

        if (node->pi_ == CONST<GenericIndex>::undefined) {
            node->SetData(bmin_, bmax_, pi);

            if (pi >= nodes_.size()) {
                nodes_.push_back(node);
            } else {
                nodes_[pi] = node;
            }
            return;
        }

        if (node->sep_ > p[dir]) {
            next = node->left_;
            bmax[dir] = node->sep_;
            lr = false;
        } else {
            next = node->right_;
            bmin[dir] = node->sep_;
            lr = true;
        }

        if (++dir == 4) dir = 0;
    }

    next = new ADTreeNode6;
    next->SetData(bmin_, bmax_, pi);
    next->sep_ = 0.5 * (bmin[dir] + bmax[dir]);

    if (pi >= nodes_.size()) {
        nodes_.push_back(next);
    } else {
        nodes_[pi] = next;
    }

    if (lr) {
        node->right_ = next;
    } else {
        node->left_ = next;
    }
}

void ADTree6::DeleteElement(GenericIndex pi)
{
    nodes_[pi]->pi_ = CONST<GenericIndex>::undefined;
}

struct inttn6
{
    ADTreeNode6* node;
    int dir;
};

void ADTree6::GetIntersecting(const Point2d& bmin, const Point2d& bmax,
                              const Point2d& pmin, const Point2d& pmax, std::vector<GenericIndex>& pis) const
{
    const double abmin[4] = {bmin.X(), bmin.Y(), pmin.X(), pmin.Y()};
    const double abmax[4] = {pmax.X(), pmax.Y(), bmax.X(), bmax.Y()};

    static inttn6 stack[ADTREE_MAX_STACK_SIZE];
    stack[0].node = root_;
    stack[0].dir = 0;
    size_t stacks = 1;

    pis.clear();
    while (stacks > 0) {
        stacks--;
        ADTreeNode6* node = stack[stacks].node;
        if (node->pi_ != CONST<PointIndex>::undefined &&
            !(node->pmin_.X() > pmax.X() || node->pmax_.X() < pmin.X() ||
              node->pmin_.Y() > pmax.Y() || node->pmax_.Y() < pmin.Y())) {
            pis.push_back(node->pi_);
        }

        int dir = stack[stacks].dir;

        if (node->left_ && abmin[dir] <= node->sep_) {
            stack[stacks].node = node->left_;
            stack[stacks].dir = (dir + 1) % 4;
            stacks++;
        }
        if (node->right_ && abmax[dir] >= node->sep_) {
            stack[stacks].node = node->right_;
            stack[stacks].dir = (dir + 1) % 4;
            stacks++;
        }
    }
}

/*************************************** Point3dTree ************************/

Point3dTree::Point3dTree(const Point2d& pmin, const Point2d& pmax)
{
    tree_ = new ADTree3(pmin, pmax);
}

Point3dTree::~Point3dTree()
{
    delete tree_;
}

void Point3dTree::Insert(const Point2d& p, PointIndex pi)
{
    tree_->Insert(p, pi);
}

void Point3dTree::GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<PointIndex>& pis) const
{
    tree_->GetIntersecting(pmin, pmax, pis);
}

/*************************************** Box3dTree **************************/

Box3dTree::Box3dTree(const Point2d& apmin, const Point2d& apmax)
    : box_pmin_{apmin}, box_pmax_{apmax}
{
    tree_ = new ADTree6(apmin, apmax);
}

Box3dTree::~Box3dTree()
{
    delete tree_;
}

void Box3dTree::Insert(const Point2d& bmin, const Point2d& bmax, GenericIndex pi)
{
    tree_->Insert(bmin, bmax, pi);
}

void Box3dTree::GetIntersecting(const Point2d& pmin, const Point2d& pmax, std::vector<GenericIndex>& pis) const
{
    tree_->GetIntersecting(box_pmin_, box_pmax_, pmin, pmax, pis);
}

}  // namespace meshit
