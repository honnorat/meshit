#include <cstring>
#include "adtree.hpp"

namespace meshit {

    /* ******************************* ADTree3 ******************************* */

    ADTreeNode3::ADTreeNode3()
    {
        pi = -1;

        left = NULL;
        right = NULL;
        father = NULL;
        nchilds = 0;
    }

    void ADTreeNode3::DeleteChilds()
    {
        if (left) {
            left->DeleteChilds();
            delete left;
            left = NULL;
        }
        if (right) {
            right->DeleteChilds();
            delete right;
            right = NULL;
        }
    }

    BlockAllocator ADTreeNode3::ball(sizeof(ADTreeNode3));

    void* ADTreeNode3::operator new(size_t s)
    {
        return ball.Alloc();
    }

    void ADTreeNode3::operator delete(void* p)
    {
        ball.Free(p);
    }

    ADTree3::ADTree3(const double* acmin,
                     const double* acmax)
            : ela(0)
    {
        memcpy(cmin, acmin, 3 * sizeof(double));
        memcpy(cmax, acmax, 3 * sizeof(double));

        root = new ADTreeNode3;
        root->sep = (cmin[0] + cmax[0]) / 2;
    }

    ADTree3::~ADTree3()
    {
        root->DeleteChilds();
        delete root;
    }

    void ADTree3::Insert(const Point<3>& p, int pi)
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
                node->data[2] = p[2];
                node->pi = pi;

                if (ela.size() < static_cast<size_t>(pi + 1)) {
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

            if (++dir == 3) {
                dir = 0;
            }
        }

        next = new ADTreeNode3;
        next->data[0] = p[0];
        next->data[1] = p[1];
        next->data[2] = p[2];
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

    void ADTree3::DeleteElement(int pi)
    {
        ADTreeNode3* node = ela[pi];

        node->pi = -1;

        node = node->father;
        while (node) {
            node->nchilds--;
            node = node->father;
        }
    }

    void ADTree3::GetIntersecting(const double* bmin,
                                  const double* bmax,
                                  std::vector<size_t>& pis) const
    {
        static Array<ADTreeNode3*> stack(1000);
        static Array<size_t> stackdir(1000);
        ADTreeNode3* node;

        stack.resize(1000);
        stackdir.resize(1000);
        pis.resize(0);

        stack[0] = root;
        stackdir[0] = 0;

        size_t stacks = 1;

        while (stacks) {
            stacks--;
            node = stack[stacks];
            size_t dir = stackdir[stacks];

            if (node->pi != -1) {
                if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
                    node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
                    node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

                    pis.push_back(node->pi);
            }

            int ndir = dir + 1;
            if (ndir == 3)
                ndir = 0;

            if (node->left && bmin[dir] <= node->sep) {
                stack[stacks] = node->left;
                stackdir[stacks] = ndir;
                stacks++;
            }
            if (node->right && bmax[dir] >= node->sep) {
                stack[stacks] = node->right;
                stackdir[stacks] = ndir;
                stacks++;
            }
        }
    }

    void ADTree3::PrintRec(std::ostream& ost, const ADTreeNode3* node) const
    {
        if (node->data) {
            ost << node->pi << ": ";
            ost << node->nchilds << " childs, ";
            for (int i = 0; i < 3; i++)
                ost << node->data[i] << " ";
            ost << std::endl;
        }
        if (node->left)
            PrintRec(ost, node->left);
        if (node->right)
            PrintRec(ost, node->right);
    }

    /* ******************************* ADTree6 ******************************* */

    ADTreeNode6::ADTreeNode6()
    {
        pi = -1;

        left = NULL;
        right = NULL;
        father = NULL;
        nchilds = 0;
    }

    void ADTreeNode6::DeleteChilds()
    {
        if (left) {
            left->DeleteChilds();
            delete left;
            left = NULL;
        }
        if (right) {
            right->DeleteChilds();
            delete right;
            right = NULL;
        }
    }

    BlockAllocator ADTreeNode6::ball(sizeof(ADTreeNode6));

    void* ADTreeNode6::operator new(size_t s)
    {
        return ball.Alloc();
    }

    void ADTreeNode6::operator delete(void* p)
    {
        ball.Free(p);
    }

    ADTree6::ADTree6(const double* acmin,
                     const double* acmax)
            : ela(0)
    {
        memcpy(cmin, acmin, 6 * sizeof(double));
        memcpy(cmax, acmax, 6 * sizeof(double));

        root = new ADTreeNode6;
        root->sep = (cmin[0] + cmax[0]) / 2;
    }

    ADTree6::~ADTree6()
    {
        root->DeleteChilds();
        delete root;
    }

    void ADTree6::Insert(const Point3d& bmin_, const Point3d& bmax_, int pi)
    {
        double p[6];

        for (int i = 0; i < 3; i++) {
            p[i] = bmin_[i];
            p[i + 3] = bmax_[i];
        }

        ADTreeNode6* node = nullptr;
        ADTreeNode6* next = root;
        int dir = 0;
        bool lr = false;

        double bmin[6];
        double bmax[6];

        memcpy(bmin, cmin, 6 * sizeof(double));
        memcpy(bmax, cmax, 6 * sizeof(double));

        while (next) {
            node = next;

            if (node->pi == -1) {
                memcpy(node->data, p, 6 * sizeof(double));
                node->pi = pi;

                if (ela.size() < static_cast<size_t>(pi + 1)) {
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

            if (++dir == 6) dir = 0;
        }

        next = new ADTreeNode6;
        memcpy(next->data, p, 6 * sizeof(double));
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

    void ADTree6::DeleteElement(int pi)
    {
        ADTreeNode6* node = ela[pi];

        node->pi = -1;

        node = node->father;
        while (node) {
            node->nchilds--;
            node = node->father;
        }
    }

    void ADTree6::PrintMemInfo(std::ostream& ost) const
    {
        ost << Elements() << " elements a " << sizeof(ADTreeNode6)
        << " Bytes = " << Elements() * sizeof(ADTreeNode6) << std::endl;
        ost << "maxind = " << ela.size() << " = " << sizeof(ADTreeNode6*) * ela.size() << " Bytes" << std::endl;
    }

    class inttn6
    {
     public:
        int dir;
        ADTreeNode6* node;
    };

    void ADTree6::GetIntersecting(const double* bmin,
                                  const double* bmax,
                                  std::vector<size_t>& pis) const
    {
        pis.resize(0);

        std::vector<inttn6> stack(100);
        stack[0].node = root;
        stack[0].dir = 0;
        size_t stacks = 1;

        while (stacks) {
            stacks--;

            ADTreeNode6* node = stack[stacks].node;
            int dir = stack[stacks].dir;

            if (node->pi != -1) {
                if (node->data[0] > bmax[0] ||
                    node->data[1] > bmax[1] ||
                    node->data[2] > bmax[2] ||
                    node->data[3] < bmin[3] ||
                    node->data[4] < bmin[4] ||
                    node->data[5] < bmin[5]) {
                    // nothing
                } else {
                    pis.push_back(node->pi);
                }
            }

            int ndir = (dir + 1) % 6;

            if (node->left && bmin[dir] <= node->sep) {
                stack[stacks].node = node->left;
                stack[stacks].dir = ndir;
                stacks++;
            }
            if (node->right && bmax[dir] >= node->sep) {
                stack[stacks].node = node->right;
                stack[stacks].dir = ndir;
                stacks++;
            }
        }
    }

    void ADTree6::PrintRec(std::ostream& ost, const ADTreeNode6* node) const
    {
        if (node->data) {
            ost << node->pi << ": ";
            ost << node->nchilds << " childs, ";
            for (int i = 0; i < 6; i++)
                ost << node->data[i] << " ";
            ost << std::endl;
        }
        if (node->left)
            PrintRec(ost, node->left);
        if (node->right)
            PrintRec(ost, node->right);
    }

    int ADTree6::DepthRec(const ADTreeNode6* node) const
    {
        int ldepth = 0;
        int rdepth = 0;

        if (node->left)
            ldepth = DepthRec(node->left);
        if (node->right)
            rdepth = DepthRec(node->right);
        return 1 + std::max(ldepth, rdepth);
    }

    int ADTree6::ElementsRec(const ADTreeNode6* node) const
    {
        int els = 1;
        if (node->left)
            els += ElementsRec(node->left);
        if (node->right)
            els += ElementsRec(node->right);
        return els;
    }

    /* ************************************* Point3dTree ********************** */


    Point3dTree::Point3dTree(const Point<3>& pmin, const Point<3>& pmax)
    {
        double pmi[3], pma[3];
        for (int i = 0; i < 3; i++) {
            pmi[i] = pmin(i);
            pma[i] = pmax(i);
        }
        tree = new ADTree3(pmi, pma);
    }

    Point3dTree::~Point3dTree()
    {
        delete tree;
    }

    void Point3dTree::Insert(const Point<3>& p, int pi)
    {
        tree->Insert(p, pi);
    }

    void Point3dTree::GetIntersecting(const Point<3>& pmin, const Point<3>& pmax, std::vector<size_t>& pis) const
    {
        double pmi[3], pma[3];
        for (int i = 0; i < 3; i++) {
            pmi[i] = pmin(i);
            pma[i] = pmax(i);
        }
        tree->GetIntersecting(pmi, pma, pis);
    }

    Box3dTree::Box3dTree(const Box<3>& abox)
            : boxpmin(abox.PMin()),
              boxpmax(abox.PMax())
    {
        double tpmin[6], tpmax[6];
        for (size_t i = 0; i < 3; i++) {
            tpmin[i] = tpmin[i + 3] = boxpmin[i];
            tpmax[i] = tpmax[i + 3] = boxpmax[i];
        }
        tree = new ADTree6(tpmin, tpmax);
    }

    Box3dTree::Box3dTree(const Point3d& apmin, const Point3d& apmax)
            : boxpmin(apmin),
              boxpmax(apmax)
    {
        double tpmin[6], tpmax[6];
        for (int i = 0; i < 3; i++) {
            tpmin[i] = tpmin[i + 3] = boxpmin[i];
            tpmax[i] = tpmax[i + 3] = boxpmax[i];
        }
        tree = new ADTree6(tpmin, tpmax);
    }

    Box3dTree::~Box3dTree()
    {
        delete tree;
    }

    void Box3dTree::Insert(const Point3d& bmin, const Point3d& bmax, int pi)
    {
        tree->Insert(bmin, bmax, pi);
    }

    void Box3dTree::GetIntersecting(const Point3d& pmin, const Point3d& pmax, std::vector<size_t>& pis) const
    {
        double tpmin[6];
        double tpmax[6];

        for (size_t i = 0; i < 3; i++) {
            tpmin[i] = boxpmin[i];
            tpmax[i] = pmax[i];

            tpmin[i + 3] = pmin[i];
            tpmax[i + 3] = boxpmax[i];
        }

        tree->GetIntersecting(tpmin, tpmax, pis);
    }
}  // namespace meshit
