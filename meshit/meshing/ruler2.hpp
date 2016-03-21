#ifndef FILE_NETRULE_HPP
#define FILE_NETRULE_HPP

#include <iostream>
#include <vector>
#include "../general/array.hpp"
#include "../linalg/densemat.hpp"
#include "../gprim/geom2d.hpp"
#include "meshtype.hpp"

namespace meshit {

    class netrule
    {
     private:
        typedef struct tf
        {
            double f1, f2, f3;
        } threefloat;

        class threeint
        {
         public:
            int i1, i2, i3;

            threeint() { }

            threeint(int ai1, int ai2, int ai3)
            {
                i1 = ai1;
                i2 = ai2;
                i3 = ai3;
            }
        };

        int quality;
        char* name;
        std::vector<Point2d> points;
        std::vector<INDEX_2> lines;
        std::vector<Point2d> freezone, freezonelimit;
        std::vector<std::vector<Point2d>*> freezone_i;
        std::vector<Point2d> transfreezone;

        std::vector<int> dellines;
        std::vector<Element2d> elements;
        std::vector<threefloat> tolerances, linetolerances;
        std::vector<threeint> orientations;
        DenseMatrix oldutonewu, oldutofreearea, oldutofreearealimit;
        std::vector<DenseMatrix*> oldutofreearea_i;
        MatrixFixWidth<3> freesetinequ;

        std::vector<Vec2d> linevecs;

        size_t noldp, noldl;
        double fzminx, fzmaxx, fzminy, fzmaxy;

        /// topological distance of line to base element
        std::vector<int> lnearness;

     public:
        netrule();
        ~netrule();

        size_t GetNP() const
        {
            return points.size();
        }

        size_t GetNL() const
        {
            return lines.size();
        }

        size_t GetNE() const
        {
            return elements.size();
        }

        size_t GetNOldP() const
        {
            return noldp;
        }

        size_t GetNOldL() const
        {
            return noldl;
        }

        size_t GetNDelL() const
        {
            return dellines.size();
        }

        size_t GetNOrientations() const
        {
            return orientations.size();
        }

        int GetQuality() const
        {
            return quality;
        }

        int GetLNearness(size_t li) const
        {
            return lnearness[li - 1];
        }

        const Point2d& GetPoint(size_t i) const
        {
            return points[i - 1];
        }

        const INDEX_2& GetLine(size_t i) const
        {
            return lines[i - 1];
        }

        const Element2d& GetElement(size_t i) const
        {
            return elements[i - 1];
        }

        const threeint& GetOrientation(size_t i) const
        {
            return orientations[i - 1];
        }

        int GetDelLine(size_t i) const
        {
            return dellines[i - 1];
        }

        double CalcPointDist(size_t pi, const Point2d& p) const
        {
            double dx = p.X() - points[pi - 1].X();
            double dy = p.Y() - points[pi - 1].Y();
            const threefloat* tfp = &tolerances[pi - 1];
            return tfp->f1 * dx * dx + tfp->f2 * dx * dy + tfp->f3 * dy * dy;
        }

        double CalcLineError(int li, const Vec2d& v) const;

        void SetFreeZoneTransformation(const Vector& u, int tolclass);

        bool IsInFreeZone(const Point2d& p) const
        {
            if (p.X() < fzminx || p.X() > fzmaxx ||
                p.Y() < fzminy || p.Y() > fzmaxy)
                return 0;

            for (size_t i = 0; i < transfreezone.size(); i++) {
                if (freesetinequ(i, 0) * p.X() +
                    freesetinequ(i, 1) * p.Y() +
                    freesetinequ(i, 2) > 0)
                    return 0;
            }
            return 1;
        }

        int IsLineInFreeZone(const Point2d& p1, const Point2d& p2) const
        {
            if ((p1.X() > fzmaxx && p2.X() > fzmaxx) ||
                (p1.X() < fzminx && p2.X() < fzminx) ||
                (p1.Y() > fzmaxy && p2.Y() > fzmaxy) ||
                (p1.Y() < fzminy && p2.Y() < fzminy))
                return 0;
            return IsLineInFreeZone2(p1, p2);
        }

        int IsLineInFreeZone2(const Point2d& p1, const Point2d& p2) const;
        int ConvexFreeZone() const;

        int GetPointNr(size_t ln, int endp) const
        {
            return lines[ln - 1].I(endp);
        }

        const DenseMatrix& GetOldUToNewU() const
        {
            return oldutonewu;
        }

        const char* Name() const
        {
            return name;
        }

        void LoadRule(std::istream& ist);
    };
}  // namespace meshit

#endif

