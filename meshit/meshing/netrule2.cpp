#include "ruler2.hpp"

namespace meshit
{
    netrule::netrule()
    {
        name = new char[1];
        name[0] = static_cast<char>(0);
        quality = 0;
    }

    netrule::~netrule()
    {
        delete[] name;
        for (size_t i = 0; i < oldutofreearea_i.size(); i++) {
            delete oldutofreearea_i[i];
        }
        for (size_t i = 0; i < freezone_i.size(); i++) {
            delete freezone_i[i];
        }
    }

    void netrule::SetFreeZoneTransformation(const Vector& devp, int tolclass)
    {
        double lam1 = 1.0 / tolclass;
        double lam2 = 1. - lam1;

        double mem1[20], mem2[20], mem3[20];

        size_t vs = oldutofreearea.Height();
        FlatVector devfree(vs, mem1);

        size_t fzs = freezone.size();
        transfreezone.resize(fzs);

        if (tolclass <= static_cast<int>(oldutofreearea_i.size())) {
            oldutofreearea_i[tolclass - 1]->Mult(devp, devfree);

            std::vector<Point2d>& fzi = *freezone_i[tolclass - 1];
            for (size_t i = 0; i < fzs; i++) {
                transfreezone[i].X() = fzi[i].X() + devfree[2 * i];
                transfreezone[i].Y() = fzi[i].Y() + devfree[2 * i + 1];
            }
        } else {
            FlatVector devfree1(vs, mem2);
            FlatVector devfree2(vs, mem3);

            oldutofreearea.Mult(devp, devfree1);
            oldutofreearealimit.Mult(devp, devfree2);
            devfree.Set2(lam1, devfree1, lam2, devfree2);

            for (size_t i = 0; i < fzs; i++) {
                transfreezone[i].X() = lam1 * freezone[i].X() + lam2 * freezonelimit[i].X() + devfree[2 * i];
                transfreezone[i].Y() = lam1 * freezone[i].Y() + lam2 * freezonelimit[i].Y() + devfree[2 * i + 1];
            }
        }

        if (fzs > 0) {
            fzmaxx = fzminx = transfreezone[0].X();
            fzmaxy = fzminy = transfreezone[0].Y();
        }

        for (size_t i = 1; i < fzs; i++) {
            if (transfreezone[i].X() > fzmaxx) fzmaxx = transfreezone[i].X();
            if (transfreezone[i].X() < fzminx) fzminx = transfreezone[i].X();
            if (transfreezone[i].Y() > fzmaxy) fzmaxy = transfreezone[i].Y();
            if (transfreezone[i].Y() < fzminy) fzminy = transfreezone[i].Y();
        }

        for (size_t i = 0; i < fzs; i++) {
            Point2d p1 = transfreezone[i];
            Point2d p2 = transfreezone[(i + 1) % fzs];

            Vec2d vn(p2.Y() - p1.Y(), p1.X() - p2.X());

            double len2 = vn.Length2();

            if (len2 < 1e-10) {
                freesetinequ(i, 0) = 0;
                freesetinequ(i, 1) = 0;
                freesetinequ(i, 2) = -1;
            } else {
                vn /= sqrt(len2);  // scaling necessary ?

                freesetinequ(i, 0) = vn.X();
                freesetinequ(i, 1) = vn.Y();
                freesetinequ(i, 2) = -(p1.X() * vn.X() + p1.Y() * vn.Y());
            }
        }
    }

    bool netrule::IsLineInFreeZone(const Point2d& p1, const Point2d& p2) const
    {
        if ((p1.X() > fzmaxx && p2.X() > fzmaxx) ||
            (p1.X() < fzminx && p2.X() < fzminx) ||
            (p1.Y() > fzmaxy && p2.Y() > fzmaxy) ||
            (p1.Y() < fzminy && p2.Y() < fzminy)) {
            return false;
        }

        for (size_t i = 0; i < transfreezone.size(); i++) {
            if (freesetinequ(i, 0) * p1.X() + freesetinequ(i, 1) * p1.Y() + freesetinequ(i, 2) > -1e-8 &&
                freesetinequ(i, 0) * p2.X() + freesetinequ(i, 1) * p2.Y() + freesetinequ(i, 2) > -1e-8) {
                return false;
            }
        }

        double nx = (p2.Y() - p1.Y());
        double ny = -(p2.X() - p1.X());
        double nl = sqrt(nx * nx + ny * ny);

        if (nl > 1e-8) {
            nx /= nl;
            ny /= nl;
            double c = -(p1.X() * nx + p1.Y() * ny);

            bool all_left = true;
            bool all_right = true;

            for (size_t i = 0; i < transfreezone.size(); i++) {
                bool left = transfreezone[i].X() * nx + transfreezone[i].Y() + c < 1e-7;
                bool right = transfreezone[i].X() * nx + transfreezone[i].Y() + c > -1e-7;

                if (!left) all_left = false;
                if (!right) all_right = false;
            }
            if (all_left || all_right) return false;
        }

        return true;
    }

    int netrule::ConvexFreeZone() const
    {
        size_t n = transfreezone.size();
        for (size_t i = 0; i < n; i++) {
            const bool counterclockwise = CCW(transfreezone[i],
                                              transfreezone[(i + 1) % n],
                                              transfreezone[(i + 2) % n], 1e-7);
            if (!counterclockwise)
                return 0;
        }
        return 1;
    }

    double netrule::CalcLineError(size_t li, const Vec2d& v) const
    {
        double dx = v.X() - linevecs[li].X();
        double dy = v.Y() - linevecs[li].Y();

        const threefloat* ltf = &linetolerances[li];
        return ltf->f1 * dx * dx +
               ltf->f2 * dx * dy +
               ltf->f3 * dy * dy;
    }
}  // namespace meshit
